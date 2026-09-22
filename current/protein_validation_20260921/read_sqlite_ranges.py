"""Read SQLite table records from a public immutable HTTP range source.

Only table b-trees are supported. Byte ranges are validated and cached exactly.
SQLite format reference: https://www.sqlite.org/fileformat.html
No source database is changed. Use --self-test for comparison with sqlite3.
"""
import argparse
import hashlib
import json
from pathlib import Path
import sqlite3
import struct
import tempfile
import urllib.request


def varint(b, at):
    n = 0
    for i in range(9):
        v = b[at + i]
        if i == 8:
            return (n << 8) | v, at + 9
        n = (n << 7) | (v & 127)
        if v < 128:
            return n, at + i + 1
    raise AssertionError


def record(payload, encoding='utf-8'):
    hs, at = varint(payload, 0)
    serials = []
    while at < hs:
        v, at = varint(payload, at)
        serials.append(v)
    assert at == hs
    out = []
    at = hs
    for v in serials:
        if v == 0:
            out.append(None)
            continue
        if v in (8, 9):
            out.append(v - 8)
            continue
        if v in (10, 11):
            raise ValueError('Reserved SQLite serial type')
        if v <= 6:
            length = [0, 1, 2, 3, 4, 6, 8][v]
            value = int.from_bytes(payload[at:at+length], 'big', signed=True)
        elif v == 7:
            length = 8
            value = struct.unpack('>d', payload[at:at+length])[0]
        else:
            length = (v - (13 if v % 2 else 12)) // 2
            value = payload[at:at+length]
            if v % 2:
                value = value.decode(encoding)
        assert len(payload[at:at+length]) == length
        out.append(value)
        at += length
    assert at == len(payload)
    return out


class Reader:
    def __init__(self, source, cache, block=65536, seed=None):
        self.source = source
        self.local = not source.startswith('https://')
        self.cache = Path(cache)
        self.cache.mkdir(parents=True, exist_ok=True)
        self.block = block
        self.downloaded = 0
        self.size = Path(source).stat().st_size if self.local else None
        if seed and not (self.cache/'000000000000.bin').exists():
            (self.cache/'000000000000.bin').write_bytes(Path(seed).read_bytes())
        h = self.read(0, 100)
        assert h[:16] == b'SQLite format 3\0'
        self.pagesize = int.from_bytes(h[16:18], 'big')
        if self.pagesize == 1:
            self.pagesize = 65536
        self.usable = self.pagesize - h[20]
        self.encoding = {1:'utf-8', 2:'utf-16le', 3:'utf-16be'}[int.from_bytes(h[56:60], 'big')]
        self.size = int.from_bytes(h[28:32], 'big') * self.pagesize
        self.schema = list(self.rows(1))

    def read(self, start, count):
        if self.local:
            with open(self.source, 'rb') as f:
                f.seek(start)
                b = f.read(count)
            assert len(b) == count
            return b
        ans = bytearray()
        end = start + count
        while start < end:
            offset = start // self.block * self.block
            p = self.cache/f'{offset:012d}.bin'
            if not p.exists():
                last = offset + self.block - 1
                if self.size is not None:
                    last = min(last, self.size-1)
                req = urllib.request.Request(self.source, headers={'Range': f'bytes={offset}-{last}'})
                with urllib.request.urlopen(req, timeout=40) as r:
                    assert r.status == 206, 'Server ignored byte range'
                    cr = r.headers['Content-Range']
                    assert cr.startswith(f'bytes {offset}-{last}/'), cr
                    b = r.read()
                    assert len(b) == last-offset+1
                    p.write_bytes(b)
                    self.downloaded += len(b)
            else:
                b = p.read_bytes()
            take = min(end-start, len(b)-(start-offset))
            assert take > 0
            ans.extend(b[start-offset:start-offset+take])
            start += take
        return bytes(ans)

    def page(self, page):
        return self.read((page-1)*self.pagesize, self.pagesize)

    def rows(self, root, with_rowid=False):
        b = self.page(root)
        h = 100 if root == 1 else 0
        kind = b[h]
        n = int.from_bytes(b[h+3:h+5], 'big')
        assert kind in (5, 13), ('Not a table b-tree', root, kind)
        head = 12 if kind == 5 else 8
        ptrs = [int.from_bytes(b[h+head+2*i:h+head+2*i+2], 'big') for i in range(n)]
        if kind == 5:
            for p in ptrs:
                yield from self.rows(int.from_bytes(b[p:p+4], 'big'), with_rowid)
            yield from self.rows(int.from_bytes(b[h+8:h+12], 'big'), with_rowid)
            return
        for p in ptrs:
            size, at = varint(b, p)
            rowid, at = varint(b, at)
            if rowid >= 2**63:
                rowid -= 2**64
            if size <= self.usable - 35:
                local = size
            else:
                minimum = ((self.usable-12)*32)//255-23
                k = minimum + (size-minimum) % (self.usable-4)
                local = k if k <= self.usable-35 else minimum
            payload = bytearray(b[at:at+local])
            left = size-local
            overflow = int.from_bytes(b[at+local:at+local+4], 'big') if left else 0
            visited = set()
            while left:
                assert overflow and overflow not in visited
                visited.add(overflow)
                ov = self.page(overflow)
                take = min(left, self.usable-4)
                payload.extend(ov[4:4+take])
                overflow = int.from_bytes(ov[:4], 'big')
                left -= take
            assert len(payload) == size
            val = record(bytes(payload), self.encoding)
            yield (rowid, val) if with_rowid else val

    def table(self, name, with_rowid=False):
        matches = [r for r in self.schema if r[0]=='table' and r[1]==name]
        assert len(matches) == 1, name
        return self.rows(matches[0][3], with_rowid)

    def manifest(self):
        return dict(source=self.source, database_size=self.size, page_size=self.pagesize,
                    ranges=[dict(start=int(p.stem),bytes=p.stat().st_size,
                                 sha256=hashlib.sha256(p.read_bytes()).hexdigest())
                            for p in sorted(self.cache.glob('*.bin'))])


def self_test():
    with tempfile.TemporaryDirectory() as t:
        p = Path(t)/'example.db'
        c = sqlite3.connect(p)
        c.execute('pragma page_size=1024')
        c.execute('create table demo(k integer, text_value text, blob_value blob, real_value real, null_value text)')
        inputs=[(i-200, '汉字'+str(i)+'x'*(i%7*700), bytes(range(256))*(i%9), i*.0123, None) for i in range(1000)]
        c.executemany('insert into demo values(?,?,?,?,?)', inputs)
        c.commit()
        expected = c.execute('select * from demo').fetchall()
        c.close()
        r = Reader(str(p), Path(t)/'cache')
        observed = [tuple(x) for x in r.table('demo')]
        assert observed == expected
        print('SQLite extraction checked against sqlite3 for 1000 records, Unicode, numbers, blobs and overflow pages')


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--self-test',action='store_true')
    ap.add_argument('--url')
    ap.add_argument('--cache',type=Path)
    ap.add_argument('--seed',type=Path)
    ap.add_argument('--schema-out',type=Path)
    args = ap.parse_args()
    if args.self_test:
        self_test()
    else:
        rr = Reader(args.url,args.cache,seed=args.seed)
        args.schema_out.write_text(json.dumps(rr.schema,ensure_ascii=False,indent=2))
        print(json.dumps({'tables':len(rr.schema),'downloaded_bytes':rr.downloaded,'schema':str(args.schema_out)}))
