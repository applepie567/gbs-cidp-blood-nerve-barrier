"""Thread-safe, verified range cache for the immutable public SQLite source."""
from pathlib import Path
import concurrent.futures
import os
import threading
import time
import urllib.error
import urllib.request

from read_sqlite_ranges import Reader

class CachedReader(Reader):
    def __init__(self,source,cache,batch_bytes=1048576):
        self.batch_bytes=batch_bytes
        self._guard=threading.Lock()
        self._locks={}
        super().__init__(source,cache,block=65536)

    def _fetch(self,offset):
        large=offset//self.batch_bytes*self.batch_bytes
        with self._guard:lock=self._locks.setdefault(large,threading.Lock())
        with lock:
            if (self.cache/f'{offset:012d}.bin').exists():return
            last=large+self.batch_bytes-1
            if self.size is not None:last=min(last,self.size-1)
            request=urllib.request.Request(self.source,headers={'Range':f'bytes={large}-{last}'})
            for attempt in range(3):
                try:
                    with urllib.request.urlopen(request,timeout=40) as result:
                        assert result.status==206
                        assert result.headers['Content-Range'].startswith(f'bytes {large}-{last}/')
                        data=result.read()
                    assert len(data)==last-large+1
                    break
                except (TimeoutError,urllib.error.URLError):
                    if attempt==2:raise
                    time.sleep(1)
            for i in range(0,len(data),self.block):
                target=self.cache/f'{large+i:012d}.bin'
                block=data[i:i+self.block]
                if target.exists():assert target.read_bytes()==block
                else:
                    temp=target.with_suffix(f'.tmp-{os.getpid()}-{threading.get_ident()}')
                    temp.write_bytes(block)
                    temp.replace(target)
            with self._guard:self.downloaded+=len(data)

    def read(self,start,count):
        if self.local:return super().read(start,count)
        out=bytearray();end=start+count
        while start<end:
            offset=start//self.block*self.block
            path=self.cache/f'{offset:012d}.bin'
            if not path.exists():self._fetch(offset)
            data=path.read_bytes()
            take=min(end-start,len(data)-(start-offset))
            assert take>0
            out.extend(data[start-offset:start-offset+take]);start+=take
        return bytes(out)

    def parallel_table(self,name,workers=6):
        root=next(r[3] for r in self.schema if r[0]=='table' and r[1]==name)
        b=self.page(root)
        if b[0]!=5:
            yield from self.table(name)
            return
        n=int.from_bytes(b[3:5],'big')
        pointers=[int.from_bytes(b[12+2*i:14+2*i],'big') for i in range(n)]
        children=[int.from_bytes(b[p:p+4],'big') for p in pointers]+[int.from_bytes(b[8:12],'big')]
        with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as pool:
            futures=[pool.submit(lambda child:list(self.rows(child)),child) for child in children]
            for i,future in enumerate(futures):
                rows=future.result()
                print(f'{name} subtree {i+1}/{len(children)}, {len(rows)} records',flush=True)
                yield from rows
