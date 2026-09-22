"""Inspect published instrument-header annotations for specimen mapping."""
from pathlib import Path
import concurrent.futures
import csv
import hashlib
import json
import re
import urllib.request

ROOT=Path(__file__).resolve().parent
BASE='https://ftp.pride.ebi.ac.uk/pride/data/archive/2025/09/PXD056286/'
OUT=ROOT/'access/raw_headers'
OUT.mkdir(parents=True,exist_ok=True)
rows=list(csv.DictReader((ROOT/'access/source_run_mapping.csv').open()))

def one(row):
    name=row['run_identifier']+'.raw'
    file=OUT/(name+'.head')
    url=BASE+name
    if not file.exists():
        req=urllib.request.Request(url,headers={'Range':'bytes=0-65535'})
        with urllib.request.urlopen(req,timeout=40) as r:
            assert r.status==206 and r.headers['Content-Range'].startswith('bytes 0-65535/')
            data=r.read()
        assert len(data)==65536
        file.write_bytes(data)
    data=file.read_bytes()
    found=[]
    for off in [0,1]:
        s=data[off:len(data)-(len(data)-off)%2].decode('utf-16le',errors='replace')
        found.extend(re.findall(r'[\x20-\x7e]{5,}',s))
    found.extend(x.decode() for x in re.findall(rb'[\x20-\x7e]{5,}',data))
    found=list(dict.fromkeys(found))
    (OUT/(name+'.strings.txt')).write_text('\n'.join(found))
    keys=[s for s in found if re.search(r'CIDP|control|patient|sample|specimen|VNP\d|N\d{2}',s,re.I)]
    return {'run_identifier':row['run_identifier'],'source':url,'bytes':len(data),'sha256':hashlib.sha256(data).hexdigest(),'relevant_strings':keys}

def main():
    results=[]
    with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
        futures={pool.submit(one,r):r['run_identifier'] for r in rows}
        for future in concurrent.futures.as_completed(futures):
            name=futures[future]
            try:result=future.result()
            except Exception as e:result={'run_identifier':name,'error':str(e)}
            results.append(result)
            print(json.dumps(result),flush=True)
    results.sort(key=lambda x:x['run_identifier'])
    (ROOT/'access/raw_header_review.json').write_text(json.dumps(results,indent=2))

if __name__=='__main__':main()
