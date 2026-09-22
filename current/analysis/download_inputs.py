"""Download public extension inputs and verify their recorded SHA256 values."""
import argparse
import csv
import hashlib
import shutil
import urllib.request
from pathlib import Path

def digest(path):
    with path.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()

def run(manifest,out):
    with manifest.open(newline='') as f:rows=list(csv.DictReader(f))
    for row in rows:
        dest=out/row['file'];dest.parent.mkdir(parents=True,exist_ok=True)
        if dest.exists() and digest(dest)==row['sha256']:
            print('Verified existing',row['file']);continue
        tmp=dest.with_suffix(dest.suffix+'.part')
        request=urllib.request.Request(row['url'],headers={'User-Agent':'GBS-CIDP-reproducibility/20260921'})
        with urllib.request.urlopen(request,timeout=120) as response,tmp.open('wb') as f:
            shutil.copyfileobj(response,f)
        if digest(tmp)!=row['sha256']:
            tmp.unlink();raise ValueError('Checksum mismatch for '+row['file'])
        if tmp.stat().st_size!=int(row['bytes']):
            tmp.unlink();raise ValueError('Size mismatch for '+row['file'])
        tmp.replace(dest);print('Downloaded and verified',row['file'])

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--manifest',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
    a=p.parse_args();run(a.manifest,a.out)
