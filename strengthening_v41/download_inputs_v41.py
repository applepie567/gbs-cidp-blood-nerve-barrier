#!/usr/bin/env python3
"""Download verified public input files using the retained manifest.

This script performs no uploads. Usage:
python download_inputs_v41.py --manifest results/new_input_manifest.csv --out local_inputs
"""
import argparse,csv,hashlib
from pathlib import Path
from urllib.request import Request,urlopen
from concurrent.futures import ThreadPoolExecutor

def download(row,out):
    relative=Path(row['file'])
    if relative.is_absolute() or '..' in relative.parts:raise ValueError('Unsafe relative input path')
    dest=out/relative;dest.parent.mkdir(parents=True,exist_ok=True)
    expected=row['sha256']
    if dest.exists() and hashlib.sha256(dest.read_bytes()).hexdigest()==expected:return str(relative)+' verified'
    temp=dest.with_name(dest.name+'.partial');sha=hashlib.sha256();size=0
    with urlopen(Request(row['url'],headers={'User-Agent':'ResearchDataAudit/1.0'}),timeout=120) as response,temp.open('wb') as target:
        while True:
            data=response.read(1024*1024)
            if not data:break
            target.write(data);sha.update(data);size+=len(data)
    if size!=int(row['bytes']) or sha.hexdigest()!=expected:
        temp.unlink(missing_ok=True);raise ValueError(f'Checksum or size mismatch for {relative}')
    temp.replace(dest);return str(relative)+' downloaded and verified'

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--manifest',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args()
    with a.manifest.open() as f:rows=list(csv.DictReader(f))
    rows=[r for r in rows if r['file'].startswith(('nerve_h5/','spatial/')) or r['file'] in ['metadata_all.csv.gz','metadata_ic.csv.gz']]
    with ThreadPoolExecutor(max_workers=4) as executor:
        for message in executor.map(lambda row:download(row,a.out),rows):print(message,flush=True)
