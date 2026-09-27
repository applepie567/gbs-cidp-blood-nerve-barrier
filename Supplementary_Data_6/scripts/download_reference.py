from pathlib import Path
import concurrent.futures, urllib.request, hashlib, json, datetime
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]
prior=ROOT/'input/prior/External_nerve_20260926'
rows=pd.read_csv(prior/'results/GSE285983_macrophage_Fc_counts.csv').to_dict('records')
dest=ROOT/'input/nerve_h5'
dest.mkdir(exist_ok=True)

def get(row):
    name=row['h5_filename']; p=dest/name
    gsm=name.split('_')[0]
    url=f'https://ftp.ncbi.nlm.nih.gov/geo/samples/{gsm[:-3]}nnn/{gsm}/suppl/{name}'
    if not p.exists():
        tmp=p.with_suffix('.h5.partial')
        with urllib.request.urlopen(url,timeout=60) as response,tmp.open('wb') as target:
            while chunk:=response.read(1024*1024): target.write(chunk)
        tmp.replace(p)
    with p.open('rb') as f: sha=hashlib.file_digest(f,'sha256').hexdigest()
    assert sha==row['h5_sha256'],name
    print(name,p.stat().st_size,flush=True)
    return dict(file=str(p.relative_to(ROOT)),url=url,bytes=p.stat().st_size,sha256=sha,retrieved_date='2026-09-27')

if __name__=='__main__':
    metadata=ROOT/'input/GSE285983_metadata_all.csv.gz'
    metadata_url='https://ftp.ncbi.nlm.nih.gov/geo/series/GSE285nnn/GSE285983/suppl/GSE285983_metadata_all.csv.gz'
    if not metadata.exists():
        with urllib.request.urlopen(metadata_url,timeout=60) as response,metadata.open('wb') as target:
            while chunk:=response.read(1024*1024): target.write(chunk)
    assert hashlib.sha256(metadata.read_bytes()).hexdigest()=='53b395600269979c9689fd63ae06fb95bc6228138e7cfe632b676faca133ec8a'
    with concurrent.futures.ThreadPoolExecutor(max_workers=4) as ex:
        results=list(ex.map(get,rows))
    p=ROOT/'input/GSE285983_metadata_all.csv.gz'
    results.append(dict(file=str(p.relative_to(ROOT)),url='https://ftp.ncbi.nlm.nih.gov/geo/series/GSE285nnn/GSE285983/suppl/GSE285983_metadata_all.csv.gz',bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest(),retrieved_date='2026-09-27'))
    pd.DataFrame(results).to_csv(ROOT/'source_manifest.csv',index=False)
    print('All reference inputs verified',len(results),flush=True)
