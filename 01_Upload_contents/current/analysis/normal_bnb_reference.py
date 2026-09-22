"""Extract anatomical expression context, not disease validation, from GSE107574."""
import argparse
import gzip
import hashlib
import io
import json
from pathlib import Path
import pandas as pd
from nerve_programs import PANELS

URL='https://ftp.ncbi.nlm.nih.gov/geo/series/GSE107nnn/GSE107574/suppl/GSE107574_fpkm_gene_expression_values.xlsx.gz'


def run(path,out):
    data=pd.read_excel(io.BytesIO(gzip.decompress(path.read_bytes())))
    data['gene']=data.gene_id.str.rsplit('_',n=1).str[-1]
    rows=[]
    for program in list(PANELS)[:3]:
        for gene in PANELS[program]:
            d=data[data.gene.eq(gene)]
            assert len(d)==1
            d=d.iloc[0]
            for sample in ['P3pHEnd_EC','P8_pHEnd_EC_basal','32P1','203P1','346P1','347P1']:
                rows.append(dict(program=program,gene=gene,original_gene_id=d.gene_id,sample=sample,
                    preparation='Cultured endothelial cells' if sample.startswith('P') else 'Laser captured normal microvessels',
                    fpkm=float(d[sample]),source_url=URL))
    result=pd.DataFrame(rows)
    result.to_csv(out/'BNB_normal_reference_expression.csv',index=False)
    summary=result.query('preparation=="Laser captured normal microvessels"').groupby(['program','gene']).agg(
        n_samples=('sample','size'),samples_fpkm_above_zero=('fpkm',lambda x:int((x>0).sum())),
        mean_fpkm=('fpkm','mean'),min_fpkm=('fpkm','min'),max_fpkm=('fpkm','max')).reset_index()
    summary.to_csv(out/'BNB_normal_reference_summary.csv',index=False)
    (out/'BNB_reference_provenance.json').write_text(json.dumps(dict(source_url=URL,
        sha256=hashlib.sha256(path.read_bytes()).hexdigest(),bytes=path.stat().st_size,
        normal_nerve_microvessel_samples=4,culture_conditions=2,target_genes=11,
        role='Anatomical expression context only. No CIDP comparison or disease replication.'),indent=2))
    print(summary.to_string(index=False))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--input',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
    a=p.parse_args();run(a.input,a.out)
