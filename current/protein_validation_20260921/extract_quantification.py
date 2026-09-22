"""Export stored PXD056286 protein intensities without assigning diagnoses.

Abundance arrays contain little-endian doubles followed by a presence byte.
The source schema and dynamic-column metadata must establish vector identities
before slot labels are interpreted as instrument runs. Missing slots are blank.
This script performs neither patient aggregation nor a disease comparison.
"""
from pathlib import Path
import base64
import csv
import gzip
import json
import math
import re
import sqlite3
import struct
import sys
from collections import Counter

ROOT=Path(__file__).resolve().parent
from public_sqlite import CachedReader

SOURCE='https://ftp.pride.ebi.ac.uk/pride/data/archive/2025/09/PXD056286/QExactiveHF02_18988_N1530.pdResult'
ACCESS=ROOT/'access'
OUT=ROOT/'results'
OUT.mkdir(exist_ok=True)
FC=['FCGR1A','FCGR2A','FCGR2B','FCGR3A','FCGR3B','FCGRT','FCER1G','TYROBP','SYK','LYN','HCK']

def get_columns(sql):
    con=sqlite3.connect(':memory:')
    con.execute(sql)
    table=con.execute("select name from sqlite_master where type='table'").fetchone()[0]
    names=[r[1] for r in con.execute('pragma table_info("'+table+'")')]
    con.close()
    return names

def decode(blob,kind='d',expected=None):
    if blob is None:return None
    width=struct.calcsize('<'+kind+'B')
    assert len(blob)%width==0
    pairs=list(struct.iter_unpack('<'+kind+'B',blob))
    assert all(flag in (0,1) for _,flag in pairs)
    if expected is not None:assert len(pairs)==expected
    assert all(math.isfinite(value) for value,flag in pairs if flag)
    return [value if flag else None for value,flag in pairs]

def raw_encode(obj):
    if isinstance(obj,bytes):return {'base64':base64.b64encode(obj).decode()}
    raise TypeError(type(obj))

def run():
    reader=CachedReader(SOURCE,ACCESS/'pdResult_blocks64')
    schema={r[1]:r[4] for r in reader.schema if r[0]=='table'}
    cols=get_columns(schema['TargetProteins'])
    meta=['UniqueSequenceID','Accession','GeneSymbol','Description','IsMasterProtein',
          'ExcludedBy','UniquePeptidesCount','GroupUniquePeptidesCount','PeptideGroupCount',
          'ProteinGroupIDs','ProteinFDRConfidence','Expqvalue']
    slots=['slot_'+str(i).zfill(2) for i in range(1,20)]
    counts=Counter();statuses=Counter();fc_hits=[]
    factors=[[] for _ in slots]
    with (OUT/'pxd056286_protein_inventory.csv').open('w',newline='') as fi, \
         (OUT/'pxd056286_abundance_raw.csv').open('w',newline='') as fr, \
         (OUT/'pxd056286_abundance_normalized.csv').open('w',newline='') as fn, \
         gzip.open(OUT/'pxd056286_protein_records.jsonl.gz','wt') as fg:
        iw=csv.writer(fi);rw=csv.writer(fr);nw=csv.writer(fn)
        iw.writerow(meta+['has_abundance_vector','valid_raw_slots','valid_normalized_slots'])
        rw.writerow(['UniqueSequenceID','Accession','GeneSymbol']+slots)
        nw.writerow(['UniqueSequenceID','Accession','GeneSymbol']+slots)
        for i,row in enumerate(reader.parallel_table('TargetProteins')):
            r=dict(zip(cols,row))
            counts['protein_records']+=1
            a=decode(r.get('Abundances'),expected=19)
            n=decode(r.get('AbundancesNormalized'),expected=19)
            conf=decode(r.get('ProteinFDRConfidence'),'i')
            q=decode(r.get('Expqvalue'))
            values=[r.get(k) for k in meta]
            values[-2]=json.dumps(conf) if conf is not None else ''
            values[-1]=json.dumps(q) if q is not None else ''
            ia=sum(x is not None for x in a) if a else 0
            inn=sum(x is not None for x in n) if n else 0
            iw.writerow(values+[int(a is not None),ia,inn])
            if a is not None or n is not None:
                counts['records_with_vector']+=1
                statuses[str(r.get('IsMasterProtein'))]+=1
                ids=[r.get(k) for k in ['UniqueSequenceID','Accession','GeneSymbol']]
                rw.writerow(ids+(a or [None]*19));nw.writerow(ids+(n or [None]*19))
                saved={k:r.get(k) for k in meta+['Abundances','AbundancesNormalized','AbundancesCounts','FoundinSamples']}
                fg.write(json.dumps(saved,default=raw_encode)+'\n')
                if inn:counts['records_with_normalized_values']+=1
                if (r.get('UniquePeptidesCount') or 0)>=2:counts['vector_records_unique_peptides_ge2']+=1
                if a and n:
                    for j,(av,nv) in enumerate(zip(a,n)):
                        if av is not None and nv is not None and av>0:
                            assert nv>=0
                            factors[j].append(nv/av)
                genes=set((r.get('GeneSymbol') or '').replace(';',' ').split())
                fasta_gene=re.search(r'\bGN=([^\s]+)',r.get('Description') or '')
                if fasta_gene:genes.add(fasta_gene.group(1))
                if genes.intersection(FC):fc_hits.append(dict(zip(meta,values))|{'valid_normalized_slots':inn})
            if (i+1)%1000==0:print(json.dumps({'records_read':i+1,'quantified':counts['records_with_normalized_values'],'new_bytes':reader.downloaded}),flush=True)
    factors_summary=[{'slot':slots[j],'n_pairs':len(v),'minimum_ratio':min(v) if v else None,'maximum_ratio':max(v) if v else None} for j,v in enumerate(factors)]
    report={'source':SOURCE,'counts':dict(counts),'master_status_codes_among_vectors':dict(statuses),
            'normalization_ratio_by_slot':factors_summary,'fc_records_with_quantification':fc_hits,
            'diagnosis_assigned':False,'biological_replicates_aggregated':False,'disease_comparison_performed':False}
    (OUT/'quantification_extraction_summary.json').write_text(json.dumps(report,indent=2))
    (ACCESS/'quantification_ranges.json').write_text(json.dumps(reader.manifest(),indent=2))
    print(json.dumps(report,indent=2),flush=True)

if __name__=='__main__':run()
