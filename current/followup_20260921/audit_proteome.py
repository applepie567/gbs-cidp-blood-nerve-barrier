"""Audit the public PXD056286 result database without assigning diagnoses.

The extraction is limited to study and workflow metadata. Protein values are
not interpreted before biological specimen and group mappings are available.
"""
from pathlib import Path
import argparse
import csv
import json
import re
import hashlib
import xml.etree.ElementTree as ET
from read_sqlite_ranges import Reader

URL='https://ftp.pride.ebi.ac.uk/pride/data/archive/2025/09/PXD056286/QExactiveHF02_18988_N1530.pdResult'

def run(access,out,remote):
    names=['StudyInformation','AnalysisDefinition','WorkflowInputFiles','StudyInformationWorkflowInputFiles',
           'Workflows','DataTypesColumnExtendedData','DataTypeColumnsCategoricalValues']
    if remote:
        rr=Reader(URL,access/'pdresult_blocks')
        (access/'pdresult_schema.json').write_text(json.dumps(rr.schema,indent=2))
        for name in names:
            (access/(name+'.json')).write_text(json.dumps(list(rr.table(name)),indent=2))
    schema=json.loads((access/'pdresult_schema.json').read_text())
    study=json.loads((access/'StudyInformation.json').read_text())
    definition=json.loads((access/'AnalysisDefinition.json').read_text())[0][1]
    xml=ET.fromstring(definition).find('StudyDefinition')
    assert xml is not None
    factors=xml.find('Factors')
    samples=xml.find('Samples')
    assert len(study)==len(samples)==19
    proteins=next(r for r in schema if r[0]=='table' and r[1]=='TargetProteins')
    rows=[]
    for row in study:
        rows.append(dict(study_id=row[0],sample=row[1],run_identifier=row[2],study_file_id=row[3],
                         sample_group=row[6],replicate_group=row[7],ratios=row[8],replicate_ratios=row[9],
                         diagnosis_mapping='Not specified',biological_replicate_mapping='Not specified',
                         source=URL))
    out.mkdir(exist_ok=True)
    with (out/'pxd056286_sample_mapping.csv').open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=rows[0]);writer.writeheader();writer.writerows(rows)
    checks=[
        ('Database readable','Yes','SQLite schema and study metadata extracted using verified HTTP byte ranges'),
        ('Instrument runs',len(study),'Run count is not a biological sample count'),
        ('Nonempty sample groups',sum(bool(r[6]) for r in study),'StudyInformation.SampleGroup'),
        ('Nonempty replicate groups',sum(bool(r[7]) for r in study),'StudyInformation.ReplicateGroup'),
        ('Study factors',len(factors),'AnalysisDefinition.StudyDefinition.Factors'),
        ('Samples with factor values',sum(bool(len(e.find('FactorValues'))) for e in samples),'AnalysisDefinition sample factor values'),
        ('Protein abundance fields','Present','TargetProteins.Abundances and AbundancesNormalized are BLOB fields'),
        ('Protein quantification decoded','No','Not attempted after the biological mapping gate failed'),
        ('Disease effect estimated','No','No verified run-to-specimen-to-diagnosis mapping in the inspected metadata'),
        ('Next usable input','Sample mapping and quantitative export','Need diagnosis, independent specimen IDs, technical replicate mapping, protein identifiers, quantitative values, and processing metadata'),
        ('Interpretation','Candidate resource remains usable in principle','Metadata availability, not demonstrated absence of target proteins, is the present obstacle')]
    with (out/'pxd056286_feasibility.csv').open('w',newline='') as f:
        writer=csv.writer(f);writer.writerow(['check','observed','interpretation']);writer.writerows(checks)
    assert all(r[6] is None and r[7] is None for r in study)
    assert len(factors)==0
    assert all(len(e.find('FactorValues'))==0 for e in samples)
    assert '"Abundances" BLOB' in proteins[4] and '"AbundancesNormalized" BLOB' in proteins[4]
    blocks=access/'pdresult_blocks'
    manifest={'accession':'PXD056286','date':'2026-09-21','source':URL,
              'scope':'Metadata audit only. No protein identification or quantitative differential analysis was performed.',
              'source_article':'https://doi.org/10.1007/s00401-025-02936-w',
              'evidence_files':[{'file':p.name,'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
                                for p in [access/'pdresult_schema.json']+[access/(n+'.json') for n in names]],
              'cached_ranges':[{'start':int(p.stem),'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
                               for p in sorted(blocks.glob('*.bin'))]}
    (out/'pxd056286_access_manifest.json').write_text(json.dumps(manifest,indent=2))
    print(json.dumps({'runs':len(study),'groups':0,'factors':0,'differential_analysis_performed':False}))

if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--access',type=Path,required=True)
    ap.add_argument('--out',type=Path,required=True)
    ap.add_argument('--remote',action='store_true',help='Fetch source metadata if a local snapshot is not already present')
    args=ap.parse_args();args.access.mkdir(parents=True,exist_ok=True)
    run(args.access,args.out,args.remote)
