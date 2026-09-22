"""Verify exported intensities, construct run mapping and audit Fc coverage."""
from pathlib import Path
import csv
import json
import re
import platform
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parent
ACCESS=ROOT/'access'
OUT=ROOT/'results'
FC=['FCGR1A','FCGR2A','FCGR2B','FCGR3A','FCGR3B','FCGRT','FCER1G','TYROBP','SYK','LYN','HCK']
SOURCE='https://ftp.pride.ebi.ac.uk/pride/data/archive/2025/09/PXD056286/QExactiveHF02_18988_N1530.pdResult'

def table(name,suffix='pdResult'):
    return json.loads((ACCESS/(suffix+'_'+name+'_with_ids.json')).read_text())

def main():
    inv=pd.read_csv(OUT/'pxd056286_protein_inventory.csv',dtype={'UniqueSequenceID':str})
    raw=pd.read_csv(OUT/'pxd056286_abundance_raw.csv',dtype={'UniqueSequenceID':str}).set_index('UniqueSequenceID')
    norm=pd.read_csv(OUT/'pxd056286_abundance_normalized.csv',dtype={'UniqueSequenceID':str}).set_index('UniqueSequenceID')
    assert inv.UniqueSequenceID.is_unique and raw.index.is_unique and norm.index.equals(raw.index)
    assert set(raw.index).issubset(set(inv.UniqueSequenceID))
    slots=['slot_'+str(i).zfill(2) for i in range(1,20)]
    boxes={map_id:sorted([r for _,r in table('DataDistributionBoxes') if r[1]==map_id],key=lambda x:x[2]) for map_id in (49,53)}
    assert len(boxes[49])==len(boxes[53])==19
    names=['F'+str(i) for i in range(1,20)]
    for i,(a,b) in enumerate(zip(boxes[49],boxes[53]),1):
        assert a[2]==b[2]==i and a[3]==b[3]==f'F{i}: Sample'
    multipliers={r[3].split(':')[0]:r[8] for _,r in table('ProcessingNodeCustomData') if r[4]=='ResultItemDataPurpose/QuanNormalizationValue'}
    a=raw[slots].to_numpy(float);b=norm[slots].to_numpy(float)
    assert np.array_equal(np.isnan(a),np.isnan(b))
    validation=[]
    for j,name in enumerate(names):
        use=np.isfinite(a[:,j])&np.isfinite(b[:,j])
        expected=a[use,j]*multipliers[name]
        assert np.allclose(expected,b[use,j],rtol=1e-12,atol=1e-8)
        relative=np.abs(expected-b[use,j])/np.maximum(np.abs(expected),1)
        validation.append({'quant_column':name,'normalization_multiplier':multipliers[name],
                           'compared_values':int(use.sum()),'max_relative_error':float(relative.max())})
    pd.DataFrame(validation).to_csv(OUT/'normalization_verification.csv',index=False)
    rename=dict(zip(slots,names))
    raw.rename(columns=rename).reset_index().to_csv(OUT/'PXD056286_raw_abundances.csv',index=False)
    norm.rename(columns=rename).reset_index().to_csv(OUT/'PXD056286_normalized_abundances.csv',index=False)

    keep=inv.IsMasterProtein.eq(0)&inv.ExcludedBy.eq(-1)&inv.GroupUniquePeptidesCount.ge(2)&inv.ProteinFDRConfidence.eq('[3]')
    native={str(r[0]) for _,r in table('TargetProteins_Subset_4','pdResultView')}
    assert set(inv.loc[keep,'UniqueSequenceID'])==native
    inv['passes_stored_filter']=keep.astype(int)
    inv['retained_quantification']=(keep&inv.valid_normalized_slots.gt(0)).astype(int)
    inv['master_status']=inv.IsMasterProtein.map({0:'Master',1:'Master candidate',2:'None',3:'Rejected master',4:'Rejected candidate'}).fillna('Not stored')
    inv['protein_confidence']=inv.ProteinFDRConfidence.map({'[3]':'High','[2]':'Medium','[1]':'Low','[0]':'Unknown'}).fillna('Not stored')
    inv['protein_q_value']=inv.Expqvalue.map(lambda v:json.loads(v)[0] if isinstance(v,str) and v else np.nan)
    inv['fasta_gene_symbol']=inv.Description.str.extract(r'\bGN=([^\s]+)',expand=False)
    quant_info=inv[inv.UniqueSequenceID.isin(raw.index)].copy()
    qc_cols=['UniqueSequenceID','Accession','GeneSymbol','fasta_gene_symbol','master_status','protein_confidence','protein_q_value',
             'UniquePeptidesCount','GroupUniquePeptidesCount','passes_stored_filter','valid_normalized_slots','ProteinGroupIDs','Description']
    quant_info[qc_cols].to_csv(OUT/'PXD056286_protein_QC.csv',index=False)
    inv.to_csv(OUT/'PXD056286_full_inventory.csv',index=False)
    retained_ids=inv.loc[inv.retained_quantification.eq(1),'UniqueSequenceID']
    for title,frame in [('raw',raw),('normalized',norm)]:
        frame.loc[retained_ids].rename(columns=rename).reset_index().to_csv(OUT/f'PXD056286_{title}_stored_filter.csv',index=False)

    def gene_tokens(row):
        tokens=set(str(row.GeneSymbol).replace(';',' ').split()) if pd.notna(row.GeneSymbol) else set()
        m=re.search(r'\bGN=([^\s]+)',str(row.Description))
        if m:tokens.add(m.group(1))
        return tokens
    gene_sets=inv.apply(gene_tokens,axis=1)
    coverage=[]
    for gene in FC:
        rows=inv[gene_sets.map(lambda s:gene in s)]
        measured=rows[rows.valid_normalized_slots.gt(0)]
        retained=rows[rows.retained_quantification.eq(1)]
        codes=' | '.join(measured.Accession.astype(str))
        if len(retained):reason='Retained. Only one component of the 11 gene Fc module.'
        elif not len(measured):reason='No stored quantitative values for this gene. This does not establish biological absence.'
        elif gene in ('FCGR3A','FCGR3B'):reason='One group-unique peptide. Stored filter requires at least two.'
        else:reason='Low confidence master candidate with one group-unique peptide. LYN and HCK share protein group 284.'
        coverage.append({'gene':gene,'stored_protein_records':len(rows),'quantified_records':len(measured),
          'retained_quantified_records':len(retained),'quantified_accessions':codes,
          'retained_accessions':' | '.join(retained.Accession.astype(str)),
          'maximum_quantified_run_count':int(measured.valid_normalized_slots.max()) if len(measured) else 0,
          'maximum_retained_run_count':int(retained.valid_normalized_slots.max()) if len(retained) else 0,'interpretation':reason})
    cov=pd.DataFrame(coverage);cov.to_csv(OUT/'PXD056286_Fc_coverage.csv',index=False)
    assert cov.retained_quantified_records.sum()==1
    assert cov.loc[cov.retained_quantified_records.gt(0),'gene'].tolist()==['FCGRT']

    old=pd.read_csv(ACCESS/'source_run_mapping.csv',keep_default_na=False)
    # Verify the inherited audit table against the newly retrieved native tables.
    study=table('StudyInformation')
    study_by_id={str(r[0]):r for _,r in study}
    workflow_by_id={r[2]:r for _,r in table('WorkflowInputFiles')}
    links={str(r[0]):r[3] for _,r in table('StudyInformationWorkflowInputFiles')}
    for _,row in old.iterrows():
        native_row=study_by_id[str(row.study_id)]
        assert [row['sample'],row.run_identifier,row.study_file_id]==native_row[1:4]
        assert native_row[4].replace('\\','/').rsplit('/',1)[-1]==row.run_identifier+'.raw'
        workflow_row=workflow_by_id[links[str(row.study_id)]]
        assert workflow_row[9]==row.study_file_id
        assert workflow_row[3].replace('\\','/').rsplit('/',1)[-1]==row.run_identifier+'.raw'
    header={r['run_identifier']:r for r in json.loads((ACCESS/'raw_header_review.json').read_text())}
    tail=old.run_identifier.str.rsplit('_',n=1).str[-1]
    counts=tail.value_counts()
    old['quant_column']=old.study_file_id
    old['normalization_multiplier']=old.study_file_id.map(multipliers)
    old['filename_code']=tail
    old['runs_with_same_filename_code']=tail.map(counts)
    old['raw_header_label']=[next(s for s in header[x]['relevant_strings'] if s.startswith('AH_CIDP_Nerv_')) for x in old.run_identifier]
    for name in ['biological_specimen_id_author','diagnosis_author','technical_replicate_group_author',
                 'publication_sample_label_author','included_in_publication_author','mapping_evidence_author']:
        old[name]=''
    old.to_csv(OUT/'PXD056286_run_mapping.csv',index=False)
    assert len(study)==19 and all(not r[6] and not r[7] for _,r in study)
    assert len(multipliers)==19 and all(np.isfinite(old.normalization_multiplier))
    summary={'protein_records_total':len(inv),'records_with_quantitative_values':len(raw),
       'stored_filter_records':int(keep.sum()),'stored_filter_records_with_values':len(retained_ids),
       'published_protein_count':1411,'published_count_exactly_reproduced':len(retained_ids)==1411,
       'instrument_runs':len(old),'distinct_filename_codes':len(counts),
       'filename_codes_repeated':counts[counts>1].to_dict(),
       'fc_genes_predefined':len(FC),'fc_genes_with_quantitative_values':int(cov.quantified_records.gt(0).sum()),
       'fc_genes_retained':int(cov.retained_quantified_records.gt(0).sum()),'fcg_rt_observed_runs':int(cov.loc[cov.gene.eq('FCGRT'),'maximum_retained_run_count'].iloc[0]),
       'native_filter_ID_set_identical':True,'normalization_values_verified':int(sum(v['compared_values'] for v in validation)),
       'run_mapping_matches_native_study_and_workflow_tables':True,
       'maximum_normalization_relative_error':max(v['max_relative_error'] for v in validation),
       'diagnosis_mapping_verified':False,'technical_replicate_mapping_verified':False,
       'independent_nerve_disease_validation_completed':False,
       'python':platform.python_version(),'pandas':pd.__version__,'numpy':np.__version__}
    (OUT/'validation_summary.json').write_text(json.dumps(summary,indent=2))
    print(json.dumps(summary,indent=2))

if __name__=='__main__':main()
