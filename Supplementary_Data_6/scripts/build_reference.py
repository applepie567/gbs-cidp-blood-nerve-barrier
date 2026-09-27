"""Aggregate every author-retained nucleus by donor and broad lineage."""
from pathlib import Path
import sys, gc, json
import numpy as np
import pandas as pd
from scipy import sparse

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'input/prior/External_nerve_20260926'))
from hdf_reader import read_cellbender

GROUPS={
 'Macrophages':['Macro1','Macro2'],
 'Schwann':['nmSC','mySC','repairSC','damageSC'],
 'Perineurial':['periC1','periC2','periC3'],
 'Endoneurial_stroma':['endoC'],
 'Epineurial_stroma':['epiC'],
 'Blood_endothelium':['ven_capEC1','ven_capEC2','artEC','venEC'],
 'Lymphatic_endothelium':['LEC'],
 'Mural':['PC1','PC2','VSMC'],
 'T_NK':['T_NK'],'B':['B'],'Mast':['Mast'],
 'Granulocytes':['Granulo'],'Adipocytes':['Adipo']}

def main():
    out=ROOT/'results'; out.mkdir(exist_ok=True)
    meta=pd.read_csv(ROOT/'input/GSE285983_metadata_all.csv.gz',usecols=['barcode','sample','cluster','level2','age','sex','center'])
    mapping={v:k for k,vs in GROUPS.items() for v in vs}
    assert set(meta.cluster)==set(mapping)
    meta['lineage']=meta.cluster.map(mapping)
    meta.to_csv(out/'atlas_nucleus_annotation.csv.gz',index=False)
    pd.DataFrame([dict(cluster=c,lineage=l) for c,l in mapping.items()]).to_csv(out/'lineage_mapping.csv',index=False)
    prior=pd.read_csv(ROOT/'input/prior/External_nerve_20260926/results/GSE285983_macrophage_Fc_counts.csv').set_index('sample_id')
    all_counts=[]; rows=[]; checks=[]; genes0=None
    for sample,d in prior.iterrows():
        path=ROOT/'input/nerve_h5'/d.h5_filename
        m=read_cellbender(path.read_bytes())
        genes=np.array([g.decode() if isinstance(g,bytes) else g for g in m['features/name']])
        barcodes=[b.decode() if isinstance(b,bytes) else b for b in m['barcodes']]
        if genes0 is None: genes0=genes
        assert np.array_equal(genes,genes0)
        sm=meta[meta['sample']==sample].set_index('barcode')
        annot={b.removeprefix(sample+'_'):mapping[c] for b,c in zip(sm.index,sm.cluster)}
        keep=[i for i,b in enumerate(barcodes) if b in annot]
        assert len(keep)==len(sm)
        group_names=list(GROUPS)
        gcode=[group_names.index(annot[barcodes[i]]) for i in keep]
        design=sparse.csr_matrix((np.ones(len(keep)),(keep,gcode)),shape=(len(barcodes),len(group_names)))
        mat=sparse.csc_matrix((m['data'],m['indices'],m['indptr']),shape=tuple(m['shape']))
        sums=np.asarray((mat@design).toarray()).T
        n=np.bincount(gcode,minlength=len(group_names))
        for j,group in enumerate(group_names):
            rows.append(dict(sample_id=sample,diagnosis=d.diagnosis,age=int(d.age),sex=d.sex,centre=d.centre,
                             lineage=group,n_nuclei=int(n[j]),library_counts=float(sums[j].sum()),
                             reference_eligible=bool(n[j]>=20)))
        gi={g:list(genes).index(g) for g in ['FCGR1A','FCGR2A','FCGR2B','FCGR3A','FCGR3B','FCGRT','FCER1G','TYROBP','SYK','LYN','HCK']}
        delta=max(abs(sums[0,idx]-float(d[g])) for g,idx in gi.items())
        assert delta==0 and sums[0].sum()==d.library_counts and n[0]==d.n_nuclei
        checks.append(dict(sample_id=sample,n_nuclei=len(keep),prior_Fc_count_max_difference=delta,prior_macrophage_library_difference=float(sums[0].sum()-d.library_counts)))
        all_counts.append(sums)
        print(sample,len(keep),'nuclei',flush=True)
        del m,mat,design,sums;gc.collect()
    counts=np.concatenate(all_counts)
    # Collapse duplicate gene symbols identically across all reference profiles.
    names,inv=np.unique(genes0,return_inverse=True)
    collapsed=np.zeros((len(rows),len(names)),dtype=np.float64)
    for j in range(len(rows)): np.add.at(collapsed[j],inv,counts[j])
    np.savez_compressed(out/'atlas_donor_lineage_counts.npz',counts=collapsed,genes=names)
    pd.DataFrame(rows).to_csv(out/'atlas_donor_lineage_metadata.csv',index=False)
    pd.DataFrame(checks).to_csv(out/'reference_reconstruction_check.csv',index=False)
    print('Reference complete',collapsed.shape,int(sum(x['n_nuclei'] for x in rows)),flush=True)

if __name__=='__main__': main()
