"""Recalculate only statistics supported by the public aggregate tables."""
from pathlib import Path
import importlib.util
import json
import tempfile
import numpy as np
import pandas as pd
from common import bh_adjust

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT/'source_data'
OUT = ROOT/'results/tables'

def load_script(filename):
    spec = importlib.util.spec_from_file_location(filename.removesuffix('.py'), ROOT/'analysis'/filename)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module

def same(a,b):
    np.testing.assert_allclose(a,b,rtol=1e-6,atol=1e-9,equal_nan=True)

def main():
    summary = load_script('09_reproduce_summary_statistics.py')
    extended = load_script('08_extended_validation.py')
    blood = pd.read_csv(SRC/'Blood_cohort_effects.csv')
    expected = pd.read_csv(SRC/'Blood_meta.csv').set_index('Program')
    loco = pd.read_csv(SRC/'Figure_2C_leave_one_cohort_out.csv').set_index(['Program','Omitted cohort'])
    for program,d in blood.groupby('Program',sort=False):
        result = summary.meta(d['Hedges g'].to_numpy(),d.Variance.to_numpy())
        for key,col in [('estimate','Summary Hedges g'),('low','95% CI low'),('high','95% CI high'),('p','Hartung–Knapp P'),('I2','I² (%)')]:
            same(result[key], expected.loc[program,col])
        for cohort in d.Cohort:
            sub=d[d.Cohort!=cohort]
            same(summary.meta(sub['Hedges g'].to_numpy(),sub.Variance.to_numpy())['estimate'],
                 loco.loc[(program,cohort),'Summary Hedges g'])
    with tempfile.TemporaryDirectory() as tmp:
        extended.OUT=Path(tmp)
        csf=extended.csf_meta()
        frozen=pd.read_csv(OUT/'csf_kmezic_adjusted_meta.csv').set_index('protein').loc[csf.protein]
        inputs=pd.read_csv(SRC/'CSF_Kmezic_published_coefficients.csv').set_index('protein').loc[csf.protein]
        for col in ['discovery_beta','discovery_se','replication_beta','replication_se']:
            same(csf[col],inputs[col])
        for col in csf.select_dtypes(include='number').columns:
            same(csf[col],frozen[col])
        nerve=extended.external_nerve_validation()
        frozen_nerve=pd.read_csv(OUT/'cidp_external_nerve_validation.csv').set_index('feature').loc[nerve.feature]
        for col in ['proportion','ci_low','ci_high']:
            same(nerve[col],frozen_nerve[col])
    modules=pd.read_csv(SRC/'CIDP_module_effects.csv')
    for _,d in modules.groupby(['comparison','cell_group']):
        same(bh_adjust(d.p_value),d.fdr_within_celltype_comparison)
    corr=pd.read_csv(OUT/'cidp_donor_module_correlations.csv')
    same(bh_adjust(corr.p_value),corr.fdr_BH_10_pairs)
    clinical=pd.read_csv(OUT/'cidp_clinical_correlations_with_fdr.csv')
    same(bh_adjust(clinical.p_value),clinical.fdr_BH_15_tests)
    report={'release':'2.1.0','manuscript':'v37','status':'passed',
        'blood_meta_modules_recalculated':len(expected),'leave_one_cohort_out_estimates_checked':len(loco),
        'published_CSF_coefficient_pairs_pooled':len(csf),'published_52_of_55_Wilson_interval_checked':True,
        'nerve_module_P_value_corrections_checked':len(modules),
        'donor_correlation_P_value_corrections_checked':len(corr),
        'clinical_P_value_corrections_checked':len(clinical),
        'raw_matrix_pipeline_rerun':False,'donor_level_tests_rerun':False,
        'donor_bootstrap_rerun':False,'correlation_coefficients_recalculated':False}
    (ROOT/'metadata/PUBLIC_SUMMARY_REPRODUCTION.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report))

if __name__=='__main__': main()
