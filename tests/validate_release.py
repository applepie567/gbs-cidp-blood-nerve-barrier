"""Check source/figure alignment without claiming a primary-matrix rerun."""
from pathlib import Path
from io import BytesIO
import json,zipfile
import numpy as np
import pandas as pd
from PIL import Image
R=Path(__file__).resolve().parents[1]
manifest=json.loads((R/'metadata/WORKBOOK_CONTENTS.json').read_text())
with pd.ExcelFile(R/'source_data/Additional_file_1_source_data_public.xlsx') as xl:
    assert xl.sheet_names==[x['sheet'] for x in manifest]
    for x in manifest:
        a=pd.read_csv(R/x['csv'],keep_default_na=False)
        b=pd.read_excel(xl,sheet_name=x['sheet'],keep_default_na=False)
        assert list(a.columns)==list(b.columns),(x['sheet'],'columns')
        assert a.shape==b.shape,(x['sheet'],'shape')
        for col in a.columns:
            an=pd.to_numeric(a[col].replace('',np.nan),errors='coerce')
            bn=pd.to_numeric(b[col].replace('',np.nan),errors='coerce')
            if (an.notna()|a[col].eq('')).all():
                np.testing.assert_allclose(an,bn,rtol=1e-10,atol=1e-12,equal_nan=True)
            else:
                assert a[col].astype(str).tolist()==b[col].astype(str).tolist(),(x['sheet'],col)
index=pd.read_csv(R/'source_data/Figure_table_index.csv')
assert all((R/p).is_file() for p in index['Data file']) and len(index)==21
names=['Figure_1_study_architecture','Figure_2_acute_GBS_blood','Figure_3_GBS_CSF','Figure_4_CIDP_nerve','Figure_5_cross_compartment_genetics','Figure_S1_genetic_localization_bootstrap','Figure_S2_published_nerve_evidence']
man=next((R/'docs').glob('GBS_CIDP_compartmentalization_v37*.docx'))
with zipfile.ZipFile(man) as z:
    for i,name in enumerate(names,1):
        a=np.asarray(Image.open(BytesIO(z.read(f'word/media/image{i}.png'))).convert('RGB'))
        b=np.asarray(Image.open(R/'figures'/f'{name}.png').convert('RGB'))
        assert np.array_equal(a,b),(name,'differs from v37')
    assert 'GSE304871' not in z.read('word/document.xml').decode()
g=pd.read_csv(R/'source_data/Figure_5D_genetic_expression_contrasts.csv')
assert len(g)==40 and (g.fdr_within_cell_group>=.05).all()
assert (pd.read_csv(R/'source_data/Blood_meta.csv')['Hartung–Knapp P']>=.05).all()
report={'release':'2.1.0','manuscript':'v37','source_workbook_sheets':len(manifest),
        'figure_table_index_rows':len(index),'embedded_images_matched':7,
        'raw_matrix_pipeline_rerun':False,
        'workbook_matches_csv_sources':True,'status':'passed'}
(R/'metadata/VALIDATION_REPORT.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report))
