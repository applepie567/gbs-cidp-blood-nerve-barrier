"""Check that the distributable tree contains aggregate tables and current artifacts."""
from pathlib import Path
import csv, json

ROOT=Path(__file__).resolve().parents[1]
distribution=json.loads((ROOT/'metadata/PUBLIC_DISTRIBUTION.json').read_text())
for rel in distribution['individual_level_files_omitted']:
    assert not (ROOT/rel).exists(), f'Individual-level file remains in the release: {rel}'
for folder in ['source_data','results','tables']:
    for path in (ROOT/folder).rglob('*.csv'):
        with path.open(newline='') as f: header=next(csv.reader(f))
        forbidden={'sample','sample_id','donor','donor_id','patient','patient_id','age','sex','center','INCAT','disease_duration'}
        assert not forbidden.intersection(header),(str(path.relative_to(ROOT)),header)
assert not list((ROOT/'results').glob('*.json')), 'Review sample-level JSON before public distribution'
assert len(list((ROOT/'docs').glob('GBS_CIDP*.docx')))==1
assert len(list((ROOT/'source_data').glob('*.xlsx')))==1
print('Public distribution: no omitted individual-level files or individual-record columns; one current manuscript and workbook.')
