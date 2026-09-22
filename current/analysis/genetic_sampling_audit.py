"""Audit the number of nuclei in source CIDP B cell donor groups.

This checks annotation counts only and does not reconstruct expression.
Output contains aggregate counts without donor identifiers.
"""
import argparse
from pathlib import Path
import pandas as pd


def audit(metadata, out):
    md = pd.read_csv(metadata, usecols=['sample', 'level2', 'cluster'], dtype=str)
    cidp = md[md.level2.eq('CIDP')]
    donors = sorted(cidp['sample'].unique())
    counts = cidp[cidp.cluster.eq('B')].groupby('sample').size().reindex(donors, fill_value=0)
    result = pd.DataFrame([dict(
        disease='CIDP', cell_group='B', n_donors=len(donors),
        n_below_20=int((counts < 20).sum()),
        min_nuclei=int(counts.min()), max_nuclei=int(counts.max()),
        total_nuclei=int(counts.sum()),
        interpretation='Sampling audit only, not a new genetic or expression test'
    )])
    out.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(out, index=False)
    print(result.to_string(index=False))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--metadata', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    audit(args.metadata, args.out)
