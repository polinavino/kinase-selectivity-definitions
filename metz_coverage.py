"""
Metz per-compound kinase coverage distribution among the 704 retained
compounds (revision for Molecular Informatics, Reviewer 1, comment on the
Metz inclusion threshold: "50 of 172 kinases is 29%... report the coverage
distribution"). Methods states the resulting range/mean/median; this script
is the saved, reproducible source for those numbers (previously computed in
an ad hoc shell command, not persisted to a script).
"""
import script_logging; script_logging.capture(__file__)
import os
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))

mz = pd.read_excel(os.path.join(HERE, "metz.xls"), sheet_name=0, header=0)
meta = ['Cmpd_ID', 'PUBCHEM_SID', 'Canonical_Smiles', 'External_Cmpd_ID', 'External_Source',
        'Cluster', 'ClusterSize', 'Cluster_MCSS', 'Molecular_Weight', 'ALogP',
        'Num_H_Acceptors', 'Num_H_Donors', 'tPSA', 'Promiscuity_1uM']
kin_cols = [c for c in mz.columns if c not in meta]
sub = mz[kin_cols].apply(pd.to_numeric, errors='coerce')
n_kin = len(kin_cols)
coverage = (~np.isnan(sub.values)).sum(axis=1)
keep = coverage >= 50

print(f"Metz: {len(sub)} total compounds, {n_kin} kinase columns")
print(f"Compounds with >= 50 of {n_kin} kinases measured (retained): {keep.sum()}")

cov_kept = coverage[keep]
pct = 100 * cov_kept / n_kin
print(f"\nCoverage among the {keep.sum()} retained compounds:")
print(f"  count (kinases): min={cov_kept.min()}, max={cov_kept.max()}, "
      f"mean={cov_kept.mean():.1f}, median={np.median(cov_kept):.0f}")
print(f"  percent of {n_kin}: min={pct.min():.1f}%, max={pct.max():.1f}%, "
      f"mean={pct.mean():.1f}%, median={np.median(pct):.1f}%")
qs = np.percentile(pct, [10, 25, 50, 75, 90])
print(f"  deciles (10/25/50/75/90): {np.round(qs, 1).tolist()}")
