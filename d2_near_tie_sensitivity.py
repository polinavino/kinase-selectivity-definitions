"""
D2 near-tie sensitivity: for Klaeger compounds whose top1-top2 affinity gap is
below 0.05 pKd units, how much does the ratio R(k) swing between k=1 and k=2?
(Revision for Molecular Informatics pre-resubmission audit: this number was
quoted in desiderata.tex's D2 motivation paragraph but had no backing script.)
"""
import script_logging; script_logging.capture(__file__)
import os
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))


def ratio_selectivity(M, top_n, floor=5.0):
    out = []
    for row in M:
        s = np.sort(row)[::-1]
        off = s[top_n] if len(s) > top_n else floor
        out.append(s[0] - max(off, floor))
    return np.array(out)


M = pd.read_csv(os.path.join(HERE, "klaeger_matrix.csv"), index_col=0).values
sorted_desc = np.sort(M, axis=1)[:, ::-1]
gap = sorted_desc[:, 0] - sorted_desc[:, 1]

near_tie = gap < 0.05
print(f"Klaeger compounds with top1-top2 affinity gap < 0.05 pKd units: {near_tie.sum()} of {M.shape[0]}")

R1 = ratio_selectivity(M, 1)
R2 = ratio_selectivity(M, 2)
print(f"Among these {near_tie.sum()} compounds:")
print(f"  mean R(1) = {R1[near_tie].mean():.3f}  (k=1: reference is the tied second target, gap ~0)")
print(f"  mean R(2) = {R2[near_tie].mean():.3f}  (k=2: reference skips past the tie)")
print(f"  individual gaps: {np.round(gap[near_tie], 4).tolist()}")
