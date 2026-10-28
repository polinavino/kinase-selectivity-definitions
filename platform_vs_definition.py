"""
Platform vs. definition: does the assay platform (Davis competition-binding vs.
Klaeger chemoproteomics) change compound rankings more or less than the choice
of selectivity definition? (Revision for Molecular Informatics, Reviewer 1,
major comment: "benchmark definitional instability against the ranking
variability introduced by assay platform.")

Uses the 46 compounds profiled in both Davis and Klaeger (resolved by PubChem
CID -> Klaeger drug-name match, davis_klaeger_overlap.py; this supersedes the
manuscript's previous, uninvestigated claim of "at most 19" shared compounds).
For this fixed set of 46 compounds, each selectivity definition is computed
separately on the Davis profile and on the Klaeger profile, and compounds are
ranked within this n=46 set (not the full 68/222-compound set) so that
platform and definition effects are compared on identical footing.
"""
import script_logging; script_logging.capture(__file__)
import os
import numpy as np
import pandas as pd
from scipy.stats import spearmanr, rankdata

HERE = os.path.dirname(os.path.abspath(__file__))


def to_ranks(s):
    return rankdata(-np.asarray(s), method='average')


def s_score(M, t):
    return -(M > t).astype(float).mean(axis=1)


def selectivity_entropy(M, b, eps=1e-10):
    sh = np.maximum(M - b, 0)
    rs = sh.sum(axis=1, keepdims=True)
    rs = np.where(rs == 0, eps, rs)
    p = sh / rs
    return -(p * np.where(p > 0, np.log2(p + eps), 0)).sum(axis=1)


def gini_selectivity(M, b):
    sh = np.maximum(M - b, 0)
    out = []
    for row in sh:
        rs = np.sort(row); n = len(rs); tot = rs.sum()
        out.append(0.0 if tot == 0 else (2 * np.sum(np.arange(1, n + 1) * rs)) / (n * tot) - (n + 1) / n)
    return np.array(out)


def ratio_selectivity(M, top_n, floor=5.0):
    out = []
    for row in M:
        s = np.sort(row)[::-1]
        off = s[top_n] if len(s) > top_n else floor
        out.append(s[0] - max(off, floor))
    return np.array(out)


DEFNS = ['s_score', 'entropy', 'gini', 'ratio']
S_THR = np.arange(5.5, 8.25, .25)
BASELINES = np.arange(5.0, 6.75, .25)


def median_ranks(M):
    ranks = {
        's_score': np.array([to_ranks(s_score(M, t)) for t in S_THR]),
        'entropy': np.array([to_ranks(-selectivity_entropy(M, b)) for b in BASELINES]),
        'gini':    np.array([to_ranks(gini_selectivity(M, b)) for b in BASELINES]),
        'ratio':   np.array([to_ranks(ratio_selectivity(M, k)) for k in range(1, 6)]),
    }
    return {d: np.median(ranks[d], axis=0) for d in DEFNS}


overlap = pd.read_csv(os.path.join(HERE, "davis_klaeger_overlap.csv"))

aff = pd.read_csv(os.path.join(HERE, "davis_affinity.csv"))
davis_full = aff.pivot(index='Drug_Index', columns='Protein_Index', values='Affinity').fillna(5.0)
klaeger_full = pd.read_csv(os.path.join(HERE, "klaeger_matrix.csv"), index_col=0)

davis_M = davis_full.loc[overlap['Davis_Drug_Index'].values].values
klaeger_M = klaeger_full.loc[overlap['Klaeger_Name'].values].values
n = len(overlap)
print(f"Shared Davis/Klaeger compounds (PubChem CID match): n = {n}")
print("(supersedes the manuscript's previous, unverified claim of 'at most 19')")

med_davis = median_ranks(davis_M)
med_klaeger = median_ranks(klaeger_M)

print(f"\nPlatform effect: Spearman r between Davis-rank and Klaeger-rank, same definition, "
      f"same {n} compounds:")
for d in DEFNS:
    r, p = spearmanr(med_davis[d], med_klaeger[d])
    print(f"  {d:8s}: r = {r:+.3f} (p = {p:.1e})")

print(f"\nDefinition effect: Spearman r between definitions, within platform, same {n} compounds:")
for platform, med in [("Davis", med_davis), ("Klaeger", med_klaeger)]:
    print(f"  {platform}:")
    for i in range(4):
        for j in range(i + 1, 4):
            r, _ = spearmanr(med[DEFNS[i]], med[DEFNS[j]])
            print(f"    {DEFNS[i]:8s} vs {DEFNS[j]:8s}: r = {r:+.3f}")

platform_rs = [spearmanr(med_davis[d], med_klaeger[d])[0] for d in DEFNS]
defn_rs = []
for med in [med_davis, med_klaeger]:
    for i in range(4):
        for j in range(i + 1, 4):
            defn_rs.append(spearmanr(med[DEFNS[i]], med[DEFNS[j]])[0])

print(f"\nSummary: platform-effect r range = {min(platform_rs):.2f}-{max(platform_rs):.2f} "
      f"(mean {np.mean(platform_rs):.2f})")
print(f"         definition-effect r range = {min(defn_rs):.2f}-{max(defn_rs):.2f} "
      f"(mean {np.mean(defn_rs):.2f})")
print("Lower r = more disagreement. Comparing the two ranges shows which source of "
      "ranking disagreement -- assay platform or scoring definition -- is larger on "
      "this shared-compound set.")
