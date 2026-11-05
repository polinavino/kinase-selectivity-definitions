"""
Figure 4: representative Davis binding profiles for the three instability types.

Selects three compounds directly from the rank-instability analysis of the Davis
dataset (same parameter sweep as selectivity_analysis.py):
  - most stable   = smallest overall rank std across all 30 configurations,
  - most unstable = largest overall rank std,
  - most ratio-sensitive = largest rank std within the ratio family, among compounds
    whose two best-bound constructs are different genes (Davis profiles mutant and
    wild-type constructs of the same kinase, e.g. BRAF(V600E) and BRAF),
and plots each compound's kinase affinities sorted in descending order.
Produces binding_profiles.png (Figure 4).
"""
import script_logging; script_logging.capture(__file__)
import os
os.environ['KMP_DUPLICATE_LIB_OK'] = 'TRUE'

import numpy as np
from scipy.stats import rankdata
import pandas as pd
import matplotlib.pyplot as plt

BASELINE = 5.0
ACTIVE_THRESHOLD = 6.0

affinity = pd.read_csv("davis_affinity.csv")
matrix = affinity.pivot(index='Drug_Index', columns='Protein_Index', values='Affinity').fillna(BASELINE)

# Davis identifies compounds only by PubChem CID; name each by its first PubChem
# synonym that is not a CAS number (cache written by davis_klaeger_overlap.py)
import json, re
_cid = pd.read_csv("davis_drugs.csv").set_index("Drug_Index")["CID"]
_syn = json.load(open("davis_pubchem_synonyms_cache.json"))
def drug_name(idx):
    cid = str(int(_cid[idx]))
    if cid in _NO_SYNONYMS:
        return _NO_SYNONYMS[cid]
    return next((s for s in _syn.get(cid, []) if not re.fullmatch(r"\d+-\d+-\d", s)), f"CID {cid}")
# PubChem records with no synonyms, named from their structure
# (CID 9809715: C31H33N5O4, InChIKey XZXHXSATPCNXJR-ZIADKAODSA-N, nintedanib)
_NO_SYNONYMS = {"9809715": "Nintedanib"}
M = matrix.values
n_drugs, n_kinases = M.shape
construct = pd.read_csv("davis_proteins.csv").set_index("Protein_Index")["Gene_Name"].reindex(matrix.columns)
gene = construct.str.replace(r"\(.*$", "", regex=True).str.replace(r"-(phosphorylated|nonphosphorylated)$", "", regex=True).values


# Definitions copied verbatim from selectivity_analysis.py so the instability
# ranking (and hence the three selected drugs) reproduces exactly.
def s_score(profiles, threshold):
    return (profiles > threshold).astype(float).mean(axis=1)


def selectivity_entropy(profiles, baseline=BASELINE, epsilon=1e-10):
    shifted = np.maximum(profiles - baseline, 0)
    row_sums = shifted.sum(axis=1, keepdims=True)
    row_sums = np.where(row_sums == 0, epsilon, row_sums)
    p = shifted / row_sums
    log_p = np.where(p > 0, np.log2(np.where(p > 0, p, 1.0)), 0)
    return -(p * log_p).sum(axis=1)


def gini_selectivity(profiles, baseline=BASELINE):
    shifted = np.maximum(profiles - baseline, 0)
    ginis = []
    for row in shifted:
        row_sorted = np.sort(row)
        n = len(row_sorted)
        total = np.cumsum(row_sorted)[-1]
        if total == 0:
            ginis.append(0.0)
            continue
        ginis.append((2 * np.sum(np.arange(1, n + 1) * row_sorted)) / (n * total) - (n + 1) / n)
    return np.array(ginis)


def ratio_selectivity(profiles, top_n=1):
    ratios = []
    for row in profiles:
        sorted_desc = np.sort(row)[::-1]
        primary = sorted_desc[0]
        if len(sorted_desc) <= top_n:
            ratios.append(0.0)
            continue
        off_target = sorted_desc[top_n]
        ratios.append(primary - 5.0 if off_target <= 5.0 else primary - off_target)
    return np.array(ratios)


def to_ranks(scores):
    return rankdata(-np.asarray(scores), method='average')


s_thresholds  = np.arange(5.5, 8.25, 0.25)
ent_baselines = np.arange(5.0, 6.75, 0.25)
gini_baselines = np.arange(5.0, 6.75, 0.25)
ratio_top_ns  = list(range(1, 6))

# Negate so rank 1 = most selective consistently across all four families,
# matching selectivity_analysis.py (otherwise the combined std is meaningless).
s_ranks     = np.array([to_ranks(-s_score(M, t))             for t in s_thresholds])
ent_ranks   = np.array([to_ranks(-selectivity_entropy(M, b)) for b in ent_baselines])
gini_ranks  = np.array([to_ranks(gini_selectivity(M, b))     for b in gini_baselines])
ratio_ranks = np.array([to_ranks(ratio_selectivity(M, n))    for n in ratio_top_ns])

all_ranks = np.vstack([s_ranks, ent_ranks, gini_ranks, ratio_ranks])
rank_std = all_ranks.std(axis=0)

most_stable    = int(rank_std.argmin())
most_unstable  = int(rank_std.argmax())
top2 = np.argsort(-M, axis=1)[:, :2]
distinct = gene[top2[:, 0]] != gene[top2[:, 1]]
most_ratio_sen = int(np.where(distinct, ratio_ranks.std(axis=0), -np.inf).argmax())


def describe(idx):
    profile = M[idx]
    s = np.sort(profile)[::-1]
    n_active = int((profile > ACTIVE_THRESHOLD).sum())
    return profile, s, n_active


for label, idx in [("most stable", most_stable), ("most unstable", most_unstable),
                   ("most ratio-sensitive", most_ratio_sen)]:
    _, s, n_active = describe(idx)
    print(f"{label:22s}: Drug {idx} ({drug_name(idx)})  sigma={rank_std[idx]:.2f}  "
          f"n_active={n_active}  top1-top2 gap={s[0]-s[1]:.3f}  "
          f"top two: {construct.iloc[top2[idx, 0]]}, {construct.iloc[top2[idx, 1]]}")
print(f"compounds whose two best-bound constructs are the same gene: {int((~distinct).sum())} of {n_drugs}")

panels = [
    (most_stable,    "Most stable",          'tab:green'),
    (most_unstable,  "Most unstable",        '#ef5350'),
    (most_ratio_sen, "Most ratio-sensitive", '#f0ad4e'),
]

fig, axes = plt.subplots(1, 3, figsize=(18, 5))
for ax, (idx, title, color) in zip(axes, panels):
    profile, s, n_active = describe(idx)
    note = f"top1-top2 gap {s[0] - s[1]:.2f}"
    sorted_desc = np.sort(profile)[::-1]
    ax.bar(range(len(sorted_desc)), sorted_desc, color=color, width=1.0)
    ax.axhline(ACTIVE_THRESHOLD, ls='--', color='black', lw=1, label='pKd=6 threshold')
    ax.axhline(BASELINE, ls=':', color='gray', lw=1, label='assay floor')
    ax.set_ylim(4.8, 11)
    ax.set_title(f"{drug_name(idx)} — {title}\n({n_active} active kinases, {note})")
    ax.set_xlabel("Kinases (sorted by affinity)")
    ax.set_ylabel("pKd")
    ax.legend(loc='upper right')

plt.tight_layout()
plt.savefig("binding_profiles.png", dpi=150, bbox_inches='tight')
print("\nSaved binding_profiles.png")
