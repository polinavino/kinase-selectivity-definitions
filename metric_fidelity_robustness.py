"""
Operationalization-sensitivity check (revision for Molecular Informatics, Reviewer 1, comment C).

The manuscript computes selectivity entropy and the Gini coefficient on
baseline-subtracted pKd (a log-affinity scale).  The original papers use
different underlying quantities:

  * Uitdehaag & Zaman (2011) normalise association constants,
    p_i = Ka_i / sum_j Ka_j  with  Ka_i = 1/Kd_i = 10^(pKd_i - 9).
  * Graczyk (2010) computes the Gini coefficient of percent-inhibition values.

This script recomputes entropy and Gini using forms closer to the originals
(association-constant weighting) and compares the resulting selectivity
rankings with the manuscript's log-affinity forms.  The point is not that one
form is correct and the other wrong: it is that the two reasonable
operationalizations produce substantially different rankings, so the choice of
operationalization is itself a source of the definitional instability the paper
characterizes.  The manuscript adopts the baseline-subtracted log-affinity forms
because they place all four metrics on a single comparable scale and expose the
baseline parameter beta whose under-specification motivates property D2.
"""
import script_logging; script_logging.capture(__file__)
import os
os.environ['KMP_DUPLICATE_LIB_OK'] = 'TRUE'
import numpy as np
import pandas as pd
from scipy.stats import spearmanr, rankdata
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))

def to_ranks(s):
    return rankdata(-np.asarray(s), method='average')          # rank 1 = most selective

FLOOR = 5.0

def ka(M, floor_weight):
    """Association constants; floor-coded entries (no binding detected) kept at the
    floor value or given zero weight."""
    Ka = np.power(10.0, M - 9.0)
    return Ka if floor_weight else np.where(M > FLOOR, Ka, 0.0)

def entropy_literal(M, floor_weight=True):          # Uitdehaag: p ~ Ka = 1/Kd
    Ka = ka(M, floor_weight)
    tot = Ka.sum(1, keepdims=True)
    p = Ka / np.where(tot == 0, 1.0, tot)
    return -(p * np.where(p > 0, np.log2(np.where(p > 0, p, 1.0)), 0)).sum(1)   # low entropy = selective

def entropy_manuscript(M, b, eps=1e-10):            # manuscript: p ~ max(pKd - b, 0)
    sh = np.maximum(M - b, 0)
    rs = np.where(sh.sum(1, keepdims=True) == 0, eps, sh.sum(1, keepdims=True))
    p = sh / rs
    return -(p * np.where(p > 0, np.log2(np.where(p > 0, p, 1.0)), 0)).sum(1)

def gini(x):
    out = []
    for row in x:
        r = np.sort(np.maximum(row, 0)); n = len(r); t = r.sum()
        out.append(0.0 if t == 0 else (2 * np.sum(np.arange(1, n+1) * r)) / (n * t) - (n+1)/n)
    return np.array(out)

datasets = {}
aff = pd.read_csv(os.path.join(HERE, "davis_affinity.csv"))
datasets['Davis'] = aff.pivot(index='Drug_Index', columns='Protein_Index',
                              values='Affinity').fillna(5.0).values
datasets['Klaeger'] = pd.read_csv(os.path.join(HERE, "klaeger_matrix.csv"), index_col=0).values

base = np.arange(5.0, 6.75, 0.25)
print(f"{'Dataset':10s}{'floor-coded entries':>22s}{'entropy: r(literal, manuscript)':>34s}{'Gini: r(literal, manuscript)':>32s}")
fig, axes = plt.subplots(2, 2, figsize=(10, 8.5))
STYLE = {True: dict(color='#888888', marker='o', label='floor-coded entries at the floor'),
         False: dict(color='#2a78d6', marker='s', label='floor-coded entries given no weight')}
for col, (name, M) in enumerate(datasets.items()):
    ent_man = np.median(np.array([to_ranks(-entropy_manuscript(M, b)) for b in base]), 0)
    gini_man = np.median(np.array([to_ranks(gini(np.maximum(M - b, 0))) for b in base]), 0)
    r = {}
    for fw in (True, False):
        ent_lit = to_ranks(-entropy_literal(M, fw))                           # rank1=selective
        gini_lit = to_ranks(gini(ka(M, fw)))                                  # Gini of Ka
        r[fw] = (spearmanr(ent_lit, ent_man)[0], spearmanr(gini_lit, gini_man)[0])
        print(f"{name:10s}{'at the floor' if fw else 'no weight':>22s}{r[fw][0]:>28.4f}{r[fw][1]:>32.4f}")
        axes[0, col].scatter(ent_lit, ent_man, s=14, alpha=0.6, edgecolors='white', linewidths=0.3, **STYLE[fw])
        axes[1, col].scatter(gini_lit, gini_man, s=14, alpha=0.6, edgecolors='white', linewidths=0.3, **STYLE[fw])
    for row, what in enumerate(("entropy", "Gini")):
        ax = axes[row, col]
        ax.set_xlabel(f"{what} rank: association constants (1/Kd)")
        ax.set_ylabel(f"{what} rank: log-affinity (pKd - baseline)")
        ax.set_title(f"{name}: {what} ranking\nSpearman r = {r[True][row]:.2f} (grey), {r[False][row]:.2f} (blue)")
axes[0, 0].legend(fontsize=8, loc='lower right', frameon=False)
plt.tight_layout()
out = os.path.join(HERE, "metric_fidelity.png")
plt.savefig(out, dpi=150, bbox_inches='tight')
print(f"\nSaved {out}")
print("Both the weighting and the treatment of floor-coded entries change the rankings;")
print("the negative Klaeger correlations come from the floor-coded entries.")
