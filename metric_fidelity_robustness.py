"""
Operationalization-sensitivity check (revision for Molecular Informatics, Reviewer C).

The manuscript computes selectivity entropy and the Gini coefficient on
baseline-subtracted pKd (a log-affinity scale).  The original papers use
different underlying quantities:

  * Uitdehaag & Zaman (2011) normalise association constants,
    p_i = Ka_i / sum_j Ka_j  with  Ka_i = 1/Kd_i = 10^(pKd_i - 9).
  * Graczyk (2007) computes the Gini coefficient of percent-inhibition values.

This script recomputes entropy and Gini using forms closer to the originals
(association-constant weighting) and compares the resulting selectivity
rankings with the manuscript's log-affinity forms.  The point is not that one
form is correct and the other wrong: it is that the two reasonable
operationalizations produce substantially different rankings, so the choice of
operationalization is itself a source of the definitional instability the paper
characterizes.  The manuscript adopts the baseline-subtracted log-affinity forms
because they place all four metrics on a single comparable scale and expose the
baseline parameter beta whose under-specification motivates property D3.
"""
import os
os.environ['KMP_DUPLICATE_LIB_OK'] = 'TRUE'
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))

def to_ranks(s):
    return len(s) - s.argsort().argsort()          # rank 1 = most selective

def entropy_literal(M, eps=1e-12):                  # Uitdehaag: p ~ Ka = 1/Kd
    Ka = np.power(10.0, M - 9.0)
    p = Ka / Ka.sum(1, keepdims=True)
    return -(p * np.where(p > 0, np.log2(p + eps), 0)).sum(1)   # low entropy = selective

def entropy_manuscript(M, b, eps=1e-10):            # manuscript: p ~ max(pKd - b, 0)
    sh = np.maximum(M - b, 0)
    rs = np.where(sh.sum(1, keepdims=True) == 0, eps, sh.sum(1, keepdims=True))
    p = sh / rs
    return -(p * np.where(p > 0, np.log2(p + eps), 0)).sum(1)

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
print(f"{'Dataset':10s}{'entropy: r(literal, manuscript)':>34s}{'Gini: r(literal, manuscript)':>32s}")
fig, axes = plt.subplots(2, 2, figsize=(10, 8.5))
for col, (name, M) in enumerate(datasets.items()):
    ent_lit = to_ranks(-entropy_literal(M))                                   # rank1=selective
    ent_man = np.median(np.array([to_ranks(-entropy_manuscript(M, b)) for b in base]), 0)
    r_ent, _ = spearmanr(ent_lit, ent_man)

    gini_lit = to_ranks(gini(np.power(10.0, M - 9.0)))                         # Gini of Ka
    gini_man = np.median(np.array([to_ranks(gini(np.maximum(M - b, 0))) for b in base]), 0)
    r_gini, _ = spearmanr(gini_lit, gini_man)

    print(f"{name:10s}{r_ent:>28.3f}{r_gini:>32.3f}")

    ax = axes[0, col]
    ax.scatter(ent_lit, ent_man, s=14, alpha=0.6, edgecolors='k', linewidths=0.2)
    ax.set_xlabel("entropy rank: literal (1/Kd)")
    ax.set_ylabel("entropy rank: manuscript (pKd-baseline)")
    ax.set_title(f"{name}: entropy selectivity ranking\nSpearman r = {r_ent:.2f}")

    ax = axes[1, col]
    ax.scatter(gini_lit, gini_man, s=14, alpha=0.6, edgecolors='k', linewidths=0.2,
               color='tab:green')
    ax.set_xlabel("Gini rank: literal (1/Kd)")
    ax.set_ylabel("Gini rank: manuscript (pKd-baseline)")
    ax.set_title(f"{name}: Gini selectivity ranking\nSpearman r = {r_gini:.2f}")
plt.tight_layout()
out = os.path.join(HERE, "metric_fidelity.png")
plt.savefig(out, dpi=150, bbox_inches='tight')
print(f"\nSaved {out}")
print("Two reasonable operationalizations of the same named metric give substantially")
print("different rankings (entropy even reverses on Klaeger): the operationalization")
print("choice is itself a source of definitional instability.")
