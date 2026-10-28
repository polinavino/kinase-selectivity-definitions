"""
Zero-active compounds and the entropy/Gini 0/0 convention (revision for
Molecular Informatics, Reviewer 1, major comment on undefined entries).

Entropy and Gini are undefined when no kinase in a compound's profile exceeds
the baseline beta (the normalized weights are 0/0). The current pipeline
resolves this by substituting a small epsilon in the denominator, which drives
the score to exactly 0 for both entropy and Gini and lets the compound be
ranked (tied with every other degenerate compound) rather than excluded. This
script states that convention explicitly and reports the alternative Reviewer
1 asked for: the zero-active vs. active rank-sigma contrast recomputed using,
for every compound, only the baseline sweep configurations in which its
entropy/Gini score is actually defined (max affinity > beta), rather than
pooling every configuration regardless of degeneracy.
"""
import script_logging; script_logging.capture(__file__)
import os
os.environ['KMP_DUPLICATE_LIB_OK'] = 'TRUE'
import numpy as np
import pandas as pd
from scipy.stats import rankdata, mannwhitneyu

HERE = os.path.dirname(os.path.abspath(__file__))


def to_ranks(scores):
    return rankdata(-np.asarray(scores), method='average')


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


M = pd.read_csv(os.path.join(HERE, "klaeger_matrix.csv"), index_col=0).values
n_drugs = M.shape[0]
baselines = np.arange(5.0, 6.75, .25)
max_aff = M.max(axis=1)
active_6 = max_aff > 6.0          # the paper's "has an active kinase" definition
print(f"Zero-active compounds (pKd>6 threshold): {(~active_6).sum()} of {n_drugs}")

# ---- convention actually used in the paper: epsilon -> score forced to 0 ----
ent_ranks_pooled = np.array([to_ranks(-selectivity_entropy(M, b)) for b in baselines])
gini_ranks_pooled = np.array([to_ranks(gini_selectivity(M, b)) for b in baselines])
for name, ranks in [('entropy', ent_ranks_pooled), ('gini', gini_ranks_pooled)]:
    sigma = ranks.std(axis=0)
    u, p = mannwhitneyu(sigma[~active_6], sigma[active_6])
    print(f"Current convention -- {name}: zero-active mean sigma={sigma[~active_6].mean():.1f}  "
          f"active mean sigma={sigma[active_6].mean():.1f}  (Mann-Whitney p={p:.2e})")
pooled = np.vstack([ent_ranks_pooled, gini_ranks_pooled])
sigma_pooled = pooled.std(axis=0)
u, p = mannwhitneyu(sigma_pooled[~active_6], sigma_pooled[active_6])
print(f"\nPooled convention (degenerate configs scored as entropy=Gini=0, included in every sweep):")
print(f"  mean rank sigma, zero-active: {sigma_pooled[~active_6].mean():.1f}  "
      f"active: {sigma_pooled[active_6].mean():.1f}  (Mann-Whitney p={p:.2e})")

# ---- alternative: only count a (compound, beta) configuration if it is DEFINED ----
# (max affinity > beta, i.e. at least one kinase contributes positive weight)
defined = max_aff[None, :] > baselines[:, None]   # shape (7, n_drugs), per-beta
n_defined_per_compound = defined.sum(axis=0)
print(f"\nPer-compound number of DEFINED configurations (of 7 baselines swept): "
      f"min={n_defined_per_compound.min()}, max={n_defined_per_compound.max()}")
print(f"  Zero-active (pKd>6) compounds have a defined score at {n_defined_per_compound[~active_6].min()}"
      f"-{n_defined_per_compound[~active_6].max()} of 7 baselines "
      f"(mean {n_defined_per_compound[~active_6].mean():.1f})")

print(f"\nAlternative convention (rank sigma computed only over each compound's own DEFINED configs):")
for name, ranks in [('entropy', ent_ranks_pooled), ('gini', gini_ranks_pooled)]:
    sigma_defined = np.full(n_drugs, np.nan)
    for i in range(n_drugs):
        rows = np.where(defined[:, i])[0]
        if len(rows) < 2:
            continue
        sigma_defined[i] = ranks[rows, i].std()
    valid = ~np.isnan(sigma_defined)
    zero_valid = valid & (~active_6)
    act_valid = valid & active_6
    u2, p2 = mannwhitneyu(sigma_defined[zero_valid], sigma_defined[act_valid])
    print(f"  {name}: {valid.sum()} of {n_drugs} compounds have >=2 defined configs "
          f"({(~valid).sum()} excluded -- defined at <2 of 7 baselines, cannot be scored at all); "
          f"zero-active mean sigma (n={zero_valid.sum()})={sigma_defined[zero_valid].mean():.1f}  "
          f"active mean sigma (n={act_valid.sum()})={sigma_defined[act_valid].mean():.1f}  (p={p2:.2e})")

sigma_defined = np.full(n_drugs, np.nan)
for i in range(n_drugs):
    rows = np.where(defined[:, i])[0]
    if len(rows) < 2:
        continue  # cannot compute a within-compound std across <2 defined configs
    vals = np.concatenate([ent_ranks_pooled[rows, i], gini_ranks_pooled[rows, i]])
    sigma_defined[i] = vals.std()

valid = ~np.isnan(sigma_defined)
zero_valid = valid & (~active_6)
act_valid = valid & active_6
if zero_valid.sum() >= 1 and act_valid.sum() >= 1:
    u2, p2 = mannwhitneyu(sigma_defined[zero_valid], sigma_defined[act_valid])
    print(f"  pooled (entropy+gini): {valid.sum()} of {n_drugs} computable; "
          f"mean rank sigma, zero-active (n={zero_valid.sum()}): {sigma_defined[zero_valid].mean():.1f}  "
          f"active (n={act_valid.sum()}): {sigma_defined[act_valid].mean():.1f}  (Mann-Whitney p={p2:.2e})")
print("\nDone.")
