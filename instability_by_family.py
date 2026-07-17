"""
Figure 3 and the instability-predictor tables (Klaeger = Table 2, Davis = Table 3).

For each definition family (S-score, entropy, Gini, ratio) we sweep the
definition's parameter, compute each active compound's rank standard deviation
across the sweep (its "instability"), and correlate that instability against a
set of binding-profile features. The single best-predicting feature per family
is plotted for Klaeger (Figure 3). The full feature x family correlation table
is printed for both Klaeger (Table 2) and Davis (Table 3), so the predictor
structure can be seen to replicate on the second pK_d dataset.
"""
import script_logging; script_logging.capture(__file__)
import os
os.environ['KMP_DUPLICATE_LIB_OK'] = 'TRUE'

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
import matplotlib.pyplot as plt

ACTIVE_THRESHOLD = 6.0   # pKd > 6 counts as an active kinase
BASELINE = 5.0           # assay detection floor (pKd datasets)


def s_score(profiles, threshold):
    return -(profiles > threshold).astype(float).mean(axis=1)


def selectivity_entropy(profiles, baseline=BASELINE, epsilon=1e-10):
    shifted = np.maximum(profiles - baseline, 0)
    row_sums = np.where(shifted.sum(axis=1, keepdims=True) == 0, epsilon,
                        shifted.sum(axis=1, keepdims=True))
    p = shifted / row_sums
    return -(-(p * np.where(p > 0, np.log2(p + epsilon), 0)).sum(axis=1))


def gini_selectivity(profiles, baseline=BASELINE):
    shifted = np.maximum(profiles - baseline, 0)
    ginis = []
    for row in shifted:
        row_sorted = np.sort(row)
        n = len(row_sorted)
        total = row_sorted.sum()
        if total == 0:
            ginis.append(0.0)
            continue
        ginis.append((2 * np.sum(np.arange(1, n + 1) * row_sorted)) / (n * total) - (n + 1) / n)
    return np.array(ginis)


def ratio_selectivity(profiles, top_n=1):
    ratios = []
    for row in profiles:
        s = np.sort(row)[::-1]
        ratios.append(s[0] - max(s[top_n] if len(s) > top_n else BASELINE, BASELINE))
    return np.array(ratios)


def to_ranks(scores):
    return len(scores) - scores.argsort().argsort()


s_thresholds  = np.arange(5.5, 8.25, 0.25)
ent_baselines = np.arange(5.0, 6.75, 0.25)
gini_baselines = np.arange(5.0, 6.75, 0.25)
ratio_top_ns  = list(range(1, 6))

families = ['s_score', 'entropy', 'gini', 'ratio']
feature_cols = ['n_active', 'top1_top2_gap', 'active_std', 'active_range', 'top2_pkd']
colors = {'s_score': 'tab:blue', 'entropy': 'tab:orange', 'gini': 'tab:green', 'ratio': 'tab:red'}


def analyze(M, label):
    """Return (instability, feature_df, active-mask, corr, best) for one dataset."""
    n_drugs, n_kinases = M.shape
    instability = {
        's_score': np.array([to_ranks(s_score(M, t))             for t in s_thresholds]).std(axis=0),
        'entropy': np.array([to_ranks(selectivity_entropy(M, b)) for b in ent_baselines]).std(axis=0),
        'gini':    np.array([to_ranks(gini_selectivity(M, b))    for b in gini_baselines]).std(axis=0),
        'ratio':   np.array([to_ranks(ratio_selectivity(M, n))   for n in ratio_top_ns]).std(axis=0),
    }
    features = {}
    for i in range(n_drugs):
        profile = M[i]
        active = profile[profile > ACTIVE_THRESHOLD]
        s = np.sort(profile)[::-1]
        features[i] = {
            'n_active':      len(active),
            'top1_top2_gap': s[0] - s[1],
            'active_std':    active.std() if len(active) > 1 else 0.0,
            'active_range':  (active.max() - active.min()) if len(active) > 1 else 0.0,
            'top2_pkd':      s[1],
        }
    feat_df = pd.DataFrame(features).T
    mask = (M > ACTIVE_THRESHOLD).sum(axis=1) > 0
    n_active_drugs = int(mask.sum())

    corr = {c: {} for c in feature_cols}
    for col in feature_cols:
        for fam in families:
            corr[col][fam] = spearmanr(feat_df[col].values[mask], instability[fam][mask])
    best = {}
    for fam in families:
        col = max(feature_cols, key=lambda c: abs(corr[c][fam][0]))
        best[fam] = (col, *corr[col][fam])

    def stars(p):
        return '***' if p < 0.001 else '**' if p < 0.01 else '*' if p < 0.05 else ''
    print(f"\n=== {label}: {n_drugs} drugs x {n_kinases} kinases; active drugs n={n_active_drugs} ===")
    print(f"{'Feature':16s} " + "".join(f"{f:>14s}" for f in families))
    for col in feature_cols:
        row = f"{col:16s} "
        for fam in families:
            r, p = corr[col][fam]
            row += f"{f'{r:+.3f}{stars(p)}':>14s}"
        print(row)
    for fam in families:
        col, r, p = best[fam]
        print(f"best predictor {fam:10s}: {col:14s} r={r:+.3f} p={p:.3f}")

    pd.DataFrame({fam: {col: corr[col][fam][0] for col in feature_cols} for fam in families}) \
        .to_csv(f"instability_predictors_{label.lower()}.csv")
    return instability, feat_df, mask, corr, best


# ---- Klaeger (Table 2) + Figure 3 ----
Mk = pd.read_csv("klaeger_matrix.csv", index_col=0).values
instability, feat_df, mask, corr, best = analyze(Mk, "Klaeger")

fig, axes = plt.subplots(1, 4, figsize=(20, 5))
for ax, fam in zip(axes, families):
    col, r, p = best[fam]
    ax.scatter(feat_df[col].values[mask], instability[fam][mask],
               color=colors[fam], edgecolors='white', linewidths=0.5, alpha=0.8, s=60)
    ax.set_title(f"{fam}\nbest feature: {col}\nr={r:.2f}, p={p:.3f}")
    ax.set_xlabel(col)
    ax.set_ylabel("Rank std")
fig.suptitle("Best predictor of instability per definition family (Klaeger)", fontsize=14)
plt.tight_layout(rect=[0, 0, 1, 0.96])
plt.savefig("instability_by_family.png", dpi=150, bbox_inches='tight')
print("\nSaved instability_by_family.png")

# ---- Davis (Table 3) ----
davis = pd.read_csv("davis_affinity.csv")
Md = davis.pivot(index='Drug_Index', columns='Protein_Index', values='Affinity').fillna(5.0).values
analyze(Md, "Davis")
