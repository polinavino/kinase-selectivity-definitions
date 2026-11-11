"""
SI Figure S3 and the profile-feature table SI Table S4 (Klaeger and Davis).

For each definition family (S-score, entropy, Gini, ratio) we sweep the
definition's parameter, compute each active compound's rank standard deviation
across the sweep (its "instability"), and correlate that instability against a
set of binding-profile features. The compounds are ranked among the compounds
with an active kinase only, those D1 allows to be scored, so a compound's
instability does not depend on how compounds with no active kinase are scored.
The same table with ranks taken among all compounds is printed as a
sensitivity check. The single best-predicting feature per definition
is plotted for Klaeger (SI Figure S3). The full feature x definition correlation table
is printed for both Klaeger and Davis (SI Table S4), so the predictor
structure can be seen to replicate on the second pK_d dataset.
"""
import script_logging; script_logging.capture(__file__)
import os
os.environ['KMP_DUPLICATE_LIB_OK'] = 'TRUE'

import numpy as np
import pandas as pd
from scipy.stats import spearmanr, rankdata, skew
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
    return -(-(p * np.where(p > 0, np.log2(np.where(p > 0, p, 1.0)), 0)).sum(axis=1))


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
    return rankdata(-np.asarray(scores), method='average')


s_thresholds  = np.arange(5.5, 8.25, 0.25)
ent_baselines = np.arange(5.0, 6.75, 0.25)
gini_baselines = np.arange(5.0, 6.75, 0.25)
ratio_top_ns  = list(range(1, 6))

families = ['s_score', 'entropy', 'gini', 'ratio']
feature_cols = ['n_active', 'top1_top2_gap', 'active_std', 'active_range', 'top2_pkd', 'active_skew']
colors = {'s_score': 'tab:blue', 'entropy': 'tab:orange', 'gini': 'tab:green', 'ratio': 'tab:red'}


def family_instability(M):
    return {
        's_score': np.array([to_ranks(s_score(M, t))             for t in s_thresholds]).std(axis=0),
        'entropy': np.array([to_ranks(selectivity_entropy(M, b)) for b in ent_baselines]).std(axis=0),
        'gini':    np.array([to_ranks(gini_selectivity(M, b))    for b in gini_baselines]).std(axis=0),
        'ratio':   np.array([to_ranks(ratio_selectivity(M, n))   for n in ratio_top_ns]).std(axis=0),
    }


def analyze(M, label):
    """Return (instability, feature_df, active-mask, corr, best) for one dataset."""
    n_drugs, n_kinases = M.shape
    act = (M > ACTIVE_THRESHOLD).any(1)
    # ranks among the compounds with an active kinase (D1-eligible); NaN elsewhere
    inst_act = family_instability(M[act])
    instability = {f: np.full(n_drugs, np.nan) for f in inst_act}
    for f in inst_act:
        instability[f][act] = inst_act[f]
    instability_all = family_instability(M)          # sensitivity: ranks among all compounds
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
            'active_skew':   skew(active) if len(active) > 2 else 0.0,
        }
    feat_df = pd.DataFrame(features).T
    mask = (M > ACTIVE_THRESHOLD).sum(axis=1) > 0
    n_active_drugs = int(mask.sum())

    corr = {c: {} for c in feature_cols}
    corr_all = {c: {} for c in feature_cols}
    for col in feature_cols:
        for fam in families:
            corr[col][fam] = spearmanr(feat_df[col].values[mask], instability[fam][mask])
            corr_all[col][fam] = spearmanr(feat_df[col].values[mask], instability_all[fam][mask])
    best = {}
    for fam in families:
        col = max(feature_cols, key=lambda c: abs(corr[c][fam][0]))
        best[fam] = (col, *corr[col][fam])

    def stars(p):
        return '***' if p < 0.001 else '**' if p < 0.01 else '*' if p < 0.05 else ''
    print(f"\n=== {label}: {n_drugs} drugs x {n_kinases} kinases; active drugs n={n_active_drugs}, "
          f"ranked among themselves ===")
    print(f"{'Feature':16s} " + "".join(f"{f:>14s}" for f in families))
    for col in feature_cols:
        row = f"{col:16s} "
        for fam in families:
            r, p = corr[col][fam]
            row += f"{f'{r:+.3f}{stars(p)}':>14s}"
        print(row)
    if n_active_drugs < n_drugs:
        # check: only compounds whose entropy and Gini are defined throughout the baseline
        # sweep (a kinase above its top, 6.5), re-ranked among themselves
        full_def = M.max(1) > ent_baselines.max()
        inst_d = family_instability(M[full_def])
        na_d = (M[full_def] > ACTIVE_THRESHOLD).sum(1)
        print(f"  check: the {int(full_def.sum())} compounds defined throughout the baseline sweep, re-ranked: "
              + "; ".join(f"n_active vs {fam} instability r = {spearmanr(na_d, inst_d[fam])[0]:+.3f} "
                          f"(p = {spearmanr(na_d, inst_d[fam])[1]:.2g})" for fam in ('entropy', 'gini')))
        print(f"  sensitivity: the same {n_active_drugs} compounds, ranked among all {n_drugs} compounds")
        for col in feature_cols:
            row = f"  {col:14s} "
            for fam in families:
                r, p = corr_all[col][fam]
                row += f"{f'{r:+.3f}{stars(p)}':>14s}"
            print(row)
    for fam in families:
        col, r, p = best[fam]
        print(f"best predictor {fam:10s}: {col:14s} r={r:+.3f} p={p:.3f}")

    pd.DataFrame({fam: {col: corr[col][fam][0] for col in feature_cols} for fam in families}) \
        .to_csv(f"instability_predictors_{label.lower()}.csv")
    return instability, feat_df, mask, corr, best


# ---- Klaeger (SI Table S4) + SI Figure S3 ----
Mk = pd.read_csv("klaeger_matrix.csv", index_col=0).values
instability, feat_df, mask, corr, best = analyze(Mk, "Klaeger")

FAM_NAME = {'s_score': 'S-score', 'entropy': 'Entropy', 'gini': 'Gini', 'ratio': 'Ratio'}
FEAT_NAME = {'n_active': 'number of active kinases', 'top1_top2_gap': 'top1-top2 gap (pKd)',
             'active_std': 'active-affinity std (pKd)', 'active_range': 'active-affinity range (pKd)',
             'top2_pkd': 'second-ranked pKd', 'active_skew': 'active-affinity skewness'}
plt.rcParams.update({'font.size': 13, 'axes.titlesize': 14})
fig, axes = plt.subplots(2, 2, figsize=(10, 8))
for ax, fam in zip(axes.ravel(), families):
    col, r, p = best[fam]
    ax.scatter(feat_df[col].values[mask], instability[fam][mask],
               color=colors[fam], edgecolors='white', linewidths=0.5, alpha=0.8, s=60)
    ax.set_title(f"{FAM_NAME[fam]}\nr={r:.2f}, " + ("p<0.001" if p < 0.001 else f"p={p:.3f}"))
    ax.set_xlabel(FEAT_NAME[col])
    ax.set_ylabel("Rank std")
fig.suptitle("Within-definition instability and profile features\n"
             "Klaeger: active compounds ranked among themselves", fontsize=14)
plt.tight_layout(rect=[0, 0, 1, 0.96])
script_logging.save_figure(fig, "instability_by_family.png", bbox_inches='tight')
print("\nSaved instability_by_family.png")

# ---- Davis (SI Table S4, right) ----
davis = pd.read_csv("davis_affinity.csv")
Md = davis.pivot(index='Drug_Index', columns='Protein_Index', values='Affinity').fillna(5.0).values
analyze(Md, "Davis")
