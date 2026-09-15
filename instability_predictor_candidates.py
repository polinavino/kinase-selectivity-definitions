"""Candidate profile-shape features considered for Table 2/3 but not added.

instability_by_family.py correlates six binding-profile features against
within-family rank instability. Before settling on that set (five original
features plus active-affinity skew, added after this script's results),
several other candidates motivated by prior selectivity-metric literature were
tested: the concentration ratio and Herfindahl-Hirschman Index (CR1, HHI --
the same economics-of-concentration measures Graczyk 2010 drew the Gini
coefficient from), the raw primary-target affinity (top1_pkd), a near-tie
count generalizing the top1-top2 gap, and the coefficient of variation of the
active-affinity distribution. This script reports why each was excluded:
CR1/HHI closely track n_active already in the table (shown here via their
direct correlation with n_active) rather than adding independent information,
and the remaining three are uniformly weaker predictors than the features
already in Table 2/3.
"""
import script_logging; script_logging.capture(__file__)
import os
os.environ['KMP_DUPLICATE_LIB_OK'] = 'TRUE'

import numpy as np
import pandas as pd
from scipy.stats import spearmanr, rankdata, skew

ACTIVE_THRESHOLD = 6.0
BASELINE = 5.0


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
            ginis.append(0.0); continue
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

CANDIDATES = ['cr1_share', 'hhi', 'top1_pkd', 'near_tie_1pkd', 'active_cv']


def analyze(M, label):
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
        shifted_active = np.maximum(active - BASELINE, 0)
        s = np.sort(profile)[::-1]
        total_mass = shifted_active.sum()
        shares = shifted_active / total_mass if total_mass > 0 else np.zeros_like(shifted_active)
        features[i] = {
            'n_active':      len(active),
            'top1_pkd':      s[0],
            'cr1_share':     shares.max() if len(shares) > 0 else 0.0,
            'hhi':           np.sum(shares ** 2) if len(shares) > 0 else 0.0,
            'near_tie_1pkd': int(np.sum(s[0] - s < 1.0)) - 1,
            'active_cv':     (active.std() / active.mean()) if len(active) > 1 and active.mean() != 0 else 0.0,
        }
    feat_df = pd.DataFrame(features).T
    mask = (M > ACTIVE_THRESHOLD).sum(axis=1) > 0
    n_active_drugs = int(mask.sum())

    def stars(p):
        return '***' if p < 0.001 else '**' if p < 0.01 else '*' if p < 0.05 else ''

    print(f"\n=== {label}: {n_drugs} drugs x {n_kinases} kinases; active drugs n={n_active_drugs} ===")
    r_hhi, p_hhi = spearmanr(feat_df['hhi'].values[mask], feat_df['n_active'].values[mask])
    r_cr1, p_cr1 = spearmanr(feat_df['cr1_share'].values[mask], feat_df['n_active'].values[mask])
    print(f"HHI vs n_active:       r={r_hhi:+.3f}{stars(p_hhi)}  (redundancy check)")
    print(f"CR1 vs n_active:       r={r_cr1:+.3f}{stars(p_cr1)}  (redundancy check)")
    print(f"{'Feature':16s} " + "".join(f"{f:>14s}" for f in families))
    for col in CANDIDATES:
        row = f"{col:16s} "
        for fam in families:
            r, p = spearmanr(feat_df[col].values[mask], instability[fam][mask])
            row += f"{f'{r:+.3f}{stars(p)}':>14s}"
        print(row)


Mk = pd.read_csv("klaeger_matrix.csv", index_col=0).values
analyze(Mk, "Klaeger")

davis = pd.read_csv("davis_affinity.csv")
Md = davis.pivot(index='Drug_Index', columns='Protein_Index', values='Affinity').fillna(5.0).values
analyze(Md, "Davis")
