"""
Sensitivity of the SI Table S4 instability predictors to the placeholder value given
to undefined profile features.

The skewness of a compound's active affinities is undefined with fewer than three
active kinases, and their standard deviation and range with fewer than two.
instability_by_family.py assigns these compounds the value 0 (the vertical
clusters at exactly zero in SI Figure S3). Here each such correlation is recomputed
over only the compounds for which the feature is defined, to confirm that the
reported associations do not depend on that placeholder.

Scoring functions and parameter grids are copied from instability_by_family.py.
"""
import script_logging; script_logging.capture(__file__)
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


s_thresholds = np.arange(5.5, 8.25, 0.25)
baselines = np.arange(5.0, 6.75, 0.25)
families = ['s_score', 'entropy', 'gini', 'ratio']
# feature -> minimum number of active kinases for it to be defined
FEATURES = {'active_std': 2, 'active_range': 2, 'active_skew': 3}


def run(M, label):
    # ranks among the compounds with an active kinase, as in SI Table S4 (instability_by_family.py)
    act = (M > ACTIVE_THRESHOLD).any(1)
    MA = M[act]
    inst_a = {
        's_score': np.array([to_ranks(s_score(MA, t)) for t in s_thresholds]).std(axis=0),
        'entropy': np.array([to_ranks(selectivity_entropy(MA, b)) for b in baselines]).std(axis=0),
        'gini':    np.array([to_ranks(gini_selectivity(MA, b)) for b in baselines]).std(axis=0),
        'ratio':   np.array([to_ranks(ratio_selectivity(MA, n)) for n in range(1, 6)]).std(axis=0),
    }
    instability = {f: np.full(len(M), np.nan) for f in inst_a}
    for f in inst_a:
        instability[f][act] = inst_a[f]
    actives = [row[row > ACTIVE_THRESHOLD] for row in M]
    n_active = np.array([len(a) for a in actives])
    feats = {
        'active_std':   np.array([a.std() if len(a) > 1 else 0.0 for a in actives]),
        'active_range': np.array([a.max() - a.min() if len(a) > 1 else 0.0 for a in actives]),
        'active_skew':  np.array([skew(a) if len(a) > 2 else 0.0 for a in actives]),
    }
    print(f"\n=== {label}: {M.shape[0]} drugs ===")
    print(f"{'feature':13s} {'family':8s} {'placeholder 0 (paper)':>26s} {'defined only':>26s}")
    for col, k in FEATURES.items():
        all_active = n_active >= 1
        defined = n_active >= k
        for fam in families:
            r0, p0 = spearmanr(feats[col][all_active], instability[fam][all_active])
            r1, p1 = spearmanr(feats[col][defined], instability[fam][defined])
            print(f"{col:13s} {fam:8s} {f'r={r0:+.3f} p={p0:.1e} n={all_active.sum()}':>26s} "
                  f"{f'r={r1:+.3f} p={p1:.1e} n={defined.sum()}':>26s}")


run(pd.read_csv("klaeger_matrix.csv", index_col=0).values, "Klaeger")
davis = pd.read_csv("davis_affinity.csv")
run(davis.pivot(index='Drug_Index', columns='Protein_Index', values='Affinity').fillna(5.0).values, "Davis")
