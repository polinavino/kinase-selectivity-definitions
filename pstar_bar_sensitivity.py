"""Sensitivity of panel-size recoverability p* to the r/rho > 0.90 bar.

The main panel-size analysis (panel_size_analysis.py) calls a subpanel size
"recovering" the full-panel ranking once mean Spearman r exceeds 0.90, a
convention borrowed from reliability research's "excellent" agreement bar
(Koo & Li 2016), not a property derived from the data. This script checks
whether the qualitative ordering of the four definitions by panel-size
requirement -- entropy and S-score fast, Gini dataset-dependent, ratio
slowest -- is an artifact of that specific choice, by recomputing p* at
r = 0.80, 0.85, 0.90, 0.95.
"""
import script_logging; script_logging.capture(__file__)
import os
os.environ['KMP_DUPLICATE_LIB_OK'] = 'TRUE'

import numpy as np
import pandas as pd
from scipy.stats import spearmanr, rankdata


def selectivity_entropy(profiles, baseline, epsilon=1e-10):
    shifted = np.maximum(profiles - baseline, 0)
    row_sums = np.where(shifted.sum(axis=1, keepdims=True) == 0, epsilon,
                         shifted.sum(axis=1, keepdims=True))
    p = shifted / row_sums
    return -(p * np.where(p > 0, np.log2(p + epsilon), 0)).sum(axis=1)


def gini_selectivity(profiles, baseline):
    shifted = np.maximum(profiles - baseline, 0)
    ginis = []
    for row in shifted:
        row_sorted = np.sort(row)
        n = len(row_sorted)
        total = row_sorted.sum()
        if total == 0:
            ginis.append(0.0); continue
        ginis.append((2 * np.sum((np.arange(1, n + 1) * row_sorted))) / (n * total) - (n + 1) / n)
    return np.array(ginis)


def s_score(profiles, threshold):
    return -(profiles > threshold).astype(float).mean(axis=1)


def ratio_selectivity(profiles, floor, top_n=1):
    ratios = []
    for row in profiles:
        s = np.sort(row)[::-1]
        ratios.append(s[0] - (max(s[top_n], floor) if len(s) > top_n else floor))
    return np.array(ratios)


def to_ranks(scores):
    return rankdata(-np.asarray(scores), method='average')


def scores_at(M, baseline, threshold, floor):
    return {
        'entropy': to_ranks(selectivity_entropy(M, baseline)),
        'gini':    to_ranks(gini_selectivity(M, baseline)),
        's_score': to_ranks(s_score(M, threshold)),
        'ratio':   to_ranks(ratio_selectivity(M, floor)),
    }


def panel_curves(M, baseline, threshold, floor, n_repeats=50, seed=42):
    n_drugs, n_kinases = M.shape
    ref_ranks = scores_at(M, baseline, threshold, floor)
    panel_sizes = list(range(50, n_kinases, 30)) + [n_kinases]
    rng = np.random.RandomState(seed)
    results = {d: {ps: [] for ps in panel_sizes} for d in ref_ranks}
    for ps in panel_sizes:
        for _ in range(n_repeats):
            idx = rng.choice(n_kinases, ps, replace=False)
            sub = scores_at(M[:, idx], baseline, threshold, floor)
            for d in results:
                results[d][ps].append(spearmanr(ref_ranks[d], sub[d])[0])
    return panel_sizes, results


def pstar_at(panel_sizes, results, d, bar):
    return next((ps for ps in panel_sizes if np.mean(results[d][ps]) > bar), None)


DATASETS = [
    dict(name='Klaeger', file='klaeger_matrix.csv', baseline=5.0, threshold=6.0, floor=5.0),
    dict(name='Metz',    file='metz_matrix.csv',    baseline=4.0, threshold=6.0, floor=4.0),
]
BARS = [0.80, 0.85, 0.90, 0.95]
order = ['entropy', 'gini', 's_score', 'ratio']

for ds in DATASETS:
    M = pd.read_csv(ds['file'], index_col=0).values
    panel_sizes, results = panel_curves(M, ds['baseline'], ds['threshold'], ds['floor'])
    print(f"\n=== {ds['name']} ({M.shape[0]} drugs x {M.shape[1]} kinases) ===")
    print(f"{'bar':>6s}" + "".join(f"{d:>10s}" for d in order))
    for bar in BARS:
        row = f"{bar:>6.2f}"
        for d in order:
            p = pstar_at(panel_sizes, results, d, bar)
            row += f"{(p if p is not None else -1):>10d}"
        print(row)
