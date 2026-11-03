"""
D3 (distributional consistency) verdicts for every measure on all four datasets.

D3 as stated in Methods: over the interval B = [x0, tau*] between the assay's
detection floor and the reliability threshold, at most a fraction eps of compound
pairs may reverse order. A pair reverses in B if f(x; b) - f(y; b) is positive at
some b in B and negative at another. Pairs are taken among the compounds D1 allows
to be scored (at least one kinase above tau*), so that D3 is assessed separately
from D1.

Parameters swept over B: the S-score's threshold tau, the entropy's and Gini's
baseline beta, and, for the candidate, either its emphasis baseline beta (gate
fixed at x0, the featured measure) or its gate x0 itself (beta pinned to x0, the
like-for-like comparison with the entropy's baseline). The ratio has no such
parameter and is not swept.

The worst-case pairwise Spearman correlation between the rankings at any two
values in B is printed alongside, for comparison with earlier versions.
"""
import script_logging; script_logging.capture(__file__)
import os
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

HERE = os.path.dirname(os.path.abspath(__file__))
EPS = 0.01


def neg_entropy(w):
    tot = w.sum(1, keepdims=True)
    p = w / np.where(tot <= 0, 1.0, tot)
    return np.where(p > 0, p * np.log2(np.where(p > 0, p, 1.0)), 0.0).sum(1)


def s_score(M, t):
    return -(M > t).mean(1)


def entropy(M, b):
    return neg_entropy(np.maximum(M - b, 0.0))


def gini(M, b):
    s = np.sort(np.maximum(M - b, 0.0), 1)
    n = s.shape[1]
    t = s.sum(1)
    g = (2 * (np.arange(1, n + 1) * s).sum(1)) / (n * np.where(t == 0, 1, t)) - (n + 1) / n
    return np.where(t == 0, 0.0, g)


def candidate(M, x0, T, b):
    return neg_entropy(np.where(M > x0, T * np.logaddexp(0.0, (M - b) / T), 0.0))


def reversals(score_rows, keep):
    """Fraction of compound pairs whose order reverses somewhere in B, and the
    worst-case pairwise Spearman correlation between the rankings in B."""
    S = np.array(score_rows)[:, keep]
    iu = np.triu_indices(S.shape[1], 1)
    sign = np.sign(S[:, iu[0]] - S[:, iu[1]])
    frac = ((sign.max(0) > 0) & (sign.min(0) < 0)).mean()
    worst = min(spearmanr(S[i], S[j])[0] for i in range(len(S)) for j in range(i + 1, len(S)))
    return frac, worst


aff = pd.read_csv(f"{HERE}/davis_affinity.csv")
DATA = {  # matrix, floor x0, reliability threshold tau*, grid step
    "Davis":         (aff.pivot(index="Drug_Index", columns="Protein_Index", values="Affinity").fillna(5.0).values,
                      5.0, 6.0, 0.25),
    "Klaeger":       (pd.read_csv(f"{HERE}/klaeger_matrix.csv", index_col=0).values, 5.0, 6.0, 0.25),
    "Anastassiadis": (pd.read_csv(f"{HERE}/anastassiadis_matrix.csv", index_col=0).values, 20.0, 50.0, 5.0),
    "Metz":          (pd.read_csv(f"{HERE}/metz_matrix.csv", index_col=0).values, 4.0, 6.0, 0.25),
}

print(f"D3 over B = [x0, tau*]: fraction of pairs reversing (rev) and worst-case pairwise Spearman (rho); "
      f"pass if rev <= {EPS:.0%}")
passing, failing = [], []
for name, (M, x0, tau, step) in DATA.items():
    B = np.arange(x0, tau + step / 2, step)
    keep = (M > tau).any(1)
    T = tau - x0
    rows = {
        "S-score (tau)":            [s_score(M, b) for b in B],
        "entropy (beta)":           [entropy(M, b) for b in B],
        "Gini (beta)":              [gini(M, b) for b in B],
        "candidate (emphasis beta)": [candidate(M, x0, T, b) for b in B],
        "candidate (gate x0)":      [candidate(M, b, T, b) for b in B],
    }
    print(f"\n{name}: B = [{x0:g}, {tau:g}] ({len(B)} values), {keep.sum()} scored compounds")
    for label, r in rows.items():
        frac, worst = reversals(r, keep)
        (passing if frac <= EPS else failing).append(frac)
        print(f"  {label:27s} rev = {100 * frac:6.2f}%   rho = {worst:.3f}   {'PASS' if frac <= EPS else 'fail'}")

print(f"\nLargest passing fraction: {100 * max(passing):.2f}%; smallest failing fraction: "
      f"{100 * min(failing):.2f}%. Every verdict is the same for any eps in between.")
