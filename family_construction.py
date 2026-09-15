"""
Beyond the single candidate: does D1-D4 pin down a family of measures, or just one?

This script builds two generalizations of the candidate in candidate_measure.py /
candidate_validation.py and checks D1-D4 (and cross-dataset ranking agreement with
the original candidate) empirically, the same way candidate_validation.py does.

  1. Renyi order alpha.  The candidate's score is -H(p), the negative Shannon
     entropy of the gated, smoothed profile.  Shannon is the alpha -> 1 case of the
     Renyi family H_alpha(p) = 1/(1-alpha) * log2(sum_i p_i^alpha).  For any alpha > 0,
     p_i = 0 still contributes 0 to sum_i p_i^alpha (0^alpha = 0), so D4's
     append-a-zero-weight-coordinate argument goes through unchanged, and D2/D3 are
     unaffected by alpha (they concern the p_i map, not the entropy functional
     applied to it). This is a genuine one-parameter family, not just a restatement.

  2. Quantile-referenced floor.  desiderata.tex's closing paragraph claims that
     replacing the fixed assay floor with a per-compound low-quantile reference
     gives a shift-invariant variant with "nearly the same ranking on both
     datasets (Spearman r ~ 1.00)". No script in the repository computed that
     number; this one does, so the claim is either confirmed or corrected.
     Caveat checked here: a per-compound floor is a statistic of the *compound's
     own row*, so appending a column (D4) can in principle move it, unlike the
     fixed floor, which cannot move by construction. Checked empirically below.
"""
import script_logging; script_logging.capture(__file__)
import os
import numpy as np
import pandas as pd
from scipy.stats import spearmanr, rankdata

HERE = os.path.dirname(os.path.abspath(__file__))


def to_ranks(s):
    return rankdata(-np.asarray(s), method='average')


def neg_renyi_H(w, alpha, eps=1e-12):
    """-H_alpha(p): the correctly-oriented candidate score (matches -Shannon-entropy
    at alpha=1, i.e. candidate_measure.py's `candidate`)."""
    tot = w.sum(1, keepdims=True)
    tot = np.where(tot <= 0, eps, tot)
    p = w / tot
    if abs(alpha - 1.0) < 1e-9:
        return (p * np.where(p > 0, np.log2(p + eps), 0.0)).sum(1)   # = -H(p)
    s = np.where(p > 0, p ** alpha, 0.0).sum(1)
    s = np.maximum(s, eps)
    return (1.0 / (alpha - 1.0)) * np.log2(s)                        # = -H_alpha(p)


def candidate_alpha(P, floor, T, alpha, baseline=None):
    b = floor if baseline is None else baseline
    g = (P > floor).astype(float)
    w = g * (T * np.logaddexp(0.0, (P - b) / T))
    return neg_renyi_H(w, alpha)


def candidate_quantile(P, q, T):
    """Per-compound floor at the q-th quantile of that compound's own profile."""
    floor_i = np.quantile(P, q, axis=1, keepdims=True)
    g = (P > floor_i).astype(float)
    w = g * (T * np.logaddexp(0.0, (P - floor_i) / T))
    return neg_renyi_H(w, 1.0)


def worst_pair(fn, grid):
    rr = np.array([to_ranks(fn(v)) for v in grid])
    return min(spearmanr(rr[i], rr[j])[0]
               for i in range(len(grid)) for j in range(i + 1, len(grid)))


def load():
    aff = pd.read_csv(f"{HERE}/davis_affinity.csv")
    davis = aff.pivot(index="Drug_Index", columns="Protein_Index", values="Affinity").fillna(5.0).values
    return {
        "Davis":   (davis, 5.0, 6.0, 1.0),
        "Klaeger": (pd.read_csv(f"{HERE}/klaeger_matrix.csv", index_col=0).values, 5.0, 6.0, 1.0),
    }


datasets = load()
BSWEEP = np.arange(4.5, 6.6, 0.25)
QSWEEP = [0.05, 0.075, 0.10, 0.125, 0.15]
ALPHAS = [0.5, 1.0, 2.0, 3.0]

print("=" * 70)
print("Renyi-order family: alpha=1 is the Shannon candidate")
print("=" * 70)
for name, (M, floor, active, T) in datasets.items():
    nact = (M > active).sum(1)
    base_rank = to_ranks(candidate_alpha(M, floor, T, alpha=1.0))
    print(f"\n{name} ({M.shape[0]} cpd x {M.shape[1]} kin):")
    sub = floor + 0.5 * (active - floor)
    for alpha in ALPHAS:
        sel = candidate_alpha(M, floor, T, alpha=alpha)
        orient = spearmanr(sel, nact)[0]
        d3 = worst_pair(lambda b: candidate_alpha(M, floor, T, alpha=alpha, baseline=b), BSWEEP)
        base = sel
        new = candidate_alpha(np.hstack([M, np.full((M.shape[0], 1), sub)]), floor, T, alpha=alpha)
        d4 = (new <= base + 1e-9).mean() * 100
        agree = spearmanr(to_ranks(sel), base_rank)[0]
        print(f"  alpha={alpha:<4} orient_r={orient:+.3f}  D3(worst-case)={d3:+.3f}  "
              f"D4(%% non-increasing)={d4:5.1f}  Spearman-vs-alpha=1={agree:.3f}")

print()
print("=" * 70)
print("Quantile-referenced floor (per-compound, replaces the fixed assay floor)")
print("=" * 70)
for name, (M, floor, active, T) in datasets.items():
    nact = (M > active).sum(1)
    fixed_rank = to_ranks(candidate_alpha(M, floor, T, alpha=1.0))
    print(f"\n{name} ({M.shape[0]} cpd x {M.shape[1]} kin):")
    q_ref = 0.10
    sel_q = candidate_quantile(M, q_ref, T)
    orient_q = spearmanr(sel_q, nact)[0]
    agree = spearmanr(to_ranks(sel_q), fixed_rank)[0]
    d3_q = worst_pair(lambda q: candidate_quantile(M, q, T), QSWEEP)
    sub = floor + 0.5 * (active - floor)
    base_q = sel_q
    new_q = candidate_quantile(np.hstack([M, np.full((M.shape[0], 1), sub)]), q_ref, T)
    d4_q = (new_q <= base_q + 1e-9).mean() * 100
    max_increase = max(0.0, (new_q - base_q).max())
    print(f"  q={q_ref}: orient_r={orient_q:+.3f}  "
          f"Spearman vs fixed-floor candidate = {agree:.3f}")
    print(f"  D3 (worst-case rank Spearman over q in {QSWEEP}) = {d3_q:+.3f}")
    print(f"  D4 (add a sub-threshold off-target): non-increasing for {d4_q:.1f}% of "
          f"compounds; max increase = {max_increase:.4f}")

print("\nDone.")
