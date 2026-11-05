"""
Variants of the candidate: which other measures of the same construction meet D1 and D3?

This script builds two generalizations of the candidate in candidate_measure.py /
candidate_validation.py and checks D1-D3 (and cross-dataset ranking agreement with
the original candidate) empirically, the same way candidate_validation.py does.

  1. Renyi order alpha.  The candidate's score is -H(p), the negative Shannon
     entropy of the gated, smoothed profile.  Shannon is the alpha -> 1 case of the
     Renyi family H_alpha(p) = 1/(1-alpha) * log2(sum_i p_i^alpha).  For any alpha > 0,
     p_i = 0 still contributes 0 to sum_i p_i^alpha (0^alpha = 0), so additions at
     or below the floor leave the score unchanged; D3 (additions in (x0, tau*]) is
     proved only for alpha = 1 and is checked here with one kinase added at tau*.

  2. Quantile-referenced floor: the fixed assay floor replaced by a per-compound
     low quantile of the profile (SI Section S4). The script also counts the
     compounds whose quantile equals the fixed floor, since on these data the
     variant then coincides with the featured measure.
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
    undefined = tot.ravel() <= 0
    tot = np.where(tot <= 0, eps, tot)
    p = w / tot
    if abs(alpha - 1.0) < 1e-9:
        out = (p * np.where(p > 0, np.log2(np.where(p > 0, p, 1.0)), 0.0)).sum(1)   # = -H(p)
    else:
        s = np.where(p > 0, p ** alpha, 0.0).sum(1)
        s = np.maximum(s, eps)
        out = (1.0 / (alpha - 1.0)) * np.log2(s)                     # = -H_alpha(p)
    # Undefined (no gated weight) scores 0, tied with a single-target profile.
    return np.where(undefined, 0.0, out)


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


def rev_frac(fn, grid, keep):
    """D3 as stated in Methods: fraction of pairs of scored compounds whose order
    reverses somewhere in the grid (see d3_pair_reversals.py)."""
    S = np.array([fn(v) for v in grid])[:, keep]
    iu = np.triu_indices(S.shape[1], 1)
    sign = np.sign(S[:, iu[0]] - S[:, iu[1]])
    return ((sign.max(0) > 0) & (sign.min(0) < 0)).mean()


def load():
    aff = pd.read_csv(f"{HERE}/davis_affinity.csv")
    davis = aff.pivot(index="Drug_Index", columns="Protein_Index", values="Affinity").fillna(5.0).values
    return {
        "Davis":   (davis, 5.0, 6.0, 1.0),
        "Klaeger": (pd.read_csv(f"{HERE}/klaeger_matrix.csv", index_col=0).values, 5.0, 6.0, 1.0),
    }


datasets = load()
BSWEEP = np.linspace(5.0, 6.0, 201)   # B = [x0, tau*] (Methods, D2), 201 points as in d3_pair_reversals.py
QSWEEP = [0.05, 0.075, 0.10, 0.125, 0.15]
ALPHAS = [0.5, 1.0, 2.0, 3.0]

print("=" * 70)
print("Renyi-order family: alpha=1 is the Shannon candidate")
print("=" * 70)
for name, (M, floor, active, T) in datasets.items():
    nact = (M > active).sum(1)
    base_rank = to_ranks(candidate_alpha(M, floor, T, alpha=1.0))
    print(f"\n{name} ({M.shape[0]} cpd x {M.shape[1]} kin):")
    sub = active   # worst case for D3: an addition exactly at tau*
    for alpha in ALPHAS:
        sel = candidate_alpha(M, floor, T, alpha=alpha)
        orient = spearmanr(sel, nact)[0]
        d3 = rev_frac(lambda b: candidate_alpha(M, floor, T, alpha=alpha, baseline=b), BSWEEP, nact > 0)
        base = sel
        new = candidate_alpha(np.hstack([M, np.full((M.shape[0], 1), sub)]), floor, T, alpha=alpha)
        d4 = (new <= base + 1e-9).mean() * 100
        agree = spearmanr(to_ranks(sel), base_rank)[0]
        print(f"  alpha={alpha:<4} orient_r={orient:+.3f}  D2(pairs reversing in B)={100 * d3:.2f}%  "
              f"D3(% non-increasing)={d4:5.1f}  Spearman-vs-alpha=1={agree:.3f}")

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
    d3_q_rev = rev_frac(lambda q: candidate_quantile(M, q, T), QSWEEP, nact > 0)
    sub = active   # worst case for D3: an addition exactly at tau*
    base_q = sel_q
    new_q = candidate_quantile(np.hstack([M, np.full((M.shape[0], 1), sub)]), q_ref, T)
    d4_q = (new_q <= base_q + 1e-9).mean() * 100
    max_increase = max(0.0, (new_q - base_q).max())
    at_floor = {q: int(np.isclose(np.quantile(M, q, axis=1), floor).sum()) for q in QSWEEP}
    print(f"  compounds whose q-th quantile equals the fixed floor (of {M.shape[0]}): "
          + ", ".join(f"q={q}: {n}" for q, n in at_floor.items()))
    print(f"  q={q_ref}: orient_r={orient_q:+.3f}  "
          f"Spearman vs fixed-floor candidate = {agree:.3f}")
    print(f"  D2 over q in {QSWEEP}: pairs reversing = {100 * d3_q_rev:.2f}%, worst-case rank Spearman = {d3_q:+.3f}")
    print(f"  D3 (add a sub-threshold off-target): non-increasing for {d4_q:.1f}% of "
          f"compounds; max increase = {max_increase:.4f}")

print("\nDone.")
