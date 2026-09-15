"""
Two more candidate-family variants beyond the two axes in family_construction.py
(Renyi order, quantile floor). Both are checked against D1-D4 and against the
alpha=1 candidate, and both are negative results kept here for the record rather
than proposed as usable alternatives:

  1. Soft gate. Replaces the hard indicator 1[x_i > floor] with a sigmoid ramp,
     to see whether D4 (monotonicity) and correct orientation survive losing the
     hard cutoff that candidate_measure.py's docstring says is essential.
  2. Rank-based weights. w_i depends only on the rank position of an active
     kinase, not its affinity value, so two profiles with the same active-set
     ranks score identically regardless of the gap between them. This is a
     deliberate stress test of D2: D2 requires bounded gap sensitivity, not zero
     gap sensitivity, and a measure that discards gap size entirely is a
     degenerate satisfier of the property as stated.
"""
import script_logging; script_logging.capture(__file__)
import os
import numpy as np
import pandas as pd
from scipy.stats import spearmanr, rankdata

HERE = os.path.dirname(os.path.abspath(__file__))


def to_ranks(s):
    return rankdata(-np.asarray(s), method='average')


def candidate(P, floor, T=1.0, eps=1e-12):
    g = (P > floor).astype(float)
    w = g * (T * np.logaddexp(0.0, (P - floor) / T))
    tot = np.where(w.sum(1, keepdims=True) <= 0, eps, w.sum(1, keepdims=True))
    p = w / tot
    return (p * np.where(p > 0, np.log2(p + eps), 0.0)).sum(1)


def candidate_softgate(P, floor, T, T_gate, eps=1e-12):
    g = 1.0 / (1.0 + np.exp(-(P - floor) / T_gate))
    w = g * (T * np.logaddexp(0.0, (P - floor) / T))
    tot = np.where(w.sum(1, keepdims=True) <= 0, eps, w.sum(1, keepdims=True))
    p = w / tot
    return (p * np.where(p > 0, np.log2(p + eps), 0.0)).sum(1)


def candidate_rankweight(P, floor, decay):
    out = np.zeros(P.shape[0])
    for i, row in enumerate(P):
        active_vals = row[row > floor]
        if len(active_vals) == 0:
            out[i] = 0.0
            continue
        w = decay ** np.arange(len(active_vals))
        p = w / w.sum()
        out[i] = (p * np.log2(p + 1e-12)).sum()
    return out


def worst_pair(fn, grid):
    rr = np.array([to_ranks(fn(v)) for v in grid])
    return min(spearmanr(rr[i], rr[j])[0]
               for i in range(len(grid)) for j in range(i + 1, len(grid)))


def d4_check(fn, M, sub_value):
    base = fn(M)
    new = fn(np.hstack([M, np.full((M.shape[0], 1), sub_value)]))
    return (new <= base + 1e-9).mean() * 100, max(0.0, (new - base).max())


def load():
    aff = pd.read_csv(f"{HERE}/davis_affinity.csv")
    davis = aff.pivot(index="Drug_Index", columns="Protein_Index", values="Affinity").fillna(5.0).values
    return {
        "Davis":   (davis, 5.0, 6.0),
        "Klaeger": (pd.read_csv(f"{HERE}/klaeger_matrix.csv", index_col=0).values, 5.0, 6.0),
    }


datasets = load()
BSWEEP = np.arange(4.5, 6.6, 0.25)

print("=" * 70)
print("1. Soft gate (sigmoid ramp instead of the hard floor indicator)")
print("=" * 70)
TGATES = [0.05, 0.1, 0.25, 0.5, 1.0]
for name, (M, floor, active) in datasets.items():
    nact = (M > active).sum(1)
    sub = floor + 0.5 * (active - floor)
    fixed_rank = to_ranks(candidate(M, floor))
    print(f"\n{name}:")
    for tg in TGATES:
        sel = candidate_softgate(M, floor, 1.0, tg)
        orient = spearmanr(sel, nact)[0]
        d3 = worst_pair(lambda b: candidate_softgate(M, b, 1.0, tg), BSWEEP)
        pct, maxinc = d4_check(lambda P: candidate_softgate(P, floor, 1.0, tg), M, sub)
        agree = spearmanr(to_ranks(sel), fixed_rank)[0]
        print(f"  T_gate={tg:<5} orient_r={orient:+.3f}  D3(worst-case)={d3:+.3f}  "
              f"D4(pct-non-incr)={pct:5.1f} maxinc={maxinc:.4f}  vs-hard-gate={agree:.3f}")

print()
print("=" * 70)
print("2. Rank-based weights (D2 stress test: same active-set ranks, different gaps)")
print("=" * 70)
DECAYS = [0.5, 0.7, 0.9]
for name, (M, floor, active) in datasets.items():
    nact = (M > active).sum(1)
    fixed_rank = to_ranks(candidate(M, floor))
    print(f"\n{name}:")
    for d in DECAYS:
        sel = candidate_rankweight(M, floor, d)
        orient = spearmanr(sel, nact)[0]
        agree = spearmanr(to_ranks(sel), fixed_rank)[0]
        pct, maxinc = d4_check(lambda P: candidate_rankweight(P, floor, d), M,
                                floor + 0.5 * (active - floor))
        print(f"  decay={d:<4} orient_r={orient:+.3f}  D4(pct-non-incr)={pct:5.1f} "
              f"maxinc={maxinc:.4f}  vs-magnitude-weighted={agree:.3f}")
    row_small_gap = np.full(M.shape[1], floor - 1.0)
    row_small_gap[:3] = [8.0, 7.9, 7.8]
    row_big_gap = np.full(M.shape[1], floor - 1.0)
    row_big_gap[:3] = [8.0, 6.2, 6.1]
    for d in DECAYS:
        s_small = candidate_rankweight(np.array([row_small_gap]), floor, d)[0]
        s_big = candidate_rankweight(np.array([row_big_gap]), floor, d)[0]
        print(f"  decay={d:<4} D2 stress: small-gap score={s_small:.4f}  "
              f"big-gap score={s_big:.4f}  diff={abs(s_small - s_big):.6f}")

print("\nDone.")
