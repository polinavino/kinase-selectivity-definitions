"""
Cross-dataset validation of the candidate selectivity measure.

For each dataset it checks orientation, the pairs reversing when the gate is
swept over B = [x0, tau*] (on the same 201-point grid as d3_pair_reversals.py)
and D3 (monotonicity, one kinase added at tau*), using the same per-dataset floor
and activity threshold used elsewhere in the analysis.  It also (a) sweeps the
scale T to show the ranking is stable across it, and (b) shows that removing the
floor gate inverts the measure, which is what makes the gate necessary.

Writes candidate_validation.csv and prints a summary.  Deterministic.
"""
import script_logging; script_logging.capture(__file__)
import os
import numpy as np
import pandas as pd
from scipy.stats import spearmanr, rankdata

HERE = os.path.dirname(os.path.abspath(__file__))


def to_ranks(s):
    return rankdata(-np.asarray(s), method='average')


def negH(w):
    tot = w.sum(1, keepdims=True)
    tot = np.where(tot <= 0, 1e-12, tot)
    p = w / tot
    return np.where(p > 0, p * np.log2(np.where(p > 0, p, 1.0)), 0.0).sum(1)


def candidate(P, floor, T, gate=True):
    """Gated entropy: exponential weights above the gate (the floor by default)."""
    g = (P > floor).astype(float) if gate else np.ones_like(P)
    return negH(g * np.exp((P - floor) / T))


def worst_pair(fn, grid):
    rr = np.array([to_ranks(fn(b)) for b in grid])
    return min(spearmanr(rr[i], rr[j])[0]
               for i in range(len(grid)) for j in range(i + 1, len(grid)))


def rev_frac(fn, grid, keep):
    """Fraction of pairs of scored compounds whose order reverses somewhere in the
    grid (as in d3_pair_reversals.py)."""
    S = np.array([fn(v) for v in grid])[:, keep]
    G, n = S.shape
    chunk = max(1, int(4e7 // (G * n)))
    count = 0
    for i0 in range(0, n, chunk):
        d = S[:, i0:i0 + chunk, None] - S[:, None, :]
        rev = (d > 0).any(0) & (d < 0).any(0)
        rows = np.arange(i0, min(i0 + chunk, n))[:, None]
        count += int((rev & (np.arange(n)[None, :] > rows)).sum())
    return count / (n * (n - 1) / 2)


# per-dataset: (matrix, floor, activity threshold, T, boundary sweep).
# The sweep is 201 points over B = [x0, tau*] (Methods), the same for the
# candidate's gate and for the entropy's baseline, as in d3_pair_reversals.py.
def load():
    aff = pd.read_csv(f"{HERE}/davis_affinity.csv")
    davis = aff.pivot(index="Drug_Index", columns="Protein_Index", values="Affinity").fillna(5.0).values
    return {
        # T = active threshold - floor in every dataset (1, 1, 30, 2): below that
        # distance the candidate can violate D3 (d4_entropy_vs_candidate.py).
        "Davis":         (davis, 5.0, 6.0, 1.0, np.linspace(5.0, 6.0, 201)),
        "Klaeger":       (pd.read_csv(f"{HERE}/klaeger_matrix.csv", index_col=0).values, 5.0, 6.0, 1.0, np.linspace(5.0, 6.0, 201)),
        # percent-inhibition assay: x0 is an operational weak-inhibition cutoff of 20%
        # (54% of entries lie between 0% and 20%), the analog of the pKd = 5 floor.
        "Anastassiadis": (pd.read_csv(f"{HERE}/anastassiadis_matrix.csv", index_col=0).values, 20.0, 50.0, 30.0, np.linspace(20.0, 50.0, 201)),
        "Metz":          (pd.read_csv(f"{HERE}/metz_matrix.csv", index_col=0).values, 4.0, 6.0, 2.0, np.linspace(4.0, 6.0, 201)),
    }


def entropy_sel(P, b, eps=1e-10):
    """Baseline-subtracted entropy (selectivity orientation), as in the main analysis."""
    sh = np.maximum(P - b, 0)
    rs = np.where(sh.sum(1, keepdims=True) == 0, eps, sh.sum(1, keepdims=True))
    p = sh / rs
    return (p * np.where(p > 0, np.log2(np.where(p > 0, p, 1.0)), 0)).sum(1)


def validate():
    rows = []
    for name, (M, floor, active, T, bsweep) in load().items():
        nact = (M > active).sum(1)
        sel = candidate(M, floor, T)
        keep = nact > 0                                        # compounds D1 allows to be scored
        orient = spearmanr(sel[keep], nact[keep])[0]           # want strongly negative
        # boundary varied: the gate moves, T stays fixed. The score jumps where the
        # gate passes an observed value: evaluate at every observed value in B and
        # just below it, as well as on the grid
        v = np.unique(M[(M > floor) & (M <= active)])
        gsweep = np.unique(np.concatenate([bsweep, v, np.nextafter(v, -np.inf)]))
        d3_bound = rev_frac(lambda b: candidate(M, b, T), gsweep, keep)
        # the entropy's own baseline swept over the same grid, for comparison
        d3_ent = rev_frac(lambda b: entropy_sel(M, b), bsweep, keep)
        sub = active                                           # worst case for D3: an addition exactly at tau*
        base = candidate(M, floor, T)
        new = candidate(np.hstack([M, np.full((M.shape[0], 1), sub)]), floor, T)
        d4 = (new <= base + 1e-9)[keep].mean() * 100           # over the compounds D1 scores
        rows.append(dict(dataset=name, n_cpd=M.shape[0], n_kin=M.shape[1],
                         orient_r=round(orient, 3),
                         gate_rev=round(d3_bound, 4), entropy_rev=round(d3_ent, 4),
                         d4_pct=round(d4, 1)))
    return pd.DataFrame(rows)


def t_sweep():
    M, floor, active, _, bsweep = load()["Klaeger"]
    nact = (M > active).sum(1)
    rows = []
    for T in [0.5, 0.75, 1.0, 1.5, 2.0]:
        sel = candidate(M, floor, T)
        orient = spearmanr(sel[nact > 0], nact[nact > 0])[0]
        added = candidate(np.hstack([M, np.full((M.shape[0], 1), active)]), floor, T)
        d4 = (added <= sel + 1e-9)[nact > 0].mean() * 100   # D3: one kinase added at tau*
        rows.append(dict(T=T, orient_r=round(orient, 3), d4_pct=round(d4, 1)))
    # cross-T ranking agreement (are rankings stable as T varies?)
    Ts = [0.5, 0.75, 1.0, 1.5, 2.0]
    rk = [to_ranks(candidate(M, floor, T)[nact > 0]) for T in Ts]
    worstT = min(spearmanr(rk[i], rk[j])[0] for i in range(len(rk)) for j in range(i + 1, len(rk)))
    ok = [i for i, T in enumerate(Ts) if T >= 1.0]   # the scales for which D3 is guaranteed
    worstT_ok = min(spearmanr(rk[i], rk[j])[0] for i in ok for j in ok if j > i)
    return pd.DataFrame(rows), round(worstT, 3), round(worstT_ok, 3)


def inversion_demo():
    M, floor, active, T, _ = load()["Klaeger"]
    M = M[(M > active).sum(1) > 0]                             # compounds D1 allows to be scored
    nact = (M > active).sum(1)
    gated = candidate(M, floor, T)
    ungated = candidate(M, floor, T, gate=False)
    # distribution-family reference ranking (S-score, higher=selective)
    fam = -(M > active).astype(float).mean(1)
    return dict(gated_vs_nactive=round(spearmanr(gated, nact)[0], 3),
                ungated_vs_nactive=round(spearmanr(ungated, nact)[0], 3),
                ungated_vs_family=round(spearmanr(ungated, fam)[0], 3),
                gated_vs_family=round(spearmanr(gated, fam)[0], 3))


if __name__ == "__main__":
    df = validate()
    df.to_csv(f"{HERE}/candidate_validation.csv", index=False)
    print("Cross-dataset validation (orientation must be negative; *_rev = fraction of scored-compound pairs reversing in B; D3 (d4_pct) = 100):")
    print(df.to_string(index=False))
    ts, worstT, worstT_ok = t_sweep()
    print("\nScale T sweep on Klaeger (D3 is guaranteed for T >= 1):")
    print(ts.to_string(index=False))
    print(f"worst-case cross-T ranking agreement (Spearman) = {worstT}; over T >= 1 = {worstT_ok}")
    inv = inversion_demo()
    print("\nWhy the gate is needed (Klaeger):")
    print(f"  gated   selectivity vs n_active = {inv['gated_vs_nactive']:+}  (correct: negative)")
    print(f"  ungated selectivity vs n_active = {inv['ungated_vs_nactive']:+}  (inverted: positive)")
    print(f"  ungated selectivity vs distribution family = {inv['ungated_vs_family']:+}  (anti-correlated)")
    print(f"  gated   selectivity vs distribution family = {inv['gated_vs_family']:+}")
    print("\nSaved candidate_validation.csv")
