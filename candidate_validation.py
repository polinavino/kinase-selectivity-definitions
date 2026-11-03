"""
Cross-dataset validation of the candidate selectivity measure.

For each dataset it checks the three properties that do not depend on panel
subsampling (orientation, D3, D4), using the same per-dataset detection floor and
active threshold used elsewhere in the analysis.  It also (a) sweeps the smoothing
width T to show the ranking is stable across it, and (b) shows that removing the
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


def candidate(P, floor, T, baseline=None, gate=True):
    """Gated smooth-hinge effective-target number. baseline defaults to floor."""
    b = floor if baseline is None else baseline
    g = (P > floor).astype(float) if gate else np.ones_like(P)
    return negH(g * (T * np.logaddexp(0.0, (P - b) / T)))


def worst_pair(fn, grid):
    rr = np.array([to_ranks(fn(b)) for b in grid])
    return min(spearmanr(rr[i], rr[j])[0]
               for i in range(len(grid)) for j in range(i + 1, len(grid)))


def rev_frac(fn, grid, keep):
    """D3 as stated in Methods: fraction of pairs of scored compounds whose order
    reverses somewhere in the grid (see d3_pair_reversals.py)."""
    S = np.array([fn(v) for v in grid])[:, keep]
    iu = np.triu_indices(S.shape[1], 1)
    sign = np.sign(S[:, iu[0]] - S[:, iu[1]])
    return ((sign.max(0) > 0) & (sign.min(0) < 0)).mean()


# per-dataset: (matrix, detection floor, active threshold, T, D3 sweep).
# The D3 sweep is B = [x0, tau*] (Methods), the same for the candidate's gate and
# emphasis and for the entropy's baseline, as in d3_pair_reversals.py.
def load():
    aff = pd.read_csv(f"{HERE}/davis_affinity.csv")
    davis = aff.pivot(index="Drug_Index", columns="Protein_Index", values="Affinity").fillna(5.0).values
    return {
        # T = active threshold - floor in every dataset (1, 1, 30, 2): below 0.58 x that
        # distance the candidate can violate D4 (d4_entropy_vs_candidate.py).
        "Davis":         (davis, 5.0, 6.0, 1.0, np.arange(5.0, 6.01, 0.25)),
        "Klaeger":       (pd.read_csv(f"{HERE}/klaeger_matrix.csv", index_col=0).values, 5.0, 6.0, 1.0, np.arange(5.0, 6.01, 0.25)),
        # percent-inhibition assay: detection floor set to the 20% noise level
        # (median value is ~4%), the analog of the pKd = 5 background floor.
        "Anastassiadis": (pd.read_csv(f"{HERE}/anastassiadis_matrix.csv", index_col=0).values, 20.0, 50.0, 30.0, np.arange(20.0, 50.1, 5.0)),
        "Metz":          (pd.read_csv(f"{HERE}/metz_matrix.csv", index_col=0).values, 4.0, 6.0, 2.0, np.arange(4.0, 6.01, 0.25)),
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
        d3_emph = rev_frac(lambda b: candidate(M, floor, T, baseline=b), bsweep, keep)  # gate pinned
        # boundary varied: gate and emphasis move together (floor = b)
        d3_bound = rev_frac(lambda b: candidate(M, b, T, baseline=b), bsweep, keep)
        # the entropy's own baseline swept over the same grid, for comparison
        d3_ent = rev_frac(lambda b: entropy_sel(M, b), bsweep, keep)
        sub = active                                           # worst case for D4: an addition exactly at tau*
        base = candidate(M, floor, T)
        new = candidate(np.hstack([M, np.full((M.shape[0], 1), sub)]), floor, T)
        d4 = (new <= base + 1e-9).mean() * 100
        rows.append(dict(dataset=name, n_cpd=M.shape[0], n_kin=M.shape[1],
                         orient_r=round(orient, 3), d3_emphasis_rev=round(d3_emph, 4),
                         d3_boundary_rev=round(d3_bound, 4), d3_entropy_rev=round(d3_ent, 4),
                         d4_pct=round(d4, 1)))
    return pd.DataFrame(rows)


def t_sweep():
    M, floor, active, _, bsweep = load()["Klaeger"]
    nact = (M > active).sum(1)
    rows = []
    for T in [0.5, 0.75, 1.0, 1.5, 2.0]:
        sel = candidate(M, floor, T)
        orient = spearmanr(sel[nact > 0], nact[nact > 0])[0]
        d3 = rev_frac(lambda b: candidate(M, floor, T, baseline=b), bsweep, nact > 0)
        rows.append(dict(T=T, orient_r=round(orient, 3), d3_emphasis_rev=round(d3, 4)))
    # cross-T ranking agreement (are rankings stable as T varies?)
    rk = [to_ranks(candidate(M, floor, T)[nact > 0]) for T in [0.5, 0.75, 1.0, 1.5, 2.0]]
    worstT = min(spearmanr(rk[i], rk[j])[0] for i in range(len(rk)) for j in range(i + 1, len(rk)))
    return pd.DataFrame(rows), round(worstT, 3)


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
    print("Cross-dataset validation (orientation must be negative; D3 = fraction of scored-compound pairs reversing in B, emphasis must be <= 0.01; D4 = 100):")
    print(df.to_string(index=False))
    ts, worstT = t_sweep()
    print("\nSmoothing-width T sweep on Klaeger:")
    print(ts.to_string(index=False))
    print(f"worst-case cross-T ranking agreement (Spearman) = {worstT}")
    inv = inversion_demo()
    print("\nWhy the gate is needed (Klaeger):")
    print(f"  gated   selectivity vs n_active = {inv['gated_vs_nactive']:+}  (correct: negative)")
    print(f"  ungated selectivity vs n_active = {inv['ungated_vs_nactive']:+}  (inverted: positive)")
    print(f"  ungated selectivity vs distribution family = {inv['ungated_vs_family']:+}  (anti-correlated)")
    print(f"  gated   selectivity vs distribution family = {inv['gated_vs_family']:+}")
    print("\nSaved candidate_validation.csv")
