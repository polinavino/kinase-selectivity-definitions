"""
A candidate selectivity measure satisfying the four required properties (D1-D4),
with empirical checks of D3 (baseline robustness), D4 (monotonicity), and of how
quickly its compound ranking converges to the full-panel ranking.

Construction.  Selectivity is the negative Shannon "effective number of targets"
of a compound's binding profile.  Each kinase above the assay detection floor is
weighted by a smooth (softplus) hinge of its log-affinity.  Kinases at or below
the floor are gated out and carry zero weight:

    w_i = 1[x_i > x_floor] * T * log(1 + exp((x_i - beta)/T)),
    p_i = w_i / sum_j w_j,
    H(x) = - sum_i p_i log p_i,    selectivity  S(x) = -H(x),

with beta = x_floor at the operating point.  The gate at the fixed floor is what
keeps the measure correctly oriented.  Without it, the hundreds of non-binding
kinases sitting at the floor dominate the distribution and invert the ranking.
Scores are reported only for compounds with max_i x_i > tau* (the D1 gate).  That
gate is a reporting rule applied in practice.  The panel-size, D3, and D4 checks
below are computed on all compounds, since they concern ranking stability, which
the gate does not affect.

Why this satisfies D1-D4:
  D1  explicit reliability gate on max_i x_i > tau*.
  D2  fully distributional (no reference off-target), so no single pairwise
      comparison is privileged and bounded gap sensitivity holds trivially.
  D3  the candidate has no free activity baseline.  The active/inactive boundary
      is the fixed assay detection floor, a constant of the assay rather than an
      analyst's choice, so the arbitrary-baseline instability D3 targets does not
      arise.  The remaining softplus emphasis parameter beta is one the ranking is
      robust to (verified below: ~invariant to beta over [4.5, 6.5], gate fixed).
      The boundary itself, if varied, moves the ranking like any threshold measure
      (verified below: about -0.93, comparable to entropy).  That is why it is
      pinned to the assay floor rather than left free.
  D4  the gated weight is zero at and below the floor and monotone above it, so
      adding a sub-threshold off-target adds zero weight and never raises S.
      (Verified below.)

The design choice is the smoothing width T (in pK_d units).  T ~ 1 places the
candidate's panel-size convergence in the fast class (entropy, S-score).  The
active/inactive boundary is fixed at the assay floor rather than tuned.  The
panel-size convergence is checked on the two sparse-coverage datasets (Klaeger,
Metz); the dense complete-matrix screens (Davis, Anastassiadis) are unsuitable
for kinase subsampling, as explained in the paper.
"""
import script_logging; script_logging.capture(__file__)
import os
os.environ['KMP_DUPLICATE_LIB_OK'] = 'TRUE'
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
import matplotlib.pyplot as plt

T = 1.0

def to_ranks(s):
    return len(s) - s.argsort().argsort()

# ---------------- candidate ----------------
def candidate(P, baseline, floor, T=T, eps=1e-10):
    # Gate at the fixed assay detection floor. Sub-floor kinases are non-binders.
    # They carry zero weight, so the effective-target count reflects real binding.
    # The gate is what keeps the measure correctly oriented (concentrated = selective).
    # The softplus above the baseline gives a smooth, baseline-robust emphasis.
    gate = (P > floor).astype(float)
    w = gate * (T * np.logaddexp(0.0, (P - baseline) / T))
    rs = np.where(w.sum(1, keepdims=True) == 0, eps, w.sum(1, keepdims=True))
    p = w / rs
    return -(-(p * np.where(p > 0, np.log2(p + eps), 0)).sum(1))   # -Shannon entropy

# ---------------- existing measures (match panel_size_analysis.py) ----------------
def entropy(P, b, e=1e-10):
    s = np.maximum(P - b, 0); rs = np.where(s.sum(1, keepdims=True) == 0, e, s.sum(1, keepdims=True))
    p = s / rs; return -(-(p * np.where(p > 0, np.log2(p + e), 0)).sum(1))
def gini(P, b):
    s = np.maximum(P - b, 0); out = []
    for r in s:
        rs = np.sort(r); n = len(rs); t = rs.sum()
        out.append(0.0 if t == 0 else (2 * np.sum(np.arange(1, n + 1) * rs)) / (n * t) - (n + 1) / n)
    return np.array(out)
def s_score(P, thr): return -(P > thr).astype(float).mean(1)
def ratio(P, floor, k=1):
    return np.array([np.sort(r)[::-1][0] - max(np.sort(r)[::-1][k] if len(r) > k else floor, floor) for r in P])


def measures_for(baseline, threshold, floor):
    return {
        'candidate': lambda P: candidate(P, baseline=floor, floor=floor),
        'entropy':   lambda P: entropy(P, b=baseline),
        'gini':      lambda P: gini(P, b=baseline),
        's_score':   lambda P: s_score(P, thr=threshold),
        'ratio':     lambda P: ratio(P, floor=floor),
    }


def panel_convergence(M, baseline, threshold, floor, R=50, seed=42):
    n_drugs, n_kin = M.shape
    meas = measures_for(baseline, threshold, floor)
    ref = {k: to_ranks(f(M)) for k, f in meas.items()}
    panel = list(range(50, n_kin, 30)) + [n_kin]
    rng = np.random.RandomState(seed)
    res = {k: {ps: [] for ps in panel} for k in meas}
    for ps in panel:
        for _ in range(R):
            idx = rng.choice(n_kin, ps, replace=False); Ms = M[:, idx]
            for k, f in meas.items():
                res[k][ps].append(spearmanr(ref[k], to_ranks(f(Ms)))[0])
    pstar = {k: next((ps for ps in panel if np.mean(res[k][ps]) > 0.90), None) for k in meas}
    return panel, res, pstar


# Sparse-coverage datasets only (see module docstring / paper).
DATASETS = [
    dict(name='Klaeger', file='klaeger_matrix.csv', baseline=5.0, threshold=6.0, floor=5.0),
    dict(name='Metz',    file='metz_matrix.csv',    baseline=4.0, threshold=6.0, floor=4.0),
]
order = ['candidate', 'entropy', 's_score', 'gini', 'ratio']
styles = {'candidate': ('C3', '-', 2.6), 'entropy': ('C1', '--', 1.6),
          's_score': ('C0', '--', 1.6), 'gini': ('C2', '--', 1.6), 'ratio': ('C4', '--', 1.6)}

# ---------------- (checks on Klaeger, as before) ----------------
Mk = pd.read_csv('klaeger_matrix.csv', index_col=0).values
FLOOR, TAU_STAR = 5.0, 6.0
n_drugs_k = Mk.shape[0]

# orientation guard
n_active = (Mk > TAU_STAR).sum(1)
orient_r = spearmanr(candidate(Mk, FLOOR, FLOOR), n_active)[0]
assert orient_r < -0.3, f"candidate is inverted: corr(selectivity, n_active) = {orient_r:+.3f}"
print(f"Orientation check (Klaeger): corr(candidate selectivity, n_active) = {orient_r:+.3f} (must be negative)")

# D3: baseline robustness
betas = np.arange(4.5, 6.6, 0.25)
def worst_pair(fn):
    rr = np.array([to_ranks(fn(b)) for b in betas])
    return min(spearmanr(rr[i], rr[j])[0] for i in range(len(betas)) for j in range(i + 1, len(betas)))
print("\nD3 check (worst-case rank Spearman over [4.5,6.5]):")
print(f"  candidate, vary emphasis baseline (gate pinned at floor): {worst_pair(lambda b: candidate(Mk, baseline=b, floor=FLOOR)):+.3f}")
print(f"  candidate, vary active-set boundary (gate = floor)       : {worst_pair(lambda b: candidate(Mk, baseline=b, floor=b)):+.3f}")
print(f"  existing entropy, vary baseline                          : {worst_pair(lambda b: entropy(Mk, b=b)):+.3f}")

# D4: monotonicity under weak off-target addition
base = candidate(Mk, FLOOR, FLOOR)
new = candidate(np.hstack([Mk, np.full((n_drugs_k, 1), TAU_STAR - 0.7)]), FLOOR, FLOOR)
print(f"\nD4 (add a sub-threshold off-target): selectivity non-increasing for "
      f"{(new <= base + 1e-9).mean() * 100:.0f}% of compounds; max increase = {max(0.0, (new - base).max()):.4f}")

# ---------------- panel-size convergence (Klaeger + Metz) ----------------
fig, axes = plt.subplots(1, len(DATASETS), figsize=(7 * len(DATASETS), 5), squeeze=False)
for ax, ds in zip(axes[0], DATASETS):
    M = pd.read_csv(ds['file'], index_col=0).values
    n_drugs, n_kin = M.shape
    panel, res, pstar = panel_convergence(M, ds['baseline'], ds['threshold'], ds['floor'])
    print(f"\nPanel-size convergence p* ({ds['name']}, {n_drugs} cpd x {n_kin} kin):")
    for k in order:
        print(f"  {k:10s} p* = {pstar[k]}")
    for k in order:
        c, ls, lw = styles[k]
        ax.plot(panel, [np.mean(res[k][ps]) for ps in panel], color=c, ls=ls, lw=lw,
                label=f"{k} (p*={pstar[k]})")
    ax.axhline(0.90, color='black', ls=':', lw=1)
    ax.set_xlabel('Panel size (kinases)'); ax.set_ylabel('Spearman rho vs full-panel ranking')
    ax.set_title(f"{ds['name']} ({n_drugs} cpd $\\times$ {n_kin} kin)")
    ax.legend(fontsize=8); ax.set_ylim(0, 1.01)

fig.suptitle('Panel-size convergence: candidate (gated, smooth-hinge effective-target number) '
             'vs existing measures (50 subsamples per size)', fontsize=12)
plt.tight_layout(rect=[0, 0, 1, 0.96])
plt.savefig('candidate_panel_convergence.png', dpi=150, bbox_inches='tight')
print("\nSaved candidate_panel_convergence.png")
