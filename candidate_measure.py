"""
The candidate selectivity measure: its orientation, its agreement with the
entropy, and checks of D2 (boundary sensitivity) and D3 (monotonicity) on
Klaeger. Property numbers follow the paper: D1 activity eligibility, D2
boundary locality, D3 monotonicity under weak off-target addition.

Construction.  Selectivity is the negative Shannon entropy of a compound's
binding profile.  Each kinase above the gate g (the assay floor x0 unless the
gate is swept) is weighted by an exponential of its readout.  Kinases at or
below the gate carry zero weight:

    w_i = 1[x_i > g] * exp((x_i - g)/T),    T = tau* - x0,
    p_i = w_i / sum_j w_j,
    S(x) = sum_i p_i log2 p_i.

The gate keeps the measure correctly oriented: without it, the non-binding
kinases at the floor dominate the distribution and invert the ranking.  Scores
are reported only for compounds with max_i x_i > tau* (the D1 gate).

  D1  explicit eligibility gate on max_i x_i > tau*.
  D2  moving the gate multiplies every retained weight by the same factor, so a
      compound's score changes only when the gate crosses one of its readouts.
      Swept over B, the gate still reverses pairs whose readouts it crosses
      (d3_pair_reversals.py).
  D3  every gated kinase has weight above 1, and a kinase added at or below tau*
      has weight at most e, so the score cannot rise when T >= tau* - x0
      (Methods, Proposition 1; d4_entropy_vs_candidate.py).  The weight jumps
      from 0 to 1 at the floor, so the measure is discontinuous there.

T = tau* - x0 on every dataset, the smallest width for which D3 holds.  Recovery
from random subpanels, compared with every other measure, is in panel_recovery.py.
"""
import script_logging; script_logging.capture(__file__)
import os
os.environ['KMP_DUPLICATE_LIB_OK'] = 'TRUE'
import numpy as np
import pandas as pd
from scipy.stats import spearmanr, rankdata
import matplotlib.pyplot as plt

T = 1.0

def to_ranks(s):
    return rankdata(-np.asarray(s), method='average')

# ---------------- candidate ----------------
def candidate(P, gate, T=T, eps=1e-10):
    # Kinases at or below the gate (the assay detection floor by default) are
    # non-binders and carry zero weight. The gate is what keeps the measure
    # correctly oriented (concentrated = selective). Above it, the weight grows
    # exponentially with the readout, so moving the gate rescales every retained
    # weight by the same factor.
    w = np.where(P > gate, np.exp((P - gate) / T), 0.0)
    rs = np.where(w.sum(1, keepdims=True) == 0, eps, w.sum(1, keepdims=True))
    p = w / rs
    return -(-(p * np.where(p > 0, np.log2(np.where(p > 0, p, 1.0)), 0)).sum(1))   # -Shannon entropy

# ---------------- existing measures (match panel_size_analysis.py) ----------------
def entropy(P, b, e=1e-10):
    s = np.maximum(P - b, 0); rs = np.where(s.sum(1, keepdims=True) == 0, e, s.sum(1, keepdims=True))
    p = s / rs; return -(-(p * np.where(p > 0, np.log2(np.where(p > 0, p, 1.0)), 0)).sum(1))
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
        # T = tau* - floor: the D3 guarantee needs T >= tau* - floor (d4_entropy_vs_candidate.py)
        'candidate': lambda P: candidate(P, gate=floor, T=threshold - floor),
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


def gated_panel_convergence(M, baseline, threshold, floor, tau, R=50, seed=42):
    """Panel-size convergence with the D1 gate actually applied.

    panel_convergence above scores every compound, including ones the candidate
    would refuse to score.  That matters at small panel sizes, because a compound
    whose only active kinases are outside the drawn subpanel falls below tau* and
    becomes undefined.  Here the Spearman is computed only over compounds that are
    D1-defined on BOTH the full panel and the subpanel, which is what a
    practitioner applying the gate would actually do, and the undefined fraction is
    reported alongside.

    Caveat, stated in the paper: the gated and ungated p* are not computed on the
    same population.  Gating removes the least-bound compounds, which are also the
    hardest to rank, so part of any improvement is selection rather than a property
    of the gate.  Both are reported for that reason.
    """
    n_drugs, n_kin = M.shape
    meas = measures_for(baseline, threshold, floor)
    full_def = M.max(1) > tau
    ref = {k: f(M) for k, f in meas.items()}
    panel = list(range(50, n_kin, 30)) + [n_kin]
    rng = np.random.RandomState(seed)
    res = {k: {ps: [] for ps in panel} for k in meas}
    undef = {ps: [] for ps in panel}
    for ps in panel:
        for _ in range(R):
            idx = rng.choice(n_kin, ps, replace=False); Ms = M[:, idx]
            sub_def = Ms.max(1) > tau
            undef[ps].append(1.0 - sub_def.mean())
            keep = full_def & sub_def
            if keep.sum() < 5:
                continue
            for k, f in meas.items():
                res[k][ps].append(spearmanr(to_ranks(ref[k][keep]), to_ranks(f(Ms)[keep]))[0])
    pstar = {k: next((ps for ps in panel if np.mean(res[k][ps]) > 0.90), None) for k in meas}
    return panel, res, pstar, undef


# All four datasets. Every measure uses the dataset's floor x0 as its baseline
# (entropy, Gini) or clipping level (ratio, k = 1), the activity threshold tau* as the
# S-score cutoff, and, for the candidate, x0 and T = tau* - x0 (Methods).
DATASETS = [
    dict(name='Davis',         file='davis_affinity.csv',       baseline=5.0,  threshold=6.0,  floor=5.0),
    dict(name='Klaeger',       file='klaeger_matrix.csv',       baseline=5.0,  threshold=6.0,  floor=5.0),
    dict(name='Anastassiadis', file='anastassiadis_matrix.csv', baseline=20.0, threshold=50.0, floor=20.0),
    dict(name='Metz',          file='metz_matrix.csv',          baseline=4.0,  threshold=6.0,  floor=4.0),
]


def load_matrix(ds):
    if ds['name'] == 'Davis':
        aff = pd.read_csv(ds['file'])
        return aff.pivot(index='Drug_Index', columns='Protein_Index', values='Affinity').fillna(5.0).values
    return pd.read_csv(ds['file'], index_col=0).values
order = ['candidate', 'entropy', 's_score', 'gini', 'ratio']
# marker/dash pairs differ per curve so near-identical p* curves (e.g. candidate,
# entropy, and gini all converge at the same panel size on Metz) stay visually
# distinguishable instead of one line hiding the others.
styles = {'candidate': ('C3', '-',  2.2, None, 1),
          'entropy':   ('C1', '--', 1.6, 'o',  1),
          's_score':   ('C0', '-.', 1.6, 's',  1),
          'gini':      ('C2', ':',  1.8, '^',  1),
          'ratio':     ('C4', '--', 1.6, 'D',  1)}

# ---------------- (checks on Klaeger, as before) ----------------
Mk = pd.read_csv('klaeger_matrix.csv', index_col=0).values
FLOOR, TAU_STAR = 5.0, 6.0
n_drugs_k = Mk.shape[0]

# orientation guard
n_active = (Mk > TAU_STAR).sum(1)
scored = n_active > 0                      # compounds D1 allows to be scored
orient_r = spearmanr(candidate(Mk, FLOOR)[scored], n_active[scored])[0]
assert orient_r < -0.3, f"candidate is inverted: corr(selectivity, n_active) = {orient_r:+.3f}"
print(f"Orientation check (Klaeger): corr(candidate selectivity, n_active) = {orient_r:+.3f} (must be negative)")

# Boundary sweep (the gate moves; T stays at tau* - x0)
betas = np.arange(5.0, 6.01, 0.25)   # B = [x0, tau*] = [5.0, 6.0]
def worst_pair(fn):
    rr = np.array([to_ranks(fn(b)) for b in betas])
    return min(spearmanr(rr[i], rr[j])[0] for i in range(len(betas)) for j in range(i + 1, len(betas)))
print("\nBoundary sweep (worst-case rank Spearman over B = [5.0,6.0]; pair-reversal fractions: d3_pair_reversals.py):")
print(f"  candidate, vary the gate        : {worst_pair(lambda b: candidate(Mk, gate=b)):+.3f}")
print(f"  existing entropy, vary baseline : {worst_pair(lambda b: entropy(Mk, b=b)):+.3f}")

# D3: monotonicity under weak off-target addition
base = candidate(Mk, FLOOR)
new = candidate(np.hstack([Mk, np.full((n_drugs_k, 1), TAU_STAR)]), FLOOR)   # worst case: at tau*
print(f"\nD3 (add a sub-threshold off-target): selectivity non-increasing for "
      f"{(new <= base + 1e-9).mean() * 100:.0f}% of compounds; max increase = {max(0.0, (new - base).max()):.4f}")

# ---------------- panel-size convergence (all four datasets) ----------------
fig, axes = plt.subplots(2, 2, figsize=(14, 10), squeeze=False)
for ax, ds in zip(axes.ravel(), DATASETS):
    M = load_matrix(ds)
    n_drugs, n_kin = M.shape
    panel, res, pstar = panel_convergence(M, ds['baseline'], ds['threshold'], ds['floor'])
    # D1 actually applied: Spearman over compounds defined on both full and subpanel.
    _, _, pstar_g, undef = gated_panel_convergence(
        M, ds['baseline'], ds['threshold'], ds['floor'], ds['threshold'])
    print(f"\nPanel-size convergence p* ({ds['name']}, {n_drugs} cpd x {n_kin} kin):")
    print(f"  {'measure':10s} {'p* (all cpd)':>13s} {'p* (D1-gated)':>14s}")
    for k in order:
        print(f"  {k:10s} {str(pstar[k]):>13s} {str(pstar_g[k]):>14s}")
    print("  fraction of compounds UNDEFINED under the D1 gate, by panel size:")
    print("   " + "  ".join(f"{ps}:{np.mean(undef[ps]) * 100:.0f}%" for ps in panel))
    print("  (gated and ungated p* are computed on different populations: gating removes")
    print("   the least-bound compounds, which are also the hardest to rank, so part of any")
    print("   improvement is selection rather than a property of the gate)")

    # Candidate vs entropy: how much do the gate and the hinge change the ranking?
    # Over the compounds both score: the candidate (D1) scores those with a kinase
    # above tau*, all of which the entropy at beta = x0 also scores.
    scored = (M > ds['threshold']).sum(1) > 0
    fns = measures_for(ds['baseline'], ds['threshold'], ds['floor'])
    r_ce = spearmanr(fns['candidate'](M)[scored], fns['entropy'](M)[scored])[0]
    print(f"  Spearman(candidate, entropy) on the full panel, {scored.sum()} scored compounds = {r_ce:+.4f}")
    for k in order:
        c, ls, lw, marker, markevery = styles[k]
        ax.plot(panel, [np.mean(res[k][ps]) for ps in panel], color=c, ls=ls, lw=lw,
                marker=marker, markevery=markevery, markersize=5,
                markerfacecolor='none' if marker else None,
                label=f"{k} (p*={pstar[k]})")
    ax.axhline(0.90, color='black', ls=':', lw=1)
    ax.set_xlabel('Panel size (kinases)'); ax.set_ylabel('Spearman rho vs full-panel ranking')
    ax.set_title(f"{ds['name']} ({n_drugs} cpd $\\times$ {n_kin} kin)")
    ax.legend(fontsize=8); ax.set_ylim(-0.2, 1.01)

fig.suptitle('Panel-size convergence: candidate (gated exponential entropy) '
             'vs existing measures (50 subsamples per size)', fontsize=12)
plt.tight_layout(rect=[0, 0, 1, 0.96])
plt.savefig('candidate_panel_convergence.png', dpi=150, bbox_inches='tight')
print("\nSaved candidate_panel_convergence.png")
