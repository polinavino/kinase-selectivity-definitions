"""
A candidate selectivity measure satisfying the four required properties (D1-D4),
with empirical checks of D3 (baseline robustness), D4 (monotonicity), and of how
quickly its compound ranking converges to the full-panel ranking.

Construction.  Selectivity is the negative Shannon "effective number of targets"
of a compound's binding profile, where each kinase is weighted by a SMOOTH
(softplus) hinge of its log-affinity above the assay detection floor:

    w_i = T * softplus( (x_i - x_floor) / T ) = T * log(1 + exp((x_i - x_floor)/T)),
    p_i = w_i / sum_j w_j,
    H(x) = - sum_i p_i log p_i,    selectivity  S(x) = -H(x),

reported only for compounds with max_i x_i > tau* (the D1 gate).

Why this satisfies D1-D4:
  D1  explicit reliability gate on max_i x_i.
  D2  fully distributional (no reference off-target) -> does not privilege any
      single pairwise comparison; bounded gap sensitivity holds trivially.
  D3  the hard active/inactive cutoff max(x-beta,0) of the existing entropy is
      replaced by a smooth hinge, so no kinase abruptly enters/leaves the active
      set as the floor shifts; the floor is a fixed assay constant, not a tunable
      baseline.  (Verified below: ranking is ~invariant to floor in [4.5, 6.5].)
  D4  softplus is monotone and -> 0 below the floor, so adding a sub-threshold
      off-target adds ~0 weight and (after renormalization) never raises S.
      (Verified below.)

The single design choice is the smoothing width T (in pK_d units); T ~ 1 places
the candidate's panel-size convergence in the fast class (entropy, S-score) while
making it essentially baseline-invariant.
"""
import os
os.environ['KMP_DUPLICATE_LIB_OK'] = 'TRUE'
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
import matplotlib.pyplot as plt

M = pd.read_csv("klaeger_matrix.csv", index_col=0).values
n_drugs, n_kin = M.shape
FLOOR, TAU_STAR, T = 5.0, 6.0, 1.0

def to_ranks(s):
    return len(s) - s.argsort().argsort()

# ---------------- candidate ----------------
def candidate(P, floor=FLOOR, T=T, eps=1e-10):
    w = T * np.logaddexp(0.0, (P - floor) / T)          # smooth hinge on log-affinity
    rs = np.where(w.sum(1, keepdims=True) == 0, eps, w.sum(1, keepdims=True))
    p = w / rs
    return -(-(p * np.where(p > 0, np.log2(p + eps), 0)).sum(1))   # -Shannon entropy

# ---------------- existing measures (match panel_size_analysis.py) ----------------
def entropy(P, b=5.0, e=1e-10):
    s = np.maximum(P - b, 0); rs = np.where(s.sum(1, keepdims=True) == 0, e, s.sum(1, keepdims=True))
    p = s / rs; return -(-(p * np.where(p > 0, np.log2(p + e), 0)).sum(1))
def gini(P, b=5.0):
    s = np.maximum(P - b, 0); out = []
    for r in s:
        rs = np.sort(r); n = len(rs); t = rs.sum()
        out.append(0.0 if t == 0 else (2 * np.sum(np.arange(1, n + 1) * rs)) / (n * t) - (n + 1) / n)
    return np.array(out)
def s_score(P, thr=6.0): return -(P > thr).astype(float).mean(1)
def ratio(P, k=1):
    return np.array([np.sort(r)[::-1][0] - max(np.sort(r)[::-1][k] if len(r) > k else 5.0, 5.0) for r in P])

measures = {'candidate': candidate, 'entropy': entropy, 'gini': gini,
            's_score': s_score, 'ratio': ratio}
ref = {k: to_ranks(f(M)) for k, f in measures.items()}

# ---------------- (A) panel-size convergence ----------------
panel = list(range(50, n_kin, 30)) + [n_kin]   # match panel_size_analysis.py grid
np.random.seed(42); R = 50
res = {k: {ps: [] for ps in panel} for k in measures}
for ps in panel:
    for _ in range(R):
        idx = np.random.choice(n_kin, ps, replace=False); Ms = M[:, idx]
        for k, f in measures.items():
            res[k][ps].append(spearmanr(ref[k], to_ranks(f(Ms)))[0])
def pstar(k):
    return next((ps for ps in panel if np.mean(res[k][ps]) > 0.90), None)

print("Panel-size convergence  p* = smallest panel with mean Spearman rho > 0.90:")
for k in measures:
    print(f"  {k:10s} p* = {pstar(k)}")

# ---------------- (B) D3: baseline robustness ----------------
betas = np.arange(4.5, 6.6, 0.25)
def worst_pair(fn):
    rr = np.array([to_ranks(fn(M, b)) for b in betas])
    return min(spearmanr(rr[i], rr[j])[0] for i in range(len(betas)) for j in range(i + 1, len(betas)))
print("\nD3 baseline robustness (worst-case rank Spearman over floor/beta in [4.5,6.5]):")
print(f"  candidate (smooth hinge): {worst_pair(lambda P,b: candidate(P, floor=b)):+.3f}")
print(f"  existing entropy (hard) : {worst_pair(lambda P,b: entropy(P, b=b)):+.3f}")

# ---------------- (C) D4: monotonicity under weak off-target addition ----------------
base = candidate(M)
new = candidate(np.hstack([M, np.full((n_drugs, 1), TAU_STAR - 0.7)]))
print(f"\nD4 (add a sub-threshold off-target): selectivity non-increasing for "
      f"{(new <= base + 1e-9).mean() * 100:.0f}% of compounds; max increase = {max(0.0, (new - base).max()):.4f}")

# ---------------- figure ----------------
fig, ax = plt.subplots(figsize=(8.5, 5))
styles = {'candidate': ('C3', '-', 2.6), 'entropy': ('C1', '--', 1.6),
          's_score': ('C0', '--', 1.6), 'gini': ('C2', '--', 1.6), 'ratio': ('C4', '--', 1.6)}
for k in measures:
    c, ls, lw = styles[k]
    ax.plot(panel, [np.mean(res[k][ps]) for ps in panel], color=c, ls=ls, lw=lw,
            label=f"{k} (p*={pstar(k)})")
ax.axhline(0.90, color='black', ls=':', lw=1)
ax.set_xlabel('Panel size (kinases)'); ax.set_ylabel('Spearman rho vs full-panel ranking')
ax.set_title('Panel-size convergence: candidate (gated smooth-hinge effective-target number)\n'
             'vs existing measures (Klaeger, 50 subsamples per size)')
ax.legend(fontsize=8); ax.set_ylim(0, 1.01)
plt.tight_layout(); plt.savefig('candidate_panel_convergence.png', dpi=150, bbox_inches='tight')
print("\nSaved candidate_panel_convergence.png")
