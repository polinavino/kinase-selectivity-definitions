"""
Illustrates the two failure modes the candidate measure is built to repair
(D1: arbitrary scores for non-binders, D3: baseline-driven instability), using
the same Klaeger computations that produced the numbers quoted in the paper
(paper/sections/desiderata.tex, D1 and D3 motivation paragraphs).

Left panel (D1): pooled rank instability (std across all four existing
definitions' parameter sweeps) for compounds with no detected binding above the
assay floor vs. compounds with at least one active kinase. The candidate does
not assign these zero-active compounds a score at all (D1 gate), rather than a
low-instability one -- shown as a hatch pattern, not a bar height.

Right panel (D3): per-compound rank instability under the free activity
baseline (entropy: baseline in [5.0, 6.5]; candidate: emphasis baseline beta in
[4.5, 6.5] with the D1 gate pinned at the assay floor) vs. n_active, for active
compounds. Entropy's instability falls off with n_active; the candidate's stays
flat near zero because it has no free activity baseline.
"""
import script_logging; script_logging.capture(__file__)
import os
os.environ['KMP_DUPLICATE_LIB_OK'] = 'TRUE'
import numpy as np
import pandas as pd
from scipy.stats import spearmanr, rankdata
import matplotlib.pyplot as plt

FLOOR, TAU_STAR = 5.0, 6.0

def to_ranks(s):
    return rankdata(-np.asarray(s), method='average')

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
def candidate(P, baseline, floor, T=1.0, eps=1e-10):
    gate = (P > floor).astype(float)
    w = gate * (T * np.logaddexp(0.0, (P - baseline) / T))
    rs = np.where(w.sum(1, keepdims=True) == 0, eps, w.sum(1, keepdims=True))
    p = w / rs
    return -(-(p * np.where(p > 0, np.log2(p + eps), 0)).sum(1))

M = pd.read_csv('klaeger_matrix.csv', index_col=0).values
n_drugs, n_kin = M.shape
n_active = (M > TAU_STAR).sum(1)
has_active = n_active > 0

# ---------------- D1: pooled rank instability, zero-active vs active ----------------
s_thresholds, ent_baselines, gini_baselines, ratio_top_ns = (
    np.arange(5.5, 8.25, 0.25), np.arange(5.0, 6.75, 0.25),
    np.arange(5.0, 6.75, 0.25), list(range(1, 6)))
all_ranks = np.vstack([
    [to_ranks(s_score(M, t)) for t in s_thresholds],
    [to_ranks(entropy(M, b)) for b in ent_baselines],
    [to_ranks(gini(M, b)) for b in gini_baselines],
    [to_ranks(ratio(M, FLOOR, n)) for n in ratio_top_ns],
])
rank_std = all_ranks.std(axis=0)
sigma_zero, sigma_active = rank_std[~has_active].mean(), rank_std[has_active].mean()
print(f"D1: rank std, zero-active = {sigma_zero:.1f}; active = {sigma_active:.1f} "
      f"(n_zero={ (~has_active).sum() }, n_active={has_active.sum()})")

# ---------------- D3: per-compound instability under the free baseline ----------------
betas = np.arange(4.5, 6.6, 0.25)
ent_ranks_b = np.array([to_ranks(entropy(M, b)) for b in np.arange(5.0, 6.6, 0.25)])
ent_instab = ent_ranks_b.std(axis=0)
cand_ranks_b = np.array([to_ranks(candidate(M, baseline=b, floor=FLOOR)) for b in betas])
cand_instab = cand_ranks_b.std(axis=0)

r_ent, p_ent = spearmanr(n_active[has_active], ent_instab[has_active])
r_cand, p_cand = spearmanr(n_active[has_active], cand_instab[has_active])
print(f"D3: corr(n_active, instability) -- entropy r={r_ent:+.3f} (p={p_ent:.1e}); "
      f"candidate r={r_cand:+.3f} (p={p_cand:.1e})")
print(f"D3: mean instability -- entropy {ent_instab[has_active].mean():.2f}; "
      f"candidate {cand_instab[has_active].mean():.2f} "
      f"({ent_instab[has_active].mean() / cand_instab[has_active].mean():.0f}x smaller)")

# ---------------- figure ----------------
fig, (axL, axR) = plt.subplots(1, 2, figsize=(11, 4.5))

bars = axL.bar(['No detected\nbinding', 'Has an active\nkinase'],
                [sigma_zero, sigma_active], color=['C1', 'C0'], width=0.55)
bars[0].set_hatch('///'); bars[0].set_edgecolor('white')
for b, v in zip(bars, [sigma_zero, sigma_active]):
    axL.text(b.get_x() + b.get_width() / 2, v + 1.5, f"{v:.0f}", ha='center', fontsize=11)
axL.set_ylabel('Rank std across existing definitions\' parameter sweeps')
axL.set_title('D1: rank instability by activity status (Klaeger)')
axL.text(0, sigma_zero + 9, 'candidate: undefined\n(D1 gate), not scored',
          ha='center', fontsize=8.5, color='C1')
axL.set_ylim(0, max(sigma_zero, sigma_active) * 1.35)

axR.scatter(n_active[has_active], ent_instab[has_active], s=16, color='C1', alpha=0.6,
            label=f'entropy (r={r_ent:+.2f})')
axR.scatter(n_active[has_active], cand_instab[has_active], s=16, color='C3', alpha=0.6,
            label=f'candidate (r={r_cand:+.2f})')
axR.set_xlabel('n_active (kinases bound)'); axR.set_ylabel('Rank std across free baseline')
axR.set_title('D3: baseline instability vs. active-set size (Klaeger)')
axR.legend(fontsize=9); axR.set_ylim(bottom=-2)
ratio = ent_instab[has_active].mean() / cand_instab[has_active].mean()
axR.text(0.97, 0.95, f"mean instability\n{ratio:.0f}$\\times$ smaller\nfor candidate",
          ha='right', va='top', fontsize=8.5, transform=axR.transAxes)

fig.suptitle('The two failure modes the candidate measure repairs', fontsize=12)
plt.tight_layout(rect=[0, 0, 1, 0.94])
plt.savefig('failure_modes.png', dpi=150, bbox_inches='tight')
print("\nSaved failure_modes.png")
