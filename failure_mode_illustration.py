"""
Illustrates two failure modes (D1: scores for compounds with no binding above the
activity threshold; D2: baseline-driven instability), using
the same Klaeger computations that produced the numbers quoted in the paper
(paper/sections/desiderata.tex, D1 and D2 paragraphs).

Left panel (D1): pooled rank instability (std across all four existing
definitions' parameter sweeps) for compounds with no binding above the pK_d > 6
activity threshold vs. compounds with at least one active kinase. 16 of the 22
such compounds on the curated panel bind some kinase above the floor of
pK_d = 5.0; "zero-active" is defined against the activity threshold, not the
floor. The candidate does not assign these
compounds a score at all (D1 gate), rather than a low-instability one -- shown as
a hatch pattern, not a bar height.

The pooled statistic is also decomposed per definition below, because it is not
uniform across the four: it is driven by the Gini, is small for the entropy and
the ratio, and reverses in sign for the S-score (whose score is exactly 0 for every
zero-active compound at thresholds of pK_d 6 or above, so its ranks are stable
by degeneracy rather than by reliability). The paper reports the decomposition rather than claiming the effect
holds "under all four definitions".

Right panel (D2): per-compound rank instability vs. n_active, for active
compounds, under three different sweeps (Reviewer 1, Molecular Informatics
revision: the original version of this panel compared entropy's sensitivity to
its gate/baseline against only the candidate's emphasis parameter, which does
not touch the gate -- an apples-to-oranges comparison the reviewer correctly
flagged). All three series are swept over B = [x0, tau*] = [5.0, 6.0]:
  (1) entropy, baseline in B -- this baseline is BOTH the entropy's
      gate (active/inactive cutoff) AND its emphasis weighting, left fully
      free in existing practice;
  (2) candidate, floor x0 in B (with emphasis beta pinned to x0 at
      each point) -- the TRUE equivalent comparison, sweeping the candidate's
      own gate exactly as entropy's gate is swept in (1). This is comparably,
      if not more, sensitive than entropy, because gate placement always
      moves kinases in and out of the active set;
  (3) candidate, emphasis beta in B with the gate FIXED at the assay
      floor (x0 = 5.0, as used throughout the rest of this paper) -- the
      candidate's only remaining free parameter once x0 is fixed. Fixing the
      gate is a protocol decision: swept as in (2), the candidate fails D2 as
      the entropy does, and only its remaining emphasis parameter is robust.
"""
import script_logging; script_logging.capture(__file__)
import os
os.environ['KMP_DUPLICATE_LIB_OK'] = 'TRUE'
import numpy as np
import pandas as pd
from scipy.stats import spearmanr, rankdata, mannwhitneyu
import matplotlib.pyplot as plt

FLOOR, TAU_STAR = 5.0, 6.0

def to_ranks(s):
    return rankdata(-np.asarray(s), method='average')

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
def candidate(P, baseline, floor, T=1.0, eps=1e-10):
    gate = (P > floor).astype(float)
    w = gate * (T * np.logaddexp(0.0, (P - baseline) / T))
    rs = np.where(w.sum(1, keepdims=True) == 0, eps, w.sum(1, keepdims=True))
    p = w / rs
    return -(-(p * np.where(p > 0, np.log2(np.where(p > 0, p, 1.0)), 0)).sum(1))

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
u_stat, u_p = mannwhitneyu(rank_std[~has_active], rank_std[has_active], alternative='two-sided')
print(f"D1: rank std, zero-active = {sigma_zero:.1f}; active = {sigma_active:.1f} "
      f"(n_zero={ (~has_active).sum() }, n_active={has_active.sum()}; "
      f"Mann-Whitney U={u_stat:.0f}, p={u_p:.2e}, two-sided)")
print(f"D1: max pK_d among the {(~has_active).sum()} zero-active compounds: "
      f"{M[~has_active].max(1).min():.2f} to {M[~has_active].max(1).max():.2f}; "
      f"{int((M[~has_active].max(1) > FLOOR).sum())} bind some kinase above the floor of {FLOOR}")

# Per-definition decomposition: the pooled figure above is NOT uniform across the four.
print("\nD1 decomposed per definition (mean rank std within that definition's own sweep):")
per_def = {
    's_score': [to_ranks(s_score(M, t)) for t in s_thresholds],
    'entropy': [to_ranks(entropy(M, b)) for b in ent_baselines],
    'gini':    [to_ranks(gini(M, b)) for b in gini_baselines],
    'ratio':   [to_ranks(ratio(M, FLOOR, n)) for n in ratio_top_ns],
}
for name, ranks in per_def.items():
    sd = np.array(ranks).std(axis=0)
    z, a = sd[~has_active].mean(), sd[has_active].mean()
    verdict = 'holds' if z > a else 'REVERSES'
    p_def = mannwhitneyu(sd[~has_active], sd[has_active]).pvalue
    print(f"  {name:8s} zero-active {z:6.2f}   active {a:6.2f}   {verdict}   (Mann-Whitney p={p_def:.2g})")
print("  (the S-score reverses because its score is exactly 0 for every zero-active")
print("   compound at thresholds of pK_d 6 or above, i.e. stable by degeneracy)")

# ---------------- D2: per-compound instability under three different sweeps ----------------
# one grid for all three sweeps: B = [x0, tau*] = [5.0, 6.0] (Methods, D2), so the
# candidate's gate is swept over exactly the entropy's range
sweep = np.arange(5.0, 6.01, 0.25)
ent_ranks_b = np.array([to_ranks(entropy(M, b)) for b in sweep])
ent_instab = ent_ranks_b.std(axis=0)
# true equivalent: sweep the candidate's own gate (floor x0), beta pinned to x0 at each point
cand_floor_ranks = np.array([to_ranks(candidate(M, baseline=f, floor=f)) for f in sweep])
cand_floor_instab = cand_floor_ranks.std(axis=0)
# the candidate's remaining free parameter once the gate is fixed by declaration (x0 = FLOOR)
cand_emph_ranks = np.array([to_ranks(candidate(M, baseline=b, floor=FLOOR)) for b in sweep])
cand_emph_instab = cand_emph_ranks.std(axis=0)

r_ent, p_ent = spearmanr(n_active[has_active], ent_instab[has_active])
r_cand_floor, p_cf = spearmanr(n_active[has_active], cand_floor_instab[has_active])
r_cand_emph, p_ce = spearmanr(n_active[has_active], cand_emph_instab[has_active])
print(f"D2: corr(n_active, instability) -- entropy r={r_ent:+.3f} (p={p_ent:.1e}); "
      f"candidate-floor-swept r={r_cand_floor:+.3f} (p={p_cf:.1e}); "
      f"candidate-emphasis-swept (gate fixed) r={r_cand_emph:+.3f} (p={p_ce:.1e})")
ratio_gate = cand_floor_instab[has_active].mean() / ent_instab[has_active].mean()
ratio_emph = ent_instab[has_active].mean() / cand_emph_instab[has_active].mean()
print(f"D2: mean instability over B=[5.0,6.0] -- entropy (baseline swept) {ent_instab[has_active].mean():.2f}; "
      f"candidate (gate/floor swept, true equivalent) {cand_floor_instab[has_active].mean():.2f} "
      f"({ratio_gate:.2f}x the entropy's); "
      f"candidate (emphasis swept, gate fixed by declaration) {cand_emph_instab[has_active].mean():.2f} "
      f"({ratio_emph:.0f}x smaller than entropy -- but only because the gate itself is no longer "
      f"a free parameter being compared)")

# ---------------- figure ----------------
fig, (axL, axR) = plt.subplots(1, 2, figsize=(11, 4.5))

bars = axL.bar(['No active\nkinase', 'Has an active\nkinase'],
                [sigma_zero, sigma_active], color=['C1', 'C0'], width=0.55)
bars[0].set_hatch('///'); bars[0].set_edgecolor('white')
for b, v in zip(bars, [sigma_zero, sigma_active]):
    axL.text(b.get_x() + b.get_width() / 2, v + 1.5, f"{v:.0f}", ha='center', fontsize=11)
axL.set_ylabel('Rank std across existing definitions\' parameter sweeps')
axL.set_title('D1: rank instability by activity status (Klaeger)')
axL.text(0, sigma_zero + 9, 'not scored by\nthe candidate',
          ha='center', fontsize=8.5, color='C1')
axL.set_ylim(0, max(sigma_zero, sigma_active) * 1.35)

axR.scatter(n_active[has_active], ent_instab[has_active], s=16, color='C1', alpha=0.6,
            label=f'entropy, baseline swept (r={r_ent:+.2f})')
axR.scatter(n_active[has_active], cand_floor_instab[has_active], s=16, color='C2', alpha=0.6,
            label=f'candidate, gate swept (r={r_cand_floor:+.2f})')
axR.scatter(n_active[has_active], cand_emph_instab[has_active], s=16, color='C3', alpha=0.6,
            label=f'candidate, gate fixed, emphasis swept (r={r_cand_emph:+.2f})')
axR.set_xlabel('n_active (kinases bound)'); axR.set_ylabel('Rank std across the swept parameter')
axR.set_title('D2: rank instability under each parameter sweep (Klaeger)')
axR.legend(fontsize=7.5, loc='upper right'); axR.set_ylim(bottom=-2)
axR.text(0.97, 0.55, f"mean: entropy {ent_instab[has_active].mean():.1f}, "
          f"gate {cand_floor_instab[has_active].mean():.1f}, "
          f"emphasis {cand_emph_instab[has_active].mean():.1f}",
          ha='right', va='top', fontsize=7.5, transform=axR.transAxes)

plt.tight_layout()
plt.savefig('failure_modes.png', dpi=150, bbox_inches='tight')
print("\nSaved failure_modes.png")
