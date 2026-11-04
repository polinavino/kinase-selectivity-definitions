"""
Checks of specific claims in the manuscript that the four datasets can test.

Each section prints the numbers cited in the text for one claim:
  1. Type 1 compounds: were their known targets quantified in Klaeger at all?
  2. Klaeger within-family correlations with and without the zero-active compounds,
     with bootstrap intervals for whether the Gini sits closer to the S-score than
     to the ratio.
  3. Table 1 with mean instead of median ranks.
  4. Anastassiadis: the top1-top2 gap vs ratio instability, without the compounds
     whose two best-inhibited kinases are saturated at the 100% ceiling.
  5. Klaeger entropy-ratio disagreement among active compounds: the ten largest.
  6. Which compounds are most unstable (Klaeger, top 20), by instability type.
  7. Mean rank instability by number of active kinases (Klaeger, Davis).
  8. Entropy instability vs the share of kinases that cross the baseline sweep.
  9. The 46 compounds shared by Davis and Klaeger: active fractions, and whether
     within-dataset rank instability predicts cross-platform rank disagreement.
 10. Metz: dependence of the rankings on how many kinases each compound was
     measured against (floor imputation of untested pairs).
 11. Within-family correlation range on each dataset, and the ratio's range against
     the family, for all compounds and for those with an active kinase
     (Introduction, Results).
 12. Compounds shared by Klaeger and Anastassiadis, by name (Discussion).
 13. Zero-active vs active rank instability per definition, on every dataset with
     zero-active compounds (Results, replication).
 14. Gini panel-size p* on every dataset when only compounds whose Gini is undefined
     on the subpanel or the full panel are excluded (Results, gated panel size).
 15. Klaeger: does the active affinity range predict instability once the number
     of active kinases is held fixed? (Results, instability predictors)
 16. Klaeger: does the weight near the threshold explain the skewness-Gini
     instability association? (Results, instability predictors)
 17. The 46 compounds shared by Davis and Klaeger: does active affinity range
     predict S-score instability on both platforms? (Results, instability predictors)
 18. The most stable compounds on each dataset: where each definition ranks them,
     and how much of the broad Davis compounds' instability comes from the ratio
     (Results, Type 3 and stability patterns)
 19. D4 under realistic panel expansion: does adding real kinases that a compound
     binds below threshold to a focused panel raise the entropy's or the
     candidate's selectivity score? (Results, "Why the Smooth Hinge Matters")
"""
import script_logging; script_logging.capture(__file__)
import os
import numpy as np
import pandas as pd
from scipy.stats import spearmanr, rankdata

HERE = os.path.dirname(os.path.abspath(__file__))
FAM = ['s_score', 'entropy', 'gini', 'ratio']


def to_ranks(s):
    return rankdata(-np.asarray(s), method='average')            # rank 1 = most selective


def ent(M, b, e=1e-10):
    s = np.maximum(M - b, 0)
    rs = np.where(s.sum(1, keepdims=True) == 0, e, s.sum(1, keepdims=True))
    p = s / rs
    return (p * np.where(p > 0, np.log2(np.where(p > 0, p, 1.0)), 0)).sum(1)        # selectivity orientation (-H)


def gini(M, b):
    s = np.sort(np.maximum(M - b, 0), 1)
    n = s.shape[1]
    t = s.sum(1)
    g = (2 * (np.arange(1, n + 1) * s).sum(1)) / (n * np.where(t == 0, 1, t)) - (n + 1) / n
    return np.where(t == 0, 0.0, g)


def s_score(M, t):
    return -(M > t).mean(1)


def ratio(M, floor, k=1):
    s = -np.sort(-M, 1)
    return s[:, 0] - np.maximum(s[:, k], floor)


def sweeps(M, s_thr, betas, floor):
    """Rank vectors for every configuration of each definition (as in the main analysis)."""
    return {'s_score': np.array([to_ranks(s_score(M, t)) for t in s_thr]),
            'entropy': np.array([to_ranks(ent(M, b)) for b in betas]),
            'gini':    np.array([to_ranks(gini(M, b)) for b in betas]),
            'ratio':   np.array([to_ranks(ratio(M, floor, k)) for k in range(1, 6)])}


def corr_cells(R, agg=np.median, idx=None):
    med = {f: agg(R[f], 0) for f in FAM}
    if idx is not None:
        med = {f: med[f][idx] for f in FAM}
    return {(a, b): spearmanr(med[a], med[b])[0] for i, a in enumerate(FAM) for b in FAM[i + 1:]}


aff = pd.read_csv(f"{HERE}/davis_affinity.csv")
Dv = aff.pivot(index="Drug_Index", columns="Protein_Index", values="Affinity").fillna(5.0)
Kl = pd.read_csv(f"{HERE}/klaeger_matrix.csv", index_col=0)
An = pd.read_csv(f"{HERE}/anastassiadis_matrix.csv", index_col=0)
Mz_df = pd.read_csv(f"{HERE}/metz_matrix.csv", index_col=0)
Md, Mk, Ma, Mz = Dv.values, Kl.values, An.values, Mz_df.values
SETTINGS = {
    'Davis':         (Md, np.arange(5.5, 8.25, .25), np.arange(5.0, 6.75, .25), 5.0, 6.0),
    'Klaeger':       (Mk, np.arange(5.5, 8.25, .25), np.arange(5.0, 6.75, .25), 5.0, 6.0),
    'Anastassiadis': (Ma, np.arange(50, 95, 5),      np.arange(0, 35, 5),       0.0, 50.0),
    'Metz':          (Mz, np.arange(5.0, 7.75, .25), np.arange(4.0, 5.75, .25), 4.0, 6.0),
}
R = {k: sweeps(*v[:4]) for k, v in SETTINGS.items()}
na_k, na_d = (Mk > 6).sum(1), (Md > 6).sum(1)
sk = -np.sort(-Mk, 1)
gap_k = sk[:, 0] - sk[:, 1]

print("1. Type 1 compounds: were their known targets quantified for any Klaeger compound?")
kb = pd.read_excel(f"{HERE}/aan4368_Table_S2.xlsx", sheet_name='Kinobeads')
kb['pKd'] = 9 - np.log10(kb['Apparent Kd'])
known = {'AZD-8055': 'MTOR', 'BMS-911543': 'JAK2', 'Vatalanib': 'KDR', 'Filgotinib': 'JAK1', 'ARRY-380': 'ERBB2'}
print(f"   ({kb['Drug'].nunique()} compounds in the Kinobeads table)")
for drug, g in known.items():
    rows = kb[kb['Gene Name'] == g]
    own = rows[rows['Drug'] == drug]
    print(f"   {drug:11s} -> {g:5s}: fitted for {rows['Drug'].nunique():2d} compounds "
          f"({(rows['Target Classification'] == 'High confidence').sum()} high-confidence rows, "
          f"max pKd {rows['pKd'].max() if len(rows) else float('nan'):.2f}); for {drug}: {'yes' if len(own) else 'no'}")

print("\n2. Klaeger within-family correlations, all compounds vs compounds with an active kinase")
act_k = na_k > 0
# active-only: re-rank among the active compounds at every configuration
R_act = sweeps(Mk[act_k], *SETTINGS['Klaeger'][1:4])
for label, Rx in (('all 222', R['Klaeger']), (f'active only ({act_k.sum()})', R_act)):
    c = corr_cells(Rx)
    print(f"   {label:18s}: S-E {c[('s_score','entropy')]:.3f}  S-G {c[('s_score','gini')]:.3f}  "
          f"E-G {c[('entropy','gini')]:.3f}  | S-R {c[('s_score','ratio')]:.3f}  E-R {c[('entropy','ratio')]:.3f}  "
          f"G-R {c[('gini','ratio')]:.3f}")
rng = np.random.default_rng(0)
for label, Rx in (('all 222', R['Klaeger']), (f'active only ({act_k.sum()})', R_act)):
    med = {f: np.median(Rx[f], 0) for f in FAM}
    pool = np.arange(len(med['gini']))
    diff = []
    for _ in range(1000):
        ix = rng.choice(pool, len(pool))
        diff.append(spearmanr(med['s_score'][ix], med['gini'][ix])[0] - spearmanr(med['gini'][ix], med['ratio'][ix])[0])
    print(f"   {label:18s}: r(S-score, Gini) - r(Gini, ratio), bootstrap 95% CI "
          f"[{np.percentile(diff, 2.5):+.6f}, {np.percentile(diff, 97.5):+.6f}]")

print("\n3. Table 1 cells with mean instead of median ranks: largest change per dataset")
for ds in SETTINGS:
    a, b = corr_cells(R[ds]), corr_cells(R[ds], agg=np.mean)
    print(f"   {ds:13s}: max |median-based - mean-based| = {max(abs(a[p] - b[p]) for p in a):.3f}")

print("\n4. Anastassiadis: top1-top2 gap vs ratio instability")
sa = -np.sort(-Ma, 1)
gap_a = sa[:, 0] - sa[:, 1]
act_a = (Ma > 50).any(1)
ceil = sa[:, 1] >= 99.5                       # both best-inhibited kinases at the ceiling
inst_r = R['Anastassiadis']['ratio'].std(0)
r1, p1 = spearmanr(gap_a[act_a], inst_r[act_a])
r2, p2 = spearmanr(gap_a[act_a & ~ceil], inst_r[act_a & ~ceil])
print(f"   active compounds {act_a.sum()}; top two both >= 99.5% inhibition: {(act_a & ceil).sum()}")
print(f"   all active: r = {r1:+.3f} (p = {p1:.2g});  excluding ceiling-saturated: r = {r2:+.3f} (p = {p2:.2g}, "
      f"n = {(act_a & ~ceil).sum()})")

print("\n5. Klaeger: ten largest entropy-ratio disagreements among active compounds")
dis = np.abs(np.median(R['Klaeger']['entropy'], 0) - np.median(R['Klaeger']['ratio'], 0))
top = [i for i in np.argsort(-dis) if act_k[i]][:10]
for i in top:
    print(f"   {Kl.index[i]:14s} n_active = {na_k[i]:3d}  top1-top2 gap = {gap_k[i]:.2f}")
print(f"   n_active range among them: {min(na_k[top])}-{max(na_k[top])}; "
      f"with gap >= 1: {sum(gap_k[i] >= 1 for i in top)} of 10")

print("\n6. Klaeger: the 20 most unstable compounds (overall rank sigma) by type")
sig_k = np.vstack([R['Klaeger'][f] for f in FAM]).std(0)
types = []
for i in np.argsort(-sig_k)[:20]:
    types.append('Type 1 (no active kinase)' if na_k[i] == 0 else
                 'Type 2 (top1-top2 gap < 0.1)' if gap_k[i] < 0.1 else
                 'Type 3 (>= 10 active, gap >= 1)' if na_k[i] >= 10 and gap_k[i] >= 1 else 'other')
print("   " + "; ".join(f"{k}: {v}" for k, v in pd.Series(types).value_counts().items()))

print("\n7. Mean overall rank sigma by number of active kinases")
sig_d = np.vstack([R['Davis'][f] for f in FAM]).std(0)
for ds, sig, na in (('Klaeger', sig_k, na_k), ('Davis', sig_d, na_d)):
    parts = []
    for a, b in ((1, 2), (3, 5), (6, 10), (11, 20), (21, 40), (41, 10**6)):
        m = (na >= a) & (na <= b)
        if m.any():
            parts.append(f"{a}-{b if b < 10**6 else '+'}: {sig[m].mean():.1f} (n={m.sum()})")
    print(f"   {ds:8s}: " + "; ".join(parts))

print("\n8. Klaeger: entropy instability vs the share of above-floor kinases crossing the baseline sweep")
churn = ((Mk > 5.0) & (Mk <= 6.5)).sum(1) / np.maximum((Mk > 5.0).sum(1), 1)
r, p = spearmanr(churn[act_k], R['Klaeger']['entropy'].std(0)[act_k])
print(f"   active compounds: r = {r:+.3f} (p = {p:.2g})")

print("\n9. The 46 compounds profiled in both Davis and Klaeger")
ov = pd.read_csv(f"{HERE}/davis_klaeger_overlap.csv")
di = [list(Dv.index).index(i) for i in ov['Davis_Drug_Index']]
ki = [list(Kl.index).index(n) for n in ov['Klaeger_Name']]
print(f"   mean active fraction of the panel: Davis {100 * (na_d[di] / Md.shape[1]).mean():.1f}%, "
      f"Klaeger {100 * (na_k[ki] / Mk.shape[1]).mean():.1f}%; paired Spearman of n_active {spearmanr(na_d[di], na_k[ki])[0]:+.2f}")
for f in FAM:
    d = np.abs(rankdata(np.median(R['Davis'][f], 0)[di]) - rankdata(np.median(R['Klaeger'][f], 0)[ki]))
    print(f"   {f:8s}: Spearman(|Davis rank - Klaeger rank|, within-dataset sigma): "
          f"Davis {spearmanr(d, sig_d[di])[0]:+.2f}, Klaeger {spearmanr(d, sig_k[ki])[0]:+.2f}")

print("\n10. Metz: dependence on per-compound coverage (untested pairs set to the floor)")
mz = pd.read_excel(f"{HERE}/metz.xls", sheet_name=0, header=0)
meta = ['Cmpd_ID', 'PUBCHEM_SID', 'Canonical_Smiles', 'External_Cmpd_ID', 'External_Source', 'Cluster',
        'ClusterSize', 'Cluster_MCSS', 'Molecular_Weight', 'ALogP', 'Num_H_Acceptors', 'Num_H_Donors', 'tPSA',
        'Promiscuity_1uM']
V = mz.set_index('Cmpd_ID')[[c for c in mz.columns if c not in meta]].apply(pd.to_numeric, errors='coerce')
V = V.reindex(Mz_df.index)
measured = V.notna().values
cov = measured.sum(1)
vals = V.values
print(f"   measured entries above the 4.0 floor: {100 * np.mean(vals[measured] > 4.0):.1f}%")
for f in FAM:
    print(f"   Spearman(coverage, median {f:8s} rank; rank 1 = most selective) = "
          f"{spearmanr(cov, np.median(R['Metz'][f], 0))[0]:+.3f}")
frac_meas = np.nansum(vals > 6, 1) / cov
print(f"   Spearman(coverage, active fraction among measured kinases) = {spearmanr(cov, frac_meas)[0]:+.3f}")
kfrac = measured.mean(0)
for thr in (0.8, 0.9):
    ks = kfrac >= thr
    print(f"   kinases measured for >= {int(100 * thr)}% of compounds: {ks.sum()}")

print("\n11. Within-family (S-score, entropy, Gini) correlations of median ranks, per dataset")
lo, hi = 1, -1
for ds in SETTINGS:
    c = corr_cells(R[ds])
    w = [c[('s_score', 'entropy')], c[('s_score', 'gini')], c[('entropy', 'gini')]]
    x = [c[(f, 'ratio')] for f in FAM[:3]]
    lo, hi = min(lo, min(w)), max(hi, max(w))
    print(f"   {ds:13s}: {min(w):.3f}-{max(w):.3f} | ratio vs family {min(x):.3f}-{max(x):.3f}")
print(f"   across all four datasets: {lo:.3f}-{hi:.3f}")
print("   compounds with an active kinase, re-ranked among themselves: within-family | ratio vs family")
for ds, (M, s_thr, betas, floor, tau) in SETTINGS.items():
    act = (M > tau).sum(1) > 0
    c = corr_cells(sweeps(M[act], s_thr, betas, floor))
    w = [c[('s_score', 'entropy')], c[('s_score', 'gini')], c[('entropy', 'gini')]]
    x = [c[(f, 'ratio')] for f in FAM[:3]]
    print(f"   {ds:13s} ({act.sum()} of {len(act)}): {min(w):.3f}-{max(w):.3f} | {min(x):.3f}-{max(x):.3f}")

print("\n12. Compounds shared by Klaeger and Anastassiadis (names, whitespace-trimmed)")
shared = set(n.strip() for n in Kl.index) & set(n.strip() for n in An.index)
print(f"   {len(shared)} shared: {sorted(shared)}")

print("\n13. Zero-active vs active mean rank sigma within each definition's own sweep")
from scipy.stats import mannwhitneyu
for ds, (M, *_rest) in SETTINGS.items():
    act = (M > _rest[-1]).sum(1) > 0
    if act.all():
        print(f"   {ds:13s}: no zero-active compounds")
        continue
    cells = []
    for f in FAM:
        sd = R[ds][f].std(0)
        cells.append(f"{f} {sd[~act].mean():.1f} vs {sd[act].mean():.1f} "
                     f"(p={mannwhitneyu(sd[~act], sd[act]).pvalue:.2g})")
    print(f"   {ds:13s} ({(~act).sum()} zero-active): " + "; ".join(cells))

print("\n14. Gini panel-size p* (beta = x0): all compounds, excluding only compounds whose Gini")
print("    is undefined on the subpanel or the full panel, and with the D1 gate (no kinase above tau*)")
for ds, M, x0, tau in (('Davis', Md, 5.0, 6.0), ('Klaeger', Mk, 5.0, 6.0),
                       ('Anastassiadis', Ma, 20.0, 50.0), ('Metz', Mz, 4.0, 6.0)):
    nk = M.shape[1]
    sizes = list(range(50, nk, 30)) + [nk]
    ref = gini(M, x0)
    pstar = {}
    for rule in ('all', 'undefined excluded', 'D1 gate'):
        rng = np.random.RandomState(42)
        curve = []
        for ps in sizes:
            rs = []
            for _ in range(50):
                idx = rng.choice(nk, ps, replace=False)
                if rule == 'all':
                    keep = np.ones(len(M), bool)
                elif rule == 'undefined excluded':
                    keep = (M[:, idx] > x0).any(1) & (M > x0).any(1)
                else:
                    keep = (M[:, idx] > tau).any(1) & (M > tau).any(1)
                rs.append(spearmanr(to_ranks(ref[keep]), to_ranks(gini(M[:, idx], x0)[keep]))[0])
            curve.append(np.mean(rs))
        pstar[rule] = next((p for p, r in zip(sizes, curve) if r > 0.90), None)
    print(f"   {ds:13s}: p* all {pstar['all']}, undefined excluded {pstar['undefined excluded']}, D1 gate {pstar['D1 gate']}")

print("\n15. Klaeger active compounds: active affinity range vs per-family instability,")
print("    holding the number of active kinases fixed (partial Spearman)")
from scipy.stats import pearsonr
act_rows = Mk[act_k]
rng_k = np.array([np.ptp(r[r > 6]) if (r > 6).sum() >= 2 else 0.0 for r in act_rows])  # 0 placeholder, as Table 2
na_act = na_k[act_k]


def partial_spearman(x, y, z):
    rx, ry, rz = rankdata(x), rankdata(y), rankdata(z)
    ex = rx - np.polyval(np.polyfit(rz, rx, 1), rz)
    ey = ry - np.polyval(np.polyfit(rz, ry, 1), rz)
    return pearsonr(ex, ey)


print(f"   Spearman(range, n_active) = {spearmanr(rng_k, na_act)[0]:+.3f} (n={act_k.sum()})")
for f in ['s_score', 'entropy', 'gini']:
    sig = R['Klaeger'][f].std(0)[act_k]
    r0 = spearmanr(rng_k, sig)
    r1 = partial_spearman(rng_k, sig, na_act)
    print(f"   {f:8s}: r(range, sigma) = {r0[0]:+.3f} (p = {r0[1]:.2g}); "
          f"partial on n_active = {r1[0]:+.3f} (p = {r1[1]:.2g})")


def partial_spearman_multi(x, y, covs):
    Z = np.column_stack([rankdata(c) for c in covs] + [np.ones(len(x))])
    res = lambda v: rankdata(v) - Z @ np.linalg.lstsq(Z, rankdata(v), rcond=None)[0]
    return pearsonr(res(x), res(y))


print("\n16. Klaeger compounds with >= 3 active kinases: skewness, the weight the beta sweep")
print("    crosses, and Gini instability")
from scipy.stats import skew
three = na_k >= 3
X3 = Mk[three]
na3 = na_k[three]
sk3 = np.array([skew(r[r > 6]) for r in X3])
w3 = np.maximum(X3 - 5.0, 0)
wcross = (w3 * (X3 <= 6.5)).sum(1) / w3.sum(1)       # share of above-floor weight on kinases the sweep [5.0, 6.5] crosses
sigG3 = R['Klaeger']['gini'].std(0)[three]
r = spearmanr(sk3, wcross)
print(f"   n = {three.sum()}; Spearman(skewness, weight share crossed) = {r[0]:+.3f} (p = {r[1]:.2g})")
r = partial_spearman_multi(wcross, sigG3, [na3])
print(f"   partial(weight share crossed, Gini sigma | n_active) = {r[0]:+.3f} (p = {r[1]:.2g})")
r = partial_spearman_multi(sk3, sigG3, [na3])
r2 = partial_spearman_multi(sk3, sigG3, [na3, wcross])
print(f"   partial(skewness, Gini sigma | n_active) = {r[0]:+.3f} (p = {r[1]:.2g}); "
      f"also holding the weight share fixed = {r2[0]:+.3f} (p = {r2[1]:.2g})")

print("\n17. The 46 shared Davis/Klaeger compounds: active affinity range vs S-score instability")
ov46 = pd.read_csv(f"{HERE}/davis_klaeger_overlap.csv")
for ds, frame, key in (('Davis', Dv, 'Davis_Drug_Index'), ('Klaeger', Kl, 'Klaeger_Name')):
    idx = frame.index.get_indexer(ov46[key])
    X = frame.values[idx]
    nact = (X > 6).sum(1)
    two = nact >= 2
    rg = np.array([np.ptp(row[row > 6]) for row in X[two]])
    sig = R[ds]['s_score'].std(0)[idx][two]
    r = spearmanr(rg, sig)
    print(f"   {ds:8s}: median n_active over all 46 = {np.median(nact):.0f}; "
          f"{two.sum()} with >= 2 active: Spearman(range, S-score sigma) = {r[0]:+.3f} (p = {r[1]:.2g})")

print("\n18. Most stable compounds: median rank under each definition (rank 1 = most selective)")
def med_ranks(ds, i):
    return "  ".join(f"{f} {np.median(R[ds][f][:, i]):5.1f}" for f in FAM)
for name in ('XL-228', 'PF-03814735'):
    i = Kl.index.get_loc(name)
    print(f"   Klaeger {name:12s}: n_active {na_k[i]:3d}, sigma {sig_k[i]:5.1f} | {med_ranks('Klaeger', i)}")
i = int(sig_d.argmin())
print(f"   Davis most stable (row {i}, Drug_Index {Dv.index[i]}): n_active {na_d[i]}, sigma {sig_d[i]:.1f} | {med_ranks('Davis', i)}")
for i in np.argsort(-na_d)[:5]:
    print(f"   Davis broad (row {i}): n_active {na_d[i]:3d}, sigma {sig_d[i]:5.1f} | {med_ranks('Davis', i)}")
broad = na_d > 40
sig_d_noratio = np.vstack([R['Davis'][f] for f in ('s_score', 'entropy', 'gini')]).std(0)
print(f"   Davis, {broad.sum()} compounds with > 40 active kinases: mean sigma {sig_d[broad].mean():.1f} pooled over all four "
      f"definitions, {sig_d_noratio[broad].mean():.1f} without the ratio's configurations")
wide = na_d > 20
print(f"   Davis, {wide.sum()} compounds with > 20 active kinases: mean sigma {sig_d[wide].mean():.1f} pooled, "
      f"{sig_d_noratio[wide].mean():.1f} without the ratio's configurations; the other {(~wide).sum()}: {sig_d[~wide].mean():.1f}")
print("   Davis, 8 most stable compounds (sigma, n_active): "
      + ", ".join(f"({sig_d[i]:.1f}, {na_d[i]})" for i in np.argsort(sig_d)[:8]))
print(f"   Klaeger compounds with more active kinases than PF-03814735: {(na_k > 57).sum()}")
narrow = (na_d >= 1) & (na_d <= 2)
print(f"   Davis, {narrow.sum()} compounds with 1-2 active kinases: mean sigma {sig_d[narrow].mean():.1f}; "
      f"median ranks " + "; ".join(f"row {i}: " + med_ranks('Davis', i) for i in np.where(narrow)[0]))


print("\n19. D4 under realistic panel expansion (selectivity score = -entropy; a rise means the")
print("    compound looks more selective after kinases it binds below threshold are added)")
import re


def cand_score(X, x0, T):
    w = (X > x0) * T * np.logaddexp(0.0, (X - x0) / T)
    t = w.sum(1, keepdims=True)
    p = w / np.where(t == 0, 1, t)
    return (p * np.log2(np.where(p > 0, p, 1.0))).sum(1)


D4SETS = (('Davis', Md, 5.0, 6.0, np.arange(5.0, 6.75, .25)),
          ('Klaeger', Mk, 5.0, 6.0, np.arange(5.0, 6.75, .25)),
          ('Anastassiadis', Ma, 20.0, 50.0, np.arange(0, 35, 5)),
          ('Metz', Mz, 4.0, 6.0, np.arange(4.0, 5.75, .25)))


def rise_rates(S, X, ok, x0, tau, betas):
    """Share of eligible compounds whose score rises from panel S to expanded panel X."""
    ent_r = {b: ((ent(X, b) - ent(S, b))[ok] > 1e-9).sum() for b in betas}
    cand_r = ((cand_score(X, x0, tau - x0) - cand_score(S, x0, tau - x0))[ok] > 1e-9).sum()
    return ent_r, cand_r, ok.sum()


for ds, M, x0, tau, betas in D4SETS:
    nk = M.shape[1]
    for P in (50, 100):
        rows = {}
        for mode in ('+1', '+5', '+10', '+30', 'all sub-threshold', 'one at tau* (worst case)'):
            rng = np.random.RandomState(42)
            tot_e = {b: 0 for b in betas}; tot_c = 0; tot_n = 0
            for _ in range(50):
                idx = rng.choice(nk, P, replace=False)
                rest = np.setdiff1d(np.arange(nk), idx)
                S = M[:, idx]
                scored = (S > tau).any(1)
                if mode.startswith('+'):
                    extra = rng.choice(rest, int(mode[1:]), replace=False)
                    X = np.hstack([S, M[:, extra]])
                    ok = scored & (M[:, extra] <= tau).all(1)
                    e, c, n = rise_rates(S, X, ok, x0, tau, betas)
                elif mode == 'all sub-threshold':
                    e = {b: 0 for b in betas}; c = 0; n = 0
                    for i in np.where(scored)[0]:
                        add = M[i, rest][M[i, rest] <= tau]
                        Xi = np.concatenate([S[i], add])[None, :]
                        Si = S[i][None, :]
                        ei, ci, ni = rise_rates(Si, Xi, np.array([True]), x0, tau, betas)
                        for b in betas: e[b] += ei[b]
                        c += ci; n += ni
                else:
                    X = np.hstack([S, np.full((len(M), 1), tau)])
                    e, c, n = rise_rates(S, X, scored, x0, tau, betas)
                for b in betas: tot_e[b] += e[b]
                tot_c += c; tot_n += n
            worst_b = max(betas, key=lambda b: tot_e[b])
            rows[mode] = (100 * tot_e[worst_b] / max(tot_n, 1), worst_b, 100 * tot_c / max(tot_n, 1), tot_n)
        print(f"   {ds:13s} {P:3d}-kinase panels (50 per size): entropy rises at its worst baseline | candidate rises | eligible pairs")
        for mode, (er, wb, cr, n) in rows.items():
            print(f"      {mode:26s}: {er:5.1f}% (beta = {wb:g}) | {cr:4.1f}% | {n}")

print("   Davis wild-type constructs, expanded by the mutant and phosphorylation-state constructs")
names = pd.read_csv(f"{HERE}/davis_proteins.csv").set_index('Protein_Index').loc[Dv.columns, 'Gene_Name']
variant = names.str.contains(r"phosphorylated|\([^)]*(?:[A-Z]\d+[A-Z]|del|ins)", regex=True).values
S, V = Md[:, ~variant], Md[:, variant]
print(f"      {(~variant).sum()} wild-type or domain constructs, {variant.sum()} mutant/phosphorylation variants")
betas = np.arange(5.0, 6.75, .25)
e = {b: 0 for b in betas}; c = 0; n = 0
for i in np.where((S > 6.0).any(1))[0]:
    add = V[i][V[i] <= 6.0]
    ei, ci, ni = rise_rates(S[i][None, :], np.concatenate([S[i], add])[None, :], np.array([True]), 5.0, 6.0, betas)
    for b in betas: e[b] += ei[b]
    c += ci; n += ni
print(f"      adding each compound's sub-threshold variants: entropy rises for "
      + ", ".join(f"beta={b:g}: {e[b]}" for b in betas) + f" of {n} compounds; candidate rises for {c}")

print("\n20. Zero-active compounds: where each definition ranks them, and the pooled rank sigma")
print("    with the undefined entropy/Gini configurations dropped (rank 1 = most selective)")
for ds, (M, s_thr, betas, floor, tau) in SETTINGS.items():
    act = (M > tau).sum(1) > 0
    if act.all():
        continue
    Rd = R[ds]
    undef = np.array([np.maximum(M - b, 0).sum(1) == 0 for b in betas])     # (n_beta, n_cpd)
    uz = undef[:, ~act]
    med = {f: np.median(np.median(Rd[f][:, ~act], 0)) for f in FAM}
    g_undef = f"{Rd['gini'][:, ~act][uz].mean():.0f}" if uz.any() else "n/a"
    g_def = Rd['gini'][:, ~act][~uz].mean()
    sig = {'all configurations': [], 'undefined entropy/Gini dropped': [], 'S-score and ratio only': []}
    for i in range(len(M)):
        d = ~undef[:, i]
        sig['all configurations'].append(np.concatenate([Rd[f][:, i] for f in FAM]).std())
        sig['undefined entropy/Gini dropped'].append(np.concatenate(
            [Rd['s_score'][:, i], Rd['entropy'][d, i], Rd['gini'][d, i], Rd['ratio'][:, i]]).std())
        sig['S-score and ratio only'].append(np.concatenate([Rd['s_score'][:, i], Rd['ratio'][:, i]]).std())
    print(f"   {ds} ({(~act).sum()} zero-active of {len(M)}): median of their median ranks: "
          + "  ".join(f"{f} {med[f]:.0f}" for f in FAM)
          + f"; entropy/Gini undefined in {100 * uz.mean():.0f}% of their configurations"
          + f" (mean Gini rank when undefined: {g_undef}; when defined: {g_def:.0f})")
    for k, v in sig.items():
        v = np.array(v)
        print(f"      pooled sigma, {k:31s}: zero-active {v[~act].mean():5.1f}  active {v[act].mean():5.1f}"
              f"  (Mann-Whitney p={mannwhitneyu(v[~act], v[act]).pvalue:.1g})")
