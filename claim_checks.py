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
 11. Within-family correlation range on each dataset (Introduction, Results).
 12. Compounds shared by Klaeger and Anastassiadis, by name (Discussion).
 13. Zero-active vs active rank instability per definition, on every dataset with
     zero-active compounds (Results, replication).
 14. Klaeger Gini panel-size p* when only compounds whose Gini is undefined on the
     subpanel or the full panel are excluded (Results, gated panel size).
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
    lo, hi = min(lo, min(w)), max(hi, max(w))
    print(f"   {ds:13s}: {min(w):.3f}-{max(w):.3f}")
print(f"   across all four datasets: {lo:.3f}-{hi:.3f}")

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

print("\n14. Klaeger Gini p* (beta = x0 = 5.0), excluding only compounds whose Gini is undefined")
rng = np.random.RandomState(42)
ref = gini(Mk, 5.0)
curve = []
for ps in list(range(50, Mk.shape[1], 30)) + [Mk.shape[1]]:
    rs = []
    for _ in range(50):
        idx = rng.choice(Mk.shape[1], ps, replace=False)
        keep = (Mk[:, idx] > 5.0).any(1) & (Mk > 5.0).any(1)
        rs.append(spearmanr(to_ranks(ref[keep]), to_ranks(gini(Mk[:, idx], 5.0)[keep]))[0])
    curve.append((ps, np.mean(rs)))
print(f"   p* = {next(p for p, r in curve if r > 0.90)}; " + " ".join(f"{p}:{r:.3f}" for p, r in curve[:5]))
