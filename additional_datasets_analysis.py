"""
Additional-dataset robustness analysis (revision for Molecular Informatics).

Replicates the implementation-dependent agreement and rank-instability analyses of the main
paper on two further kinase profiling datasets that use different assay
technologies from Davis (pKd binding) and Klaeger (pKd chemoproteomics):

  * Anastassiadis et al., Nat. Biotechnol. 2011 -- 178 inhibitors x 300 kinases
    (291 after curate_panels.py),
    single-concentration (0.5 uM) functional activity assay reported as percent
    residual kinase activity.  Converted here to percent inhibition
    (100 - residual, clipped to [0,100]); higher = more potent.
  * Metz et al., Nat. Chem. Biol. 2011 -- pKi affinity data; restricted to the
    704 compounds profiled against >= 50 of the 172 kinases.  Same -log10
    affinity scale as pKd.

Raw supplementary files (downloaded programmatically; URLs below) are converted
to compound x kinase potency matrices and written next to this script.

  Anastassiadis: https://static-content.springer.com/esm/art%3A10.1038%2Fnbt.2017/MediaObjects/41587_2011_BFnbt2017_MOESM23_ESM.xls
  Metz:          https://static-content.springer.com/esm/art%3A10.1038%2Fnchembio.530/MediaObjects/41589_2011_BFnchembio530_MOESM137_ESM.xls

Each metric is computed on the native potency scale of its dataset with a
scale-appropriate parameter sweep. Only the resulting within-dataset rankings
are compared across datasets. A nonlinear transformation of the raw readout
can change a score's ranking, particularly for entropy and Gini weights.
"""
import script_logging; script_logging.capture(__file__)
import os
os.environ['KMP_DUPLICATE_LIB_OK'] = 'TRUE'
import sys
import numpy as np
import pandas as pd
from scipy.stats import spearmanr, rankdata
import matplotlib.pyplot as plt
from panel_curation import load_anastassiadis_panel

HERE = os.path.dirname(os.path.abspath(__file__))

# --------------------------- selectivity metrics ---------------------------
def s_score(M, threshold):
    return -(M > threshold).astype(float).mean(axis=1)

def selectivity_entropy(M, baseline, eps=1e-10):
    sh = np.maximum(M - baseline, 0)
    rs = sh.sum(axis=1, keepdims=True)
    rs = np.where(rs == 0, eps, rs)
    p = sh / rs
    return -(-(p * np.where(p > 0, np.log2(np.where(p > 0, p, 1.0)), 0)).sum(axis=1))

def gini_selectivity(M, baseline):
    sh = np.maximum(M - baseline, 0)
    out = []
    for row in sh:
        rs = np.sort(row); n = len(rs); tot = rs.sum()
        out.append(0.0 if tot == 0
                   else (2 * np.sum(np.arange(1, n + 1) * rs)) / (n * tot) - (n + 1) / n)
    return np.array(out)

def ratio_selectivity(M, top_n, floor):
    out = []
    for row in M:
        s = np.sort(row)[::-1]
        off = s[top_n] if len(s) > top_n else floor
        out.append(s[0] - max(off, floor))
    return np.array(out)

def to_ranks(scores):
    return rankdata(-np.asarray(scores), method='average')

DEFNS = ['s_score', 'entropy', 'gini', 'ratio']

def analyze(M, s_thr, baselines, ratio_floor, active_thr, label):
    ranks = {
        's_score': np.array([to_ranks(s_score(M, t)) for t in s_thr]),
        'entropy': np.array([to_ranks(selectivity_entropy(M, b)) for b in baselines]),
        'gini':    np.array([to_ranks(gini_selectivity(M, b)) for b in baselines]),
        'ratio':   np.array([to_ranks(ratio_selectivity(M, k, ratio_floor)) for k in range(1, 6)]),
    }
    all_ranks = np.vstack(list(ranks.values()))
    rank_std = all_ranks.std(axis=0)
    # feature associations: compounds with an active kinase, ranked among themselves (as SI Table S4)
    act_c = (M > active_thr).any(1)
    Ma_ = M[act_c]
    ranks_a = {
        's_score': np.array([to_ranks(s_score(Ma_, t)) for t in s_thr]),
        'entropy': np.array([to_ranks(selectivity_entropy(Ma_, b)) for b in baselines]),
        'gini':    np.array([to_ranks(gini_selectivity(Ma_, b)) for b in baselines]),
        'ratio':   np.array([to_ranks(ratio_selectivity(Ma_, k, ratio_floor)) for k in range(1, 6)]),
    }
    instab = {d: np.full(M.shape[0], np.nan) for d in DEFNS}
    for d in DEFNS:
        instab[d][act_c] = ranks_a[d].std(axis=0)
    med_a = {d: np.median(ranks_a[d], axis=0) for d in DEFNS}
    corr_act = np.array([[spearmanr(med_a[a], med_a[b])[0] for b in DEFNS] for a in DEFNS])
    med = {d: np.median(ranks[d], axis=0) for d in DEFNS}
    corr = np.array([[spearmanr(med[a], med[b])[0] for b in DEFNS] for a in DEFNS])

    n_active = (M > active_thr).sum(axis=1)
    has_active = n_active > 0
    feats = {}
    for i in range(M.shape[0]):
        prof = M[i]; act = prof[prof > active_thr]; s = np.sort(prof)[::-1]
        feats[i] = {'n_active': len(act), 'top1_top2_gap': s[0] - s[1],
                    'active_range': (act.max() - act.min()) if len(act) > 1 else 0}
    fdf = pd.DataFrame(feats).T
    mask = has_active

    def fc(feature, defn):
        return spearmanr(fdf[feature].values[mask], instab[defn][mask])

    return dict(label=label, n_cpd=M.shape[0], n_kin=M.shape[1], corr=corr,
                corr_act=corr_act, n_act=int(act_c.sum()),
                rank_std=rank_std, has_active=has_active, n_active=n_active,
                fc=fc, med=med)

# ------------------------------- datasets ----------------------------------
def load_davis():
    aff = pd.read_csv(os.path.join(HERE, "davis_affinity.csv"))
    M = aff.pivot(index='Drug_Index', columns='Protein_Index', values='Affinity').fillna(5.0).values
    return analyze(M, np.arange(5.5, 8.25, .25), np.arange(5.0, 6.75, .25), 5.0, 6.0,
                   "Davis (pKd binding)")

def load_klaeger():
    km = pd.read_csv(os.path.join(HERE, "klaeger_matrix.csv"), index_col=0)
    return analyze(km.values, np.arange(5.5, 8.25, .25), np.arange(5.0, 6.75, .25), 5.0, 6.0,
                   "Klaeger (pKd chemoproteomics)")

def load_anastassiadis():
    M = load_anastassiadis_panel(HERE).values
    return analyze(M, np.arange(50, 95, 5), np.arange(0, 35, 5), 0.0, 50.0,
                   "Anastassiadis (% inhibition)")

def load_metz():
    csv = os.path.join(HERE, "metz_matrix.csv")
    if os.path.exists(csv):
        M = pd.read_csv(csv, index_col=0).values
    else:
        mz = pd.read_excel(os.path.join(HERE, "metz.xls"), sheet_name=0, header=0)
        meta = ['Cmpd_ID', 'PUBCHEM_SID', 'Canonical_Smiles', 'External_Cmpd_ID', 'External_Source',
                'Cluster', 'ClusterSize', 'Cluster_MCSS', 'Molecular_Weight', 'ALogP',
                'Num_H_Acceptors', 'Num_H_Donors', 'tPSA', 'Promiscuity_1uM']
        kin_cols = [c for c in mz.columns if c not in meta]
        sub = mz[kin_cols].apply(pd.to_numeric, errors='coerce')
        keep = (~np.isnan(sub.values)).sum(axis=1) >= 50
        M = np.nan_to_num(sub.values[keep], nan=4.0)
        pd.DataFrame(M, index=mz['Cmpd_ID'].values[keep], columns=kin_cols).to_csv(csv)
    return analyze(M, np.arange(5.0, 7.75, .25), np.arange(4.0, 5.75, .25), 4.0, 6.0,
                   "Metz (pKi)")

# --------------------------------- run -------------------------------------
results = [load_davis(), load_klaeger(), load_anastassiadis(), load_metz()]

rows = []
pvals = []
print(f"{'Dataset':32s}{'n_cpd':>7s}{'n_kin':>7s}{'distr.r':>10s}{'ratio-distr.r':>15s}"
      f"{'gap~ratio':>12s}{'nact~entropy':>14s}")
for r in results:
    c = r['corr']
    # within distribution family = mean of (s,ent),(s,gini),(ent,gini)
    distr = np.mean([c[0, 1], c[0, 2], c[1, 2]])
    # ratio vs each distribution metric
    rd_lo, rd_hi = c[3, :3].min(), c[3, :3].max()
    gap_r, gap_p = r['fc']('top1_top2_gap', 'ratio')
    na_r, na_p = r['fc']('n_active', 'entropy')
    print(f"{r['label']:32s}{r['n_cpd']:>7d}{r['n_kin']:>7d}{distr:>10.3f}"
          f"{rd_lo:>8.2f}-{rd_hi:<5.2f}{gap_r:>10.4f}  {na_r:>10.4f}")
    pvals.append(f"{r['label'].split(' ')[0]} gap~ratio p = {gap_p:.2g}, nact~entropy p = {na_p:.2g}")
    rows.append(dict(dataset=r['label'], n_cpd=r['n_cpd'], n_kin=r['n_kin'],
                     within_distr_r=round(distr, 3),
                     ratio_vs_distr_min=round(rd_lo, 3), ratio_vs_distr_max=round(rd_hi, 3),
                     gap_ratio_r=round(gap_r, 4), gap_ratio_p=gap_p,
                     nactive_entropy_r=round(na_r, 4), nactive_entropy_p=na_p,
                     zero_active=int((~r['has_active']).sum()),
                     rankstd_zeroactive=round(float(r['rank_std'][~r['has_active']].mean()), 1)
                         if (~r['has_active']).any() else np.nan,
                     rankstd_active=round(float(r['rank_std'][r['has_active']].mean()), 1)))
print("p-values: " + "; ".join(pvals))
pd.DataFrame(rows).to_csv(os.path.join(HERE, "cross_dataset_summary.csv"), index=False)
print("\nSaved cross_dataset_summary.csv")

# ------------------------------- figure ------------------------------------
fig, axes = plt.subplots(4, 2, figsize=(9, 10.8))
labels = ['S', 'Ent', 'Gini', 'Ratio']
for row, r in enumerate(results):
    for col, (key, what) in enumerate((("corr_act", "active"), ("corr", "all compounds"))):
        ax = axes[row, col]
        C = r[key]
        im = ax.imshow(C, vmin=0, vmax=1, cmap='RdYlGn')
        ax.set_xticks(range(4)); ax.set_yticks(range(4))
        ax.set_xticklabels(labels, fontsize=11); ax.set_yticklabels(labels, fontsize=11)
        for i in range(4):
            for j in range(4):
                ax.text(j, i, f"{C[i, j]:.2f}", ha='center', va='center', fontsize=11)
        n = r['n_act'] if key == "corr_act" else r['n_cpd']
        ax.set_title(f"{r['label'].split(' ')[0]} ({n}×{r['n_kin']})", fontsize=11)
fig.subplots_adjust(left=0.10, right=0.88, bottom=0.04, top=0.88, hspace=0.58, wspace=0.4)
color_axis = fig.add_axes([0.91, 0.30, 0.015, 0.40])
fig.colorbar(im, cax=color_axis, label='Spearman r')
fig.suptitle("Spearman correlations between median rankings", fontsize=13, y=0.995)
fig.text(0.31, 0.95, "Active compounds", ha='center', fontsize=13)
fig.text(0.70, 0.95, "All compounds", ha='center', fontsize=13)
out = os.path.join(HERE, "cross_dataset_correlations.png")
script_logging.save_figure(fig, out, bbox_inches='tight')
print(f"Saved {out}")
