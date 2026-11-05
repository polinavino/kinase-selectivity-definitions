"""
Cross-definition correlation matrix under the ORIGINAL metric operationalizations
(revision for Molecular Informatics, Reviewer 1, major comment #1).

The manuscript computes entropy and Gini on a common baseline-subtracted
log-affinity scale so that all four definitions are comparable and so that the
baseline is an explicit, swept parameter (motivating D2). Reviewer 1 asked for
a robustness check: does the two-family clustering (S-score/entropy/Gini vs.
ratio) survive if entropy and Gini are instead computed in forms closer to
their original publications?

  * Entropy (Uitdehaag & Zaman 2011): association-constant weighting,
    p_i = Ka_i / sum_j Ka_j with Ka_i = 1/Kd_i, computed per compound with the
    floor-coded entries (no binding detected on Davis and Klaeger, not tested on
    Metz) excluded from the sum, and, as a sensitivity analysis, with them included
    at the floor value.
  * Gini (Graczyk 2010): computed directly on native percent-inhibition
    values (no baseline subtraction), which is only meaningful on the one
    percent-inhibition dataset, Anastassiadis.

S-score and ratio are already computed in their original form in the main
analysis (Methods, "Faithfulness to the original definitions"), so they are
reused unchanged here.
"""
import script_logging; script_logging.capture(__file__)
import os
os.environ['KMP_DUPLICATE_LIB_OK'] = 'TRUE'
import numpy as np
import pandas as pd
from scipy.stats import spearmanr, rankdata

HERE = os.path.dirname(os.path.abspath(__file__))


def to_ranks(s):
    return rankdata(-np.asarray(s), method='average')


def s_score(M, t):
    return -(M > t).astype(float).mean(axis=1)


def gini_selectivity(M, b=0.0):
    sh = np.maximum(M - b, 0)
    out = []
    for row in sh:
        rs = np.sort(row); n = len(rs); tot = rs.sum()
        out.append(0.0 if tot == 0 else (2 * np.sum(np.arange(1, n + 1) * rs)) / (n * tot) - (n + 1) / n)
    return np.array(out)


def ratio_selectivity(M, top_n, floor):
    out = []
    for row in M:
        s = np.sort(row)[::-1]
        off = s[top_n] if len(s) > top_n else floor
        out.append(s[0] - max(off, floor))
    return np.array(out)


def selectivity_entropy(M, b, eps=1e-10):
    sh = np.maximum(M - b, 0)
    rs = sh.sum(axis=1, keepdims=True)
    rs = np.where(rs == 0, eps, rs)
    p = sh / rs
    return -(p * np.where(p > 0, np.log2(np.where(p > 0, p, 1.0)), 0)).sum(axis=1)


def entropy_literal_dense(M, eps=1e-12):
    """Ka-weighted entropy over the full (complete) matrix -- Davis."""
    Ka = np.power(10.0, M - 9.0)
    p = Ka / Ka.sum(1, keepdims=True)
    return -(p * np.where(p > 0, np.log2(np.where(p > 0, p, 1.0)), 0)).sum(1)


def entropy_literal_sparse(pkd_rows, eps=1e-12):
    """Ka-weighted entropy per compound, each over only its own measured
    targets (a ragged list of pKd arrays, one per compound)."""
    out = []
    for pkds in pkd_rows:
        Ka = np.power(10.0, np.asarray(pkds) - 9.0)
        p = Ka / Ka.sum()
        out.append(-(p * np.where(p > 0, np.log2(np.where(p > 0, p, 1.0)), 0)).sum())
    return np.array(out)


DEFNS = ['s_score', 'entropy', 'gini', 'ratio']
S_THR = np.arange(5.5, 8.25, .25)
BASELINES = np.arange(5.0, 6.75, .25)


def corr_with_literal_entropy(M, ent_lit_scores, ratio_floor=5.0, s_thr=S_THR, baselines=BASELINES):
    # s_thr/baselines default to the pKd grids; Metz passes its own grids, as in
    # additional_datasets_analysis.py, so its S-score/Gini match the main analysis.
    ranks = {
        's_score': np.array([to_ranks(s_score(M, t)) for t in s_thr]),
        'gini':    np.array([to_ranks(gini_selectivity(M, b)) for b in baselines]),
        'ratio':   np.array([to_ranks(ratio_selectivity(M, k, ratio_floor)) for k in range(1, 6)]),
    }
    med = {d: np.median(ranks[d], axis=0) for d in ['s_score', 'gini', 'ratio']}
    med['entropy'] = to_ranks(-ent_lit_scores)   # rank 1 = most selective = lowest entropy
    return np.array([[spearmanr(med[a], med[b])[0] for b in DEFNS] for a in DEFNS])


def print_matrix(title, corr):
    print(f"\n{title}")
    print(f"{'':8s}" + "".join(f"{d:8s}" for d in DEFNS))
    for i, d in enumerate(DEFNS):
        print(f"{d:8s}" + "".join(f"{corr[i, j]:8.3f}" for j in range(4)))


# --------------------------------- Davis ------------------------------------
aff = pd.read_csv(os.path.join(HERE, "davis_affinity.csv"))
M_davis = aff.pivot(index='Drug_Index', columns='Protein_Index', values='Affinity').fillna(5.0).values
ent_lit_davis = entropy_literal_dense(M_davis)
corr_davis = corr_with_literal_entropy(M_davis, ent_lit_davis)
print_matrix("Davis (n=68): manuscript S/Gini/ratio vs. literal Ka-weighted entropy "
             "(complete matrix, floor-coded entries included)", corr_davis)
# Davis codes non-binding pairs at the floor (pKd 5.0); a floor-coded entry is a
# censored value, not a measured affinity. Sensitivity: exclude them per compound.
ent_lit_davis_meas = entropy_literal_sparse([r[r > 5.0] for r in M_davis])
print(f"Davis floor-coded entries: {100 * (M_davis == 5.0).mean():.2f}%; Ka entropy with vs without them: "
      f"Spearman {spearmanr(ent_lit_davis, ent_lit_davis_meas)[0]:.3f}")
print_matrix("Davis (n=68): manuscript S/Gini/ratio vs. literal Ka-weighted entropy "
             "(floor-coded entries excluded)", corr_with_literal_entropy(M_davis, ent_lit_davis_meas))

# -------------------------------- Klaeger ------------------------------------
# curated protein-kinase panel (curate_panels.py); entries above the floor are the
# measured high-confidence affinities, the rest are pairs with no detected binding
M_klaeger = pd.read_csv(os.path.join(HERE, "klaeger_matrix.csv"), index_col=0)
pkd_rows = [r[r > 5.0] for r in M_klaeger.values]
# a compound with no measured curated target has an undefined entropy, scored 0 as elsewhere
ent_lit_klaeger = np.array([entropy_literal_sparse([r])[0] if len(r) else 0.0 for r in pkd_rows])
print(f"Klaeger compounds with no measured curated target (entropy undefined, set to 0): "
      f"{sum(len(r) == 0 for r in pkd_rows)}")
corr_klaeger = corr_with_literal_entropy(M_klaeger.values, ent_lit_klaeger)
print_matrix("Klaeger (n=222): manuscript S/Gini/ratio vs. literal Ka-weighted entropy "
             "(floor-coded entries excluded)", corr_klaeger)
print_matrix("Klaeger (n=222): manuscript S/Gini/ratio vs. literal Ka-weighted entropy "
             "(floor-coded entries included)",
             corr_with_literal_entropy(M_klaeger.values, entropy_literal_dense(M_klaeger.values)))
act_k = (M_klaeger.values > 6.0).any(1)
print_matrix(f"Klaeger, the {int(act_k.sum())} compounds with an active kinase: manuscript S/Gini/ratio vs. "
             "literal Ka-weighted entropy (floor-coded entries excluded)",
             corr_with_literal_entropy(M_klaeger.values[act_k], ent_lit_klaeger[act_k]))
print_matrix(f"Klaeger, the {int(act_k.sum())} compounds with an active kinase: manuscript S/Gini/ratio vs. "
             "literal Ka-weighted entropy (floor-coded entries included)",
             corr_with_literal_entropy(M_klaeger.values[act_k], entropy_literal_dense(M_klaeger.values[act_k])))

# ---------------------------------- Metz --------------------------------------
mz = pd.read_excel(os.path.join(HERE, "metz.xls"), sheet_name=0, header=0)
meta = ['Cmpd_ID', 'PUBCHEM_SID', 'Canonical_Smiles', 'External_Cmpd_ID', 'External_Source',
        'Cluster', 'ClusterSize', 'Cluster_MCSS', 'Molecular_Weight', 'ALogP',
        'Num_H_Acceptors', 'Num_H_Donors', 'tPSA', 'Promiscuity_1uM']
kin_cols = [c for c in mz.columns if c not in meta]
sub = mz[kin_cols].apply(pd.to_numeric, errors='coerce')
coverage = (~np.isnan(sub.values)).sum(axis=1)
keep = coverage >= 50
sub_kept = sub.values[keep]
M_metz = np.nan_to_num(sub_kept, nan=4.0)
ent_lit_metz = entropy_literal_sparse([row[~np.isnan(row)] for row in sub_kept])
corr_metz = corr_with_literal_entropy(M_metz, ent_lit_metz, ratio_floor=4.0,
                                      s_thr=np.arange(5.0, 7.75, .25),
                                      baselines=np.arange(4.0, 5.75, .25))
print_matrix("Metz (n=704): manuscript S/Gini/ratio vs. literal Ka-weighted entropy "
             "(floor-imputed entries excluded: measured entries only)", corr_metz)
print_matrix("Metz (n=704): manuscript S/Gini/ratio vs. literal Ka-weighted entropy "
             "(floor-imputed entries included)",
             corr_with_literal_entropy(M_metz, entropy_literal_dense(M_metz), ratio_floor=4.0,
                                       s_thr=np.arange(5.0, 7.75, .25), baselines=np.arange(4.0, 5.75, .25)))

# ----------------------------- Anastassiadis -----------------------------------
csv = os.path.join(HERE, "anastassiadis_matrix.csv")
M_ana = pd.read_csv(csv, index_col=0).values   # percent inhibition, 0-100
ranks_ana = {
    's_score': np.array([to_ranks(s_score(M_ana, t)) for t in np.arange(50, 95, 5)]),
    'entropy': np.array([to_ranks(-selectivity_entropy(M_ana, b)) for b in np.arange(0, 35, 5)]),
    'ratio':   np.array([to_ranks(ratio_selectivity(M_ana, k, 0.0)) for k in range(1, 6)]),
}
med_ana = {d: np.median(ranks_ana[d], axis=0) for d in ['s_score', 'entropy', 'ratio']}
gini_native = gini_selectivity(M_ana, b=0.0)   # Graczyk's native %inhibition Gini, no baseline
med_ana['gini'] = to_ranks(gini_native)
corr_ana = np.array([[spearmanr(med_ana[a], med_ana[b])[0] for b in DEFNS] for a in DEFNS])
print_matrix("Anastassiadis (n=178): manuscript S/entropy/ratio vs. native %inhibition Gini "
             "(no baseline subtraction)", corr_ana)

print("\nDone: original-operationalization robustness check complete for all four datasets.")
