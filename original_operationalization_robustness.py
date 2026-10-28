"""
Cross-definition correlation matrix under the ORIGINAL metric operationalizations
(revision for Molecular Informatics, Reviewer 1, major comment #1).

The manuscript computes entropy and Gini on a common baseline-subtracted
log-affinity scale so that all four definitions are comparable and so that the
baseline is an explicit, swept parameter (motivating D3). Reviewer 1 asked for
a robustness check: does the two-family clustering (S-score/entropy/Gini vs.
ratio) survive if entropy and Gini are instead computed in forms closer to
their original publications?

  * Entropy (Uitdehaag & Zaman 2011): association-constant weighting,
    p_i = Ka_i / sum_j Ka_j with Ka_i = 1/Kd_i, computed PER COMPOUND ONLY
    OVER ITS MEASURED TARGETS (Davis is a complete matrix, so "measured" is
    every cell; Klaeger and Metz are sparse, so unmeasured kinase-compound
    pairs are excluded from the sum rather than treated as floor measurements).
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
    return -(p * np.where(p > 0, np.log2(p + eps), 0)).sum(axis=1)


def entropy_literal_dense(M, eps=1e-12):
    """Ka-weighted entropy over the full (complete) matrix -- Davis."""
    Ka = np.power(10.0, M - 9.0)
    p = Ka / Ka.sum(1, keepdims=True)
    return -(p * np.where(p > 0, np.log2(p + eps), 0)).sum(1)


def entropy_literal_sparse(pkd_rows, eps=1e-12):
    """Ka-weighted entropy per compound, each over only its own measured
    targets (a ragged list of pKd arrays, one per compound)."""
    out = []
    for pkds in pkd_rows:
        Ka = np.power(10.0, np.asarray(pkds) - 9.0)
        p = Ka / Ka.sum()
        out.append(-(p * np.where(p > 0, np.log2(p + eps), 0)).sum())
    return np.array(out)


DEFNS = ['s_score', 'entropy', 'gini', 'ratio']
S_THR = np.arange(5.5, 8.25, .25)
BASELINES = np.arange(5.0, 6.75, .25)


def corr_with_literal_entropy(M, ent_lit_scores, ratio_floor=5.0):
    ranks = {
        's_score': np.array([to_ranks(s_score(M, t)) for t in S_THR]),
        'gini':    np.array([to_ranks(gini_selectivity(M, b)) for b in BASELINES]),
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
             "(complete matrix, all entries 'measured')", corr_davis)

# -------------------------------- Klaeger ------------------------------------
raw = pd.read_excel(os.path.join(HERE, "aan4368_Table_S2.xlsx"), sheet_name='Kinobeads')
raw['pKd'] = 9 - np.log10(raw['Apparent Kd'])
raw = raw.dropna(subset=['pKd'])
high = raw[raw['Target Classification'] == 'High confidence']
M_klaeger = high.groupby(['Drug', 'Gene Name'])['pKd'].mean().unstack(fill_value=5.0)
drugs = M_klaeger.index.tolist()
by_drug = high.groupby('Drug')['pKd'].apply(list)
pkd_rows = [by_drug.get(d, []) for d in drugs]
assert all(len(r) > 0 for r in pkd_rows), "every Klaeger drug must have >=1 measured target"
ent_lit_klaeger = entropy_literal_sparse(pkd_rows)
corr_klaeger = corr_with_literal_entropy(M_klaeger.values, ent_lit_klaeger)
print_matrix("Klaeger (n=222): manuscript S/Gini/ratio vs. literal Ka-weighted entropy "
             "(measured high-confidence entries only, per compound)", corr_klaeger)

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
corr_metz = corr_with_literal_entropy(M_metz, ent_lit_metz, ratio_floor=4.0)
print_matrix("Metz (n=704): manuscript S/Gini/ratio vs. literal Ka-weighted entropy "
             "(measured entries only, per compound)", corr_metz)

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
