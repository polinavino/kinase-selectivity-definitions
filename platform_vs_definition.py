"""
Platform vs. definition: does the assay platform (Davis competition-binding vs.
Klaeger chemoproteomics) change compound rankings more or less than the choice
of selectivity definition? (Revision for Molecular Informatics, Reviewer 1,
major comment: "benchmark definitional instability against the ranking
variability introduced by assay platform.")

Uses the 46 compounds profiled in both Davis and Klaeger (resolved by PubChem
CID -> Klaeger drug-name match, davis_klaeger_overlap.py; this supersedes the
manuscript's previous, uninvestigated claim of "at most 19" shared compounds).
For this fixed set of 46 compounds, each selectivity definition is computed
separately on the Davis profile and on the Klaeger profile, and compounds are
ranked within this n=46 set (not the full 68/222-compound set) so that
platform and definition effects are compared on identical footing.
"""
import script_logging; script_logging.capture(__file__)
import os
import numpy as np
import pandas as pd
from scipy.stats import spearmanr, rankdata

HERE = os.path.dirname(os.path.abspath(__file__))


def to_ranks(s):
    return rankdata(-np.asarray(s), method='average')


def s_score(M, t):
    return -(M > t).astype(float).mean(axis=1)


def selectivity_entropy(M, b, eps=1e-10):
    sh = np.maximum(M - b, 0)
    rs = sh.sum(axis=1, keepdims=True)
    rs = np.where(rs == 0, eps, rs)
    p = sh / rs
    return -(p * np.where(p > 0, np.log2(np.where(p > 0, p, 1.0)), 0)).sum(axis=1)


def gini_selectivity(M, b):
    sh = np.maximum(M - b, 0)
    out = []
    for row in sh:
        rs = np.sort(row); n = len(rs); tot = rs.sum()
        out.append(0.0 if tot == 0 else (2 * np.sum(np.arange(1, n + 1) * rs)) / (n * tot) - (n + 1) / n)
    return np.array(out)


def ratio_selectivity(M, top_n, floor=5.0):
    out = []
    for row in M:
        s = np.sort(row)[::-1]
        off = s[top_n] if len(s) > top_n else floor
        out.append(s[0] - max(off, floor))
    return np.array(out)


DEFNS = ['s_score', 'entropy', 'gini', 'ratio']
S_THR = np.arange(5.5, 8.25, .25)
BASELINES = np.arange(5.0, 6.75, .25)


def median_ranks(M):
    ranks = {
        's_score': np.array([to_ranks(s_score(M, t)) for t in S_THR]),
        'entropy': np.array([to_ranks(-selectivity_entropy(M, b)) for b in BASELINES]),
        'gini':    np.array([to_ranks(gini_selectivity(M, b)) for b in BASELINES]),
        'ratio':   np.array([to_ranks(ratio_selectivity(M, k)) for k in range(1, 6)]),
    }
    return {d: np.median(ranks[d], axis=0) for d in DEFNS}


overlap = pd.read_csv(os.path.join(HERE, "davis_klaeger_overlap.csv"))

aff = pd.read_csv(os.path.join(HERE, "davis_affinity.csv"))
davis_full = aff.pivot(index='Drug_Index', columns='Protein_Index', values='Affinity').fillna(5.0)
klaeger_full = pd.read_csv(os.path.join(HERE, "klaeger_matrix.csv"), index_col=0)

davis_M = davis_full.loc[overlap['Davis_Drug_Index'].values].values
klaeger_M = klaeger_full.loc[overlap['Klaeger_Name'].values].values
n = len(overlap)
print(f"Shared Davis/Klaeger compounds (PubChem CID match): n = {n}")
print("(supersedes the manuscript's previous, unverified claim of 'at most 19')")

med_davis = median_ranks(davis_M)
med_klaeger = median_ranks(klaeger_M)

print(f"\nPlatform effect: Spearman r between Davis-rank and Klaeger-rank, same definition, "
      f"same {n} compounds:")
for d in DEFNS:
    r, p = spearmanr(med_davis[d], med_klaeger[d])
    print(f"  {d:8s}: r = {r:+.3f} (p = {p:.1e})")

print(f"\nDefinition effect: Spearman r between definitions, within platform, same {n} compounds:")
for platform, med in [("Davis", med_davis), ("Klaeger", med_klaeger)]:
    print(f"  {platform}:")
    for i in range(4):
        for j in range(i + 1, 4):
            r, _ = spearmanr(med[DEFNS[i]], med[DEFNS[j]])
            print(f"    {DEFNS[i]:8s} vs {DEFNS[j]:8s}: r = {r:+.3f}")

platform_rs = [spearmanr(med_davis[d], med_klaeger[d])[0] for d in DEFNS]
defn_rs = []
for med in [med_davis, med_klaeger]:
    for i in range(4):
        for j in range(i + 1, 4):
            defn_rs.append(spearmanr(med[DEFNS[i]], med[DEFNS[j]])[0])

print(f"\nSummary: platform-effect r range = {min(platform_rs):.2f}-{max(platform_rs):.2f} "
      f"(mean {np.mean(platform_rs):.2f})")
print(f"         definition-effect r range = {min(defn_rs):.2f}-{max(defn_rs):.2f} "
      f"(mean {np.mean(defn_rs):.2f})")
print("Lower r = more disagreement. Comparing the two ranges shows which source of "
      "ranking disagreement -- assay platform or scoring definition -- is larger on "
      "this shared-compound set.")

# ---------------------------------------------------------------------------
# The two platforms profile different targets (433 Davis constructs, 343 Klaeger
# proteins), so the comparison above mixes assay and panel. Repeat it on the
# genes both panels share: Davis constructs collapsed to genes (best affinity per
# gene), Davis names mapped to HGNC symbols, Klaeger columns matched on any
# ';'-separated symbol.
ALIAS = {
    'AMPK-alpha1': 'PRKAA1', 'AMPK-alpha2': 'PRKAA2', 'ARK5': 'NUAK1', 'ASK1': 'MAP3K5', 'ASK2': 'MAP3K6',
    'BIKE': 'BMP2K', 'BRK': 'PTK6', 'CDC2L5': 'CDK13', 'CTK': 'MATK', 'DLK': 'MAP3K12', 'DMPK2': 'CDC42BPG',
    'ERK1': 'MAPK3', 'ERK2': 'MAPK1', 'ERK3': 'MAPK6', 'ERK4': 'MAPK4', 'ERK5': 'MAPK7', 'ERK8': 'MAPK15',
    'FAK': 'PTK2', 'GCN2': 'EIF2AK4', 'HPK1': 'MAP4K1', 'IKK-alpha': 'CHUK', 'IKK-beta': 'IKBKB',
    'IKK-epsilon': 'IKBKE', 'JNK1': 'MAPK8', 'JNK2': 'MAPK9', 'JNK3': 'MAPK10', 'LKB1': 'STK11',
    'LOK': 'STK10', 'LZK': 'MAP3K13', 'MEK1': 'MAP2K1', 'MEK2': 'MAP2K2', 'MEK3': 'MAP2K3', 'MEK4': 'MAP2K4',
    'MEK5': 'MAP2K5', 'MEK6': 'MAP2K6', 'MKK7': 'MAP2K7', 'MINK': 'MINK1', 'MLK1': 'MAP3K9',
    'MLK2': 'MAP3K10', 'MLK3': 'MAP3K11', 'MRCKA': 'CDC42BPA', 'MRCKB': 'CDC42BPB', 'MST1': 'STK4',
    'MST2': 'STK3', 'MST3': 'STK24', 'MST4': 'STK26', 'NDR1': 'STK38', 'NDR2': 'STK38L', 'PCTK1': 'CDK16',
    'PCTK2': 'CDK17', 'PCTK3': 'CDK18', 'PFTAIRE2': 'CDK15', 'PFTK1': 'CDK14', 'PAK7': 'PAK5',
    'PKAC-alpha': 'PRKACA', 'PKAC-beta': 'PRKACB', 'PRKR': 'EIF2AK2', 'PRP4': 'PRPF4B', 'PYK2': 'PTK2B',
    'QSK': 'SIK3', 'RSK1': 'RPS6KA1', 'RSK2': 'RPS6KA3', 'RSK3': 'RPS6KA2', 'RSK4': 'RPS6KA6',
    'S6K1': 'RPS6KB1', 'SIK': 'SIK1', 'SNARK': 'NUAK2', 'TAK1': 'MAP3K7', 'TIE2': 'TEK', 'TRKA': 'NTRK1',
    'TRKB': 'NTRK2', 'TRKC': 'NTRK3', 'VEGFR2': 'KDR', 'YANK1': 'STK32A', 'YES': 'YES1', 'YSK1': 'STK25',
    'YSK4': 'MAP3K19', 'p38-alpha': 'MAPK14', 'p38-beta': 'MAPK11', 'p38-delta': 'MAPK13',
    'p38-gamma': 'MAPK12',
}
prot = pd.read_csv(os.path.join(HERE, "davis_proteins.csv")).set_index("Protein_Index")["Gene_Name"]
names = prot.reindex(davis_full.columns)
gene = names.str.replace(r"\(.*$", "", regex=True).str.replace(r"-(phosphorylated|nonphosphorylated)$", "", regex=True)
gene = gene.map(lambda g: ALIAS.get(g, g))
davis_gene = davis_full.T.groupby(gene.values).max().T
kcol = {sym: c for c in klaeger_full.columns for sym in c.split(';')}
common = sorted(g for g in davis_gene.columns if g in kcol)
davis_C = davis_gene.loc[overlap['Davis_Drug_Index'].values, common].values
klaeger_C = klaeger_full.loc[overlap['Klaeger_Name'].values, [kcol[g] for g in common]].values
print(f"\nShared panel: {len(common)} genes profiled on both platforms "
      f"(of {davis_gene.shape[1]} Davis genes and {klaeger_full.shape[1]} Klaeger targets)")
print(f"  compounds with no active target on the shared panel: Davis {int((~(davis_C > 6).any(1)).sum())}, "
      f"Klaeger {int((~(klaeger_C > 6).any(1)).sum())}")
mdC, mkC = median_ranks(davis_C), median_ranks(klaeger_C)
plat_C = [spearmanr(mdC[d], mkC[d])[0] for d in DEFNS]
defn_C = [spearmanr(med[DEFNS[i]], med[DEFNS[j]])[0] for med in (mdC, mkC)
          for i in range(4) for j in range(i + 1, 4)]
for d, r in zip(DEFNS, plat_C):
    print(f"  platform effect, {d:8s}: r = {r:+.3f} (p = {spearmanr(mdC[d], mkC[d])[1]:.1e})")
print(f"Summary on the shared panel: platform-effect r range = {min(plat_C):.2f}-{max(plat_C):.2f} "
      f"(mean {np.mean(plat_C):.2f}); definition-effect r range = {min(defn_C):.2f}-{max(defn_C):.2f} "
      f"(mean {np.mean(defn_C):.2f})")
