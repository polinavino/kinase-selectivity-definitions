"""
Klaeger: the released 343-target panel against the curated protein-kinase panel.

The 343 Klaeger targets are the proteins the Kinobeads pulldown quantified, not
only kinases: they include cyclins and other subunits of kinase complexes,
co-captured with their kinase partners, as well as small-molecule and lipid
kinases and unrelated proteins. A CDK hit can therefore also count its cyclin.
This script groups the columns that curate_panels.py removes, by kind, and
recomputes the main Klaeger results on the released panel, printed beside
the curated panel (217 protein kinases, no protein-group columns; curate_panels.py),
which is the panel every other analysis uses.
"""
import script_logging; script_logging.capture(__file__)
import os
os.environ['KMP_DUPLICATE_LIB_OK'] = 'TRUE'
import numpy as np
import pandas as pd
from scipy.stats import spearmanr, rankdata, skew, mannwhitneyu

X0, TAU, T = 5.0, 6.0, 1.0

# Columns that are not protein kinases, by category.
NON_PK = {
    "kinase-complex subunits (cyclins, regulatory and adaptor subunits)": [
        "CCNA2", "CCNB1", "CCNB2", "CCNE1", "CCNH", "CCNI", "CCNK", "CCNT1", "CCNT2",
        "FAM58A;FAM58BP", "CKS1B", "CABLES1", "MNAT1", "CSNK2B", "PRKAB1", "PRKAB2",
        "PRKAG1", "PRKAG2", "PRKAR2A", "PHKA2", "CAB39", "MOB1A;MOB1B", "SAV1",
        "PDCD10", "RASSF2", "RASSF3", "RASSF5", "TPRKB", "LAGE3", "OSGEP", "INCENP",
        "GSKIP", "DCAF7", "FAM83A", "FAM83B", "AZI2", "TANK", "TBKBP1", "AFF1", "AFF4",
        "CCAR2", "AP1G1", "AP2A1", "AP2A2", "AP2B1", "AP2M1"],
    "small-molecule and lipid kinases": [
        "ADK", "AK2", "CMPK1", "DCK", "PDXK", "PFKP", "NME2",
        "PIP4K2A", "PIP4K2A;PIP4K2B", "PIP4K2C"],
    "other enzymes": [
        "ACAD10", "ACAD11", "ACOX3", "CARS", "CERS2", "COMT", "DCTPP1", "DHCR24", "FECH",
        "GARS", "MAT2A", "NQO2", "PLA2G2A", "PYGB", "PYGL", "TOP2B", "PTPN1", "INPPL1",
        "ERCC2", "UBASH3B"],
    "other proteins": [
        "ACTR3", "AIMP1", "C2CD5", "CKAP5", "DDB1", "DNAJA1", "DYNLL1", "DYNLL2", "EIF3J",
        "EIF4G1", "EXOC4", "FANCI", "FIBP", "GAPVD1", "GRAMD1A", "GRB2", "HSP90AB2P",
        "IPO7", "KIAA0195", "KLHL6", "LASP1", "MYH10", "NPM1", "OSBPL3", "PAG1", "PDCD11",
        "PSMG1", "RPL10", "RPL35A", "RPL4", "RPS4X", "S100A8", "SLC16A3",
        "SLC25A5", "SLC5A1", "SMARCC1", "SMC2", "SNRNP200", "SNRPF", "SYNGR2", "TEX264",
        "U2AF1;U2AF1L4", "UNC119"],
    "not annotated as kinases by Klaeger et al. and flagged as possibly co-enriched with a kinase": [
        "BCR", "TRRAP"],
    "annotated as a kinase by Klaeger et al. but with no kinase domain (UniProt: encoded in the 3'-UTR of AAK1)": [
        "Q6ZSR9"],
}


def to_ranks(s):
    return rankdata(-np.asarray(s), method="average")


def neg_entropy(w):
    tot = w.sum(1, keepdims=True)
    p = w / np.where(tot <= 0, 1e-10, tot)
    return (p * np.where(p > 0, np.log2(np.where(p > 0, p, 1.0)), 0.0)).sum(1)


def s_score(M, t):
    return -(M > t).mean(1)


def entropy(M, b):
    return neg_entropy(np.maximum(M - b, 0.0))


def gini(M, b):
    s = np.sort(np.maximum(M - b, 0.0), 1)
    n, t = s.shape[1], s.sum(1)
    g = (2 * (np.arange(1, n + 1) * s).sum(1)) / (n * np.where(t == 0, 1, t)) - (n + 1) / n
    return np.where(t == 0, 0.0, g)


def ratio(M, k):
    s = -np.sort(-M, 1)
    return s[:, 0] - np.maximum(s[:, k], X0)


def candidate(M, x0=X0, b=X0):
    return neg_entropy(np.where(M > x0, T * np.logaddexp(0.0, (M - b) / T), 0.0))


SWEEPS = {
    "s_score": (s_score, np.arange(5.5, 8.25, 0.25)),
    "entropy": (entropy, np.arange(5.0, 6.75, 0.25)),
    "gini":    (gini,    np.arange(5.0, 6.75, 0.25)),
    "ratio":   (ratio,   range(1, 6)),
}
FAMS = list(SWEEPS)


def reversals(rows, keep):
    S = np.array(rows)[:, keep]
    iu = np.triu_indices(S.shape[1], 1)
    sign = np.sign(S[:, iu[0]] - S[:, iu[1]])
    return ((sign.max(0) > 0) & (sign.min(0) < 0)).mean()


def pstars(M, R=50, seed=42):
    meas = {"candidate": candidate, "entropy": lambda P: entropy(P, X0), "gini": lambda P: gini(P, X0),
            "s_score": lambda P: s_score(P, TAU), "ratio": lambda P: ratio(P, 1)}
    ref = {k: to_ranks(f(M)) for k, f in meas.items()}
    n = M.shape[1]
    panel = list(range(50, n, 30)) + [n]
    rng = np.random.RandomState(seed)
    res = {k: {ps: [] for ps in panel} for k in meas}
    for ps in panel:
        for _ in range(R):
            Ms = M[:, rng.choice(n, ps, replace=False)]
            for k, f in meas.items():
                res[k][ps].append(spearmanr(ref[k], to_ranks(f(Ms)))[0])
    return {k: next((ps for ps in panel if np.mean(res[k][ps]) > 0.90), None) for k in meas}


def summarize(M):
    out = {}
    n_act = (M > TAU).sum(1)
    act = n_act > 0
    out["kinases"] = M.shape[1]
    out["zero-active compounds"] = int((~act).sum())
    out["active share of matrix (%)"] = 100 * (M > TAU).mean()
    out["mean active kinases per compound"] = n_act.mean()

    ranks = {f: np.array([to_ranks(fn(M, p)) for p in grid]) for f, (fn, grid) in SWEEPS.items()}
    med = {f: np.median(r, 0) for f, r in ranks.items()}
    for i, a in enumerate(FAMS):
        for b in FAMS[i + 1:]:
            out[f"Table 1: {a} vs {b}"] = spearmanr(med[a], med[b])[0]
    # active compounds: re-ranked among themselves at every configuration (as in claim_checks.py, 2)
    med_act = {f: np.median([to_ranks(fn(M[act], p)) for p in grid], 0) for f, (fn, grid) in SWEEPS.items()}
    for i, a in enumerate(FAMS):
        for b in FAMS[i + 1:]:
            out[f"Table 1, active compounds: {a} vs {b}"] = spearmanr(med_act[a], med_act[b])[0]
    out["active compounds"] = int(act.sum())
    rng = np.random.default_rng(0)
    for lab, m in (("all", med), ("active", med_act)):
        diff = []
        for _ in range(1000):
            ix = rng.choice(len(m["gini"]), len(m["gini"]))
            diff.append(spearmanr(m["s_score"][ix], m["gini"][ix])[0] - spearmanr(m["gini"][ix], m["ratio"][ix])[0])
        out[f"r(S,G) - r(G,R), {lab}: 2.5%"] = np.percentile(diff, 2.5)
        out[f"r(S,G) - r(G,R), {lab}: 97.5%"] = np.percentile(diff, 97.5)
    sd = np.vstack(list(ranks.values())).std(0)
    out["pooled rank sd, zero-active"] = sd[~act].mean() if (~act).any() else np.nan
    out["pooled rank sd, active"] = sd[act].mean()
    out["Mann-Whitney p"] = mannwhitneyu(sd[~act], sd[act]).pvalue if (~act).any() else np.nan

    inst = {f: r.std(0) for f, r in ranks.items()}
    s = -np.sort(-M, 1)
    gap = s[:, 0] - s[:, 1]
    rng_ = np.array([np.ptp(r[r > TAU]) if (r > TAU).sum() > 1 else 0.0 for r in M])
    skw = np.array([skew(r[r > TAU]) if (r > TAU).sum() > 2 else 0.0 for r in M])
    for lab, feat, fam in [("n_active vs entropy instability", n_act, "entropy"),
                           ("top1-top2 gap vs ratio instability", gap, "ratio"),
                           ("active range vs S-score instability", rng_, "s_score"),
                           ("skewness vs Gini instability", skw, "gini")]:
        r, p = spearmanr(feat[act], inst[fam][act])
        out[f"Table 2: {lab}"] = r
        out[f"Table 2: {lab} (p)"] = p


    c = candidate(M)
    out["candidate vs n_active (orientation)"] = spearmanr(c, n_act)[0]
    out["candidate vs entropy, scored compounds"] = spearmanr(c[act], entropy(M, X0)[act])[0]
    added = candidate(np.hstack([M, np.full((M.shape[0], 1), TAU)]))
    out["D3 at tau*: % non-increasing (scored)"] = 100 * (added[act] <= c[act] + 1e-12).mean()
    for k, v in pstars(M).items():
        out[f"p* {k}"] = v
    return out


matrix = pd.read_csv("klaeger_matrix_released.csv", index_col=0)
cols = matrix.columns.tolist()
non_pk = [c for group in NON_PK.values() for c in group]
assert len(non_pk) == len(set(non_pk)) and set(non_pk) <= set(cols)
pk = pd.read_csv("klaeger_matrix.csv", index_col=0).columns.tolist()   # curated panel
M_all = matrix.values
M_pk = matrix[pk].values

print(f"Klaeger panel: {len(cols)} targets, {len(cols) - len(non_pk)} protein kinases ({len(pk)} curated, "
      f"without protein-group columns), {len(non_pk)} other targets")
for group, names in NON_PK.items():
    print(f"  {group} ({len(names)}): {', '.join(names)}")

act_all = M_all > TAU
act_np = matrix[non_pk].values > TAU
print(f"\nActive entries (pKd > {TAU:g}) on non-kinase columns: {act_np.sum()} of {act_all.sum()} "
      f"({100 * act_np.sum() / act_all.sum():.1f}%), in {int(act_np.any(1).sum())} of {len(matrix)} compounds")
only_np = act_np.any(1) & ~(M_pk > TAU).any(1)
print(f"Compounds whose only active targets are non-kinases: {int(only_np.sum())}"
      + (f" ({', '.join(matrix.index[only_np])})" if only_np.any() else ""))

full, prot = summarize(M_all), summarize(M_pk)
print(f"\n{'Result':48s} {'all 343':>10s} {'curated ' + str(len(pk)):>16s}")
for k in full:
    a, b = full[k], prot[k]
    fmt = (lambda v: f"{v:.3g}") if "(p)" in k or "Mann" in k else \
          (lambda v: f"{v}" if isinstance(v, (int, np.integer)) or v is None else f"{v:.3f}")
    print(f"{k:48s} {fmt(a):>10s} {fmt(b):>16s}")

# D2 example (top1-top2 gap < 0.05): how many near-ties involve a non-kinase column,
# and the R(1)/R(2) contrast on the protein-kinase panel
s_all = -np.sort(-M_all, 1)
tied = np.where((s_all[:, 0] - s_all[:, 1]) < 0.05)[0]
cols_arr = np.array(cols)
with_np = [i for i in tied if set(cols_arr[np.argsort(-M_all[i])[:2]]) & set(non_pk)]
print(f"\nCompounds with top1-top2 gap < 0.05: {len(tied)} on the full panel; "
      f"{len(with_np)} of them tie a non-kinase column ("
      + ", ".join(f"{matrix.index[i]}: {'/'.join(cols_arr[np.argsort(-M_all[i])[:2]])}" for i in with_np) + ")")
for lab, M in (("full panel", M_all), ("curated panel", M_pk)):
    s = -np.sort(-M, 1)
    t = (s[:, 0] - s[:, 1]) < 0.05
    print(f"  {lab}: {t.sum()} compounds with gap < 0.05; mean R(1) = {ratio(M, 1)[t].mean():.2f}, "
          f"mean R(2) = {ratio(M, 2)[t].mean():.2f}")
