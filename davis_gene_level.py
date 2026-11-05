"""
Davis: released constructs against curated reference targets.

The 433 released Davis constructs include mutant, phosphorylation-state and
single-domain constructs of the same gene (e.g. BRAF and BRAF(V600E); 12 EGFR
constructs), lipid kinases and three non-human kinases, so a compound's two
best-bound "targets" can be the wild-type and mutant forms of one kinase. The
analyses use one prespecified reference construct per human protein kinase gene
(curate_panels.py). This script recomputes the main Davis results on the released
constructs and, as an alternative to the reference construct, on the best
affinity per gene, printed beside the curated panel.
"""
import script_logging; script_logging.capture(__file__)
import os
os.environ['KMP_DUPLICATE_LIB_OK'] = 'TRUE'
import numpy as np
import pandas as pd
from scipy.stats import spearmanr, rankdata, skew

X0, TAU, T = 5.0, 6.0, 1.0


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
    first = {k: next((ps for ps in panel if np.mean(res[k][ps]) > 0.90), None) for k in meas}
    # reached only on the full panel: report "more than the largest subpanel tested"
    return {k: (f">{panel[-2]}" if v == n else v) for k, v in first.items()}


def summarize(M):
    out = {}
    n_act = (M > TAU).sum(1)
    act = n_act > 0
    out["targets"] = M.shape[1]
    out["mean active targets per compound"] = n_act.mean()
    ranks = {f: np.array([to_ranks(fn(M, p)) for p in grid]) for f, (fn, grid) in SWEEPS.items()}
    med = {f: np.median(r, 0) for f, r in ranks.items()}
    for i, a in enumerate(FAMS):
        for b in FAMS[i + 1:]:
            out[f"Table 1: {a} vs {b}"] = spearmanr(med[a], med[b])[0]
    inst = {f: r.std(0) for f, r in ranks.items()}
    s = -np.sort(-M, 1)
    feats = {"n_active": n_act, "top1-top2 gap": s[:, 0] - s[:, 1], "top2 pKd": s[:, 1],
             "active range": np.array([np.ptp(r[r > TAU]) if (r > TAU).sum() > 1 else 0.0 for r in M]),
             "active skew": np.array([skew(r[r > TAU]) if (r > TAU).sum() > 2 else 0.0 for r in M])}
    for feat, fam in [("n_active", "entropy"), ("n_active", "gini"), ("n_active", "ratio"),
                      ("top1-top2 gap", "ratio"), ("top2 pKd", "ratio")]:
        r, p = spearmanr(feats[feat][act], inst[fam][act])
        out[f"Table 2: {feat} vs {fam} instability"] = r
        out[f"Table 2: {feat} vs {fam} instability (p)"] = p
    out["compounds with top1-top2 gap < 0.05"] = int((feats["top1-top2 gap"] < 0.05).sum())
    out["mean overall rank sd"] = np.vstack(list(ranks.values())).std(0).mean()
    c = candidate(M)
    out["candidate vs entropy, scored compounds"] = spearmanr(c[act], entropy(M, X0)[act])[0]
    for k, v in pstars(M).items():
        out[f"p* {k}"] = v
    return out


prot = pd.read_csv("davis_proteins.csv").set_index("Protein_Index")["Gene_Name"]
drugs_cid = pd.read_csv("davis_drugs.csv")["CID"].values
rel = pd.read_csv("davis_affinity_released.csv").pivot(index="Drug_Index", columns="Protein_Index", values="Affinity").fillna(5.0)
cur = pd.read_csv("davis_affinity.csv").pivot(index="Drug_Index", columns="Protein_Index", values="Affinity").fillna(5.0)
names = prot.reindex(rel.columns)
gene = names.str.replace(r"\(.*$", "", regex=True).str.replace(r"-(phosphorylated|nonphosphorylated)$", "", regex=True)
cur_genes = set(gene[cur.columns])
hk = [c for c in rel.columns if gene[c] in cur_genes]          # human protein kinase constructs
Mg = rel[hk].T.groupby(gene[hk].values).max().T                 # best affinity per gene

print(f"Davis: {rel.shape[1]} released constructs ({gene.nunique()} distinct names); curated panel "
      f"{cur.shape[1]} reference targets; {Mg.shape[1]} genes for the best-affinity alternative")
same = [gene[t[0]] == gene[t[1]] for t in (rel.iloc[r].sort_values(ascending=False).index[:2] for r in range(len(rel)))]
print(f"Released panel: compounds whose two best-bound constructs are the same gene: {sum(same)} of {len(rel)}")
print(f"Released panel: EGFR constructs: {int((gene == 'EGFR').sum())}")
r_rank = np.array([to_ranks(ratio(rel.values, k)) for k in range(1, 6)])
worst = int(np.argmax(r_rank.std(0)))
print(f"Released panel: most ratio-unstable compound (rank sd over k = 1-5): row {worst}, "
      f"CID {drugs_cid[worst]}, sd {r_rank.std(0)[worst]:.1f}")

drugs = pd.read_csv("davis_drugs.csv")
gdc = int(drugs.index[drugs["CID"] == 11717001][0])   # GDC-0879
row = rel.iloc[gdc].sort_values(ascending=False)
crow = cur.iloc[gdc].sort_values(ascending=False)
print(f"GDC-0879 (row {gdc}): released top constructs "
      + ", ".join(f"{names[c]} {v:.3f}" for c, v in row.iloc[:3].items())
      + f" (gap {row.iloc[0] - row.iloc[1]:.3f}); curated top "
      + ", ".join(f"{prot[c]} {v:.3f}" for c, v in crow.iloc[:3].items())
      + f" (gap {crow.iloc[0] - crow.iloc[1]:.3f})")

res = {"released 433": summarize(rel.values), "curated": summarize(cur.values), "best per gene": summarize(Mg.values)}
labs = list(res)
print(f"\n{'Result':48s}" + "".join(f"{l:>15s}" for l in labs))
for k in res[labs[0]]:
    vals = [res[l][k] for l in labs]
    fmt = (lambda v: f"{v:.2g}") if "(p)" in k else \
          (lambda v: f"{v}" if isinstance(v, (int, np.integer, str)) or v is None else f"{v:.3f}")
    print(f"{k:48s}" + "".join(f"{fmt(v):>15s}" for v in vals))
