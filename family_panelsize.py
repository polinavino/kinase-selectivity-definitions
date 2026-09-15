"""
Does the Renyi order (family_construction.py's first axis) affect how quickly the
ranking converges to the full-panel ranking? candidate_measure.py's panel_convergence
established p*=110 (Klaeger) / 80 (Metz) for the alpha=1 (Shannon) candidate. This
script repeats that check across alpha and asks two follow-up questions before
treating any difference as a usable finding:

  1. Convergence table. p* by Renyi order, same panel grid and resampling scheme
     as candidate_measure.py's panel_convergence (R=50, seed=42), plus full
     convergence curves (mean of 3 seeds, R=100 each) to check the p* gap is not
     an artifact of one seed landing near a grid boundary.
  2. Tie diagnostic. A faster-converging measure that discriminates between fewer
     compounds would be uninteresting (less information, not more stability). This
     checks whether the number of distinct scores and the size of any tied group
     changes with alpha.
  3. Sparse/dense split. Whether agreement with the alpha=1 ranking holds up on
     compounds with few active kinases (n_active < 5), i.e. whether a
     faster-converging order buys its speed by moving the ranking most on
     precisely the least-stable compounds.
"""
import script_logging; script_logging.capture(__file__)
import os
import numpy as np
import pandas as pd
from scipy.stats import spearmanr, rankdata

HERE = os.path.dirname(os.path.abspath(__file__))


def to_ranks(s):
    return rankdata(-np.asarray(s), method='average')


def neg_renyi_H(w, alpha, eps=1e-12):
    """-H_alpha(p), matching family_construction.py's neg_renyi_H, extended with
    numerical limit proxies for alpha -> 0 (richness) and alpha -> inf
    (min-entropy, dominated by the single largest p_i)."""
    tot = w.sum(1, keepdims=True)
    tot = np.where(tot <= 0, eps, tot)
    p = w / tot
    if abs(alpha - 1.0) < 1e-9:
        return (p * np.where(p > 0, np.log2(p + eps), 0.0)).sum(1)
    if alpha > 30:
        return np.log2(np.maximum(p.max(1), eps))
    if alpha < 1e-3:
        return -np.log2(np.maximum((p > 0).sum(1), 1))
    s = np.where(p > 0, p ** alpha, 0.0).sum(1)
    s = np.maximum(s, eps)
    return (1.0 / (alpha - 1.0)) * np.log2(s)


def candidate_alpha(P, floor, T, alpha):
    g = (P > floor).astype(float)
    w = g * (T * np.logaddexp(0.0, (P - floor) / T))
    return neg_renyi_H(w, alpha)


def panel_convergence(M, floor, alpha, R=50, seed=42):
    n_drugs, n_kin = M.shape
    ref = to_ranks(candidate_alpha(M, floor, 1.0, alpha))
    panel = list(range(50, n_kin, 30)) + [n_kin]
    rng = np.random.RandomState(seed)
    res = {ps: [] for ps in panel}
    for ps in panel:
        for _ in range(R):
            idx = rng.choice(n_kin, ps, replace=False)
            res[ps].append(spearmanr(ref, to_ranks(candidate_alpha(M[:, idx], floor, 1.0, alpha)))[0])
    pstar = next((ps for ps in panel if np.mean(res[ps]) > 0.90), None)
    return panel, res, pstar


DATASETS = [
    dict(name='Klaeger', file='klaeger_matrix.csv', floor=5.0, active=6.0),
    dict(name='Metz',    file='metz_matrix.csv',    floor=4.0, active=6.0),
]
ALPHAS = [0.01, 0.25, 0.5, 0.75, 1.0, 2.0, 3.0]

print("=" * 78)
print("1. Panel-size convergence p* by Renyi order (grid and R as in candidate_measure.py)")
print("=" * 78)
curves = {}
for ds in DATASETS:
    M = pd.read_csv(f"{HERE}/{ds['file']}", index_col=0).values
    print(f"\n{ds['name']} ({M.shape[0]} cpd x {M.shape[1]} kin):")
    for alpha in ALPHAS:
        _, res, pstar = panel_convergence(M, ds['floor'], alpha)
        curves[(ds['name'], alpha)] = res
        print(f"  alpha={alpha:<5} p* = {pstar}")

print()
print("Full convergence curves (mean +/- sd Spearman rho vs full panel), 3 seeds x R=100,")
print("to confirm the p* gap above is not a single-seed artifact at a grid boundary:")
for ds in DATASETS:
    M = pd.read_csv(f"{HERE}/{ds['file']}", index_col=0).values
    print(f"\n{ds['name']}:")
    for alpha in ALPHAS:
        seed_curves = [panel_convergence(M, ds['floor'], alpha, R=100, seed=s)[1] for s in (1, 2, 3)]
        panel = sorted(seed_curves[0].keys())
        line = f"  alpha={alpha:<5}"
        for ps in panel:
            means = [np.mean(c[ps]) for c in seed_curves]
            line += f" n={ps}:{np.mean(means):.3f}"
        print(line)

print()
print("=" * 78)
print("2. Tie diagnostic: does a faster-converging order discriminate between fewer compounds?")
print("=" * 78)
for ds in DATASETS:
    M = pd.read_csv(f"{HERE}/{ds['file']}", index_col=0).values
    n_zero_active = int(((M > ds['active']).sum(1) == 0).sum())
    print(f"\n{ds['name']} ({M.shape[0]} compounds, {n_zero_active} with zero active kinases):")
    for alpha in ALPHAS + [50.0]:
        sel = candidate_alpha(M, ds['floor'], 1.0, alpha)
        _, counts = np.unique(np.round(sel, 8), return_counts=True)
        print(f"  alpha={alpha:<5} distinct scores = {len(counts):4d} / {M.shape[0]}   largest tie group = {counts.max()}")

print()
print("=" * 78)
print("3. Agreement with the alpha=1 (Shannon) ranking, split by active-set size")
print("=" * 78)
for ds in DATASETS:
    M = pd.read_csv(f"{HERE}/{ds['file']}", index_col=0).values
    nact = (M > ds['active']).sum(1)
    sparse = nact < 5
    shannon_rank = to_ranks(candidate_alpha(M, ds['floor'], 1.0, 1.0))
    print(f"\n{ds['name']}: {sparse.sum()} compounds with n_active<5 (sparse), "
          f"{(~sparse).sum()} with n_active>=5 (dense):")
    for alpha in [0.01, 0.5, 2.0, 50.0]:
        sel_rank = to_ranks(candidate_alpha(M, ds['floor'], 1.0, alpha))
        agree_sparse = spearmanr(sel_rank[sparse], shannon_rank[sparse])[0]
        agree_dense = spearmanr(sel_rank[~sparse], shannon_rank[~sparse])[0]
        print(f"  alpha={alpha:<5} sparse={agree_sparse:.3f}  dense={agree_dense:.3f}")

print("\nDone.")
