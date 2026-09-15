"""
Reproduces the Bosc et al. / CATDS comparison referenced in the Discussion
(response to Reviewer 1): are Bosc et al.'s window score and ranking score
distinct from the ratio family R(k), or the same underlying construct; and
does a CATDS-style proxy correlate more with the distribution family
(entropy, Gini) than R(k) does.

Window score (WS) and ranking score (RS) are implemented exactly as defined
in Bosc, Meyer & Bonnet 2017 (BMC Bioinformatics 18:17):
  WS(x) = (# affinities within `window` log-units [or %] of the top hit) / n_kinases
          -- lower is more selective in the original paper; negated here so
          that, like every other score in this repository, higher = more selective.
  RS(x) = (top affinity) - (affinity at rank position k)
          -- higher is more selective, as in the original paper.
Parameter grids follow the original paper: windows of 2/1/0.5 log-units for
the pKd/pKi datasets (Davis, Klaeger, Metz) and 20/10/5 percentage points for
the percent-inhibition dataset (Anastassiadis); rank positions 20/10/5 for all
four datasets.

CATDS (Klaeger et al. 2017) divides each target's concentration-dependent
binding reduction by the summed reduction across the whole panel, which
requires raw dose-response curves this repository does not have (only the
derived apparent-Kd matrices). The proxy used here applies the same
normalization idea to the baseline-subtracted affinity profile already used
for entropy/Gini: catds_proxy(x) = max_i(shifted_i) / sum_j(shifted_j), i.e.
the top target's share of the total panel-wide binding-reduction budget.
Higher share = more concentrated = more selective.
"""
import script_logging; script_logging.capture(__file__)
import os
import numpy as np
import pandas as pd
from scipy.stats import spearmanr, rankdata

HERE = os.path.dirname(os.path.abspath(__file__))


def to_ranks(scores):
    return rankdata(-np.asarray(scores), method='average')


def ratio_selectivity(M, top_n, floor):
    out = []
    for row in M:
        s = np.sort(row)[::-1]
        off = s[top_n] if len(s) > top_n else floor
        out.append(s[0] - max(off, floor))
    return np.array(out)


def window_score(M, window, floor):
    top = M.max(axis=1, keepdims=True)
    within = (M >= np.maximum(top - window, floor)).sum(axis=1)
    return -(within / M.shape[1])  # negated: higher = more selective


def ranking_score(M, k, floor):
    out = []
    for row in M:
        s = np.sort(row)[::-1]
        at_k = s[k - 1] if len(s) >= k else floor
        out.append(s[0] - max(at_k, floor))
    return np.array(out)


def catds_proxy(M, baseline):
    shifted = np.maximum(M - baseline, 0)
    rs = shifted.sum(axis=1)
    top = shifted.max(axis=1)
    return np.where(rs == 0, 0.0, top / np.where(rs == 0, 1.0, rs))


def selectivity_entropy(M, baseline, eps=1e-10):
    sh = np.maximum(M - baseline, 0)
    rs = sh.sum(axis=1, keepdims=True)
    rs = np.where(rs == 0, eps, rs)
    p = sh / rs
    return -(-(p * np.where(p > 0, np.log2(p + eps), 0)).sum(axis=1))


def gini_selectivity(M, baseline):
    sh = np.maximum(M - baseline, 0)
    out = []
    for row in sh:
        rs = np.sort(row); n = len(rs); tot = rs.sum()
        out.append(0.0 if tot == 0
                   else (2 * np.sum(np.arange(1, n + 1) * rs)) / (n * tot) - (n + 1) / n)
    return np.array(out)


# per-dataset: (matrix, floor, window sweep, catds/entropy/gini baseline sweep, label)
def load():
    davis = pd.read_csv(f"{HERE}/davis_affinity.csv").pivot(
        index="Drug_Index", columns="Protein_Index", values="Affinity").fillna(5.0).values
    return {
        "Davis":         (davis, 5.0, [2.0, 1.0, 0.5], np.arange(5.0, 6.75, .25)),
        "Klaeger":       (pd.read_csv(f"{HERE}/klaeger_matrix.csv", index_col=0).values, 5.0, [2.0, 1.0, 0.5], np.arange(5.0, 6.75, .25)),
        "Anastassiadis": (pd.read_csv(f"{HERE}/anastassiadis_matrix.csv", index_col=0).values, 0.0, [20.0, 10.0, 5.0], np.arange(0, 35, 5)),
        "Metz":          (pd.read_csv(f"{HERE}/metz_matrix.csv", index_col=0).values, 4.0, [2.0, 1.0, 0.5], np.arange(4.0, 5.75, .25)),
    }


def analyze(name, M, floor, windows, baselines):
    med_ratio = np.median([to_ranks(ratio_selectivity(M, k, floor)) for k in range(1, 6)], axis=0)
    med_ws = np.median([to_ranks(window_score(M, w, floor)) for w in windows], axis=0)
    med_rs = np.median([to_ranks(ranking_score(M, k, floor)) for k in (20, 10, 5)], axis=0)
    med_catds = np.median([to_ranks(catds_proxy(M, b)) for b in baselines], axis=0)
    med_entropy = np.median([to_ranks(selectivity_entropy(M, b)) for b in baselines], axis=0)
    med_gini = np.median([to_ranks(gini_selectivity(M, b)) for b in baselines], axis=0)
    return dict(
        dataset=name, n_cpd=M.shape[0], n_kin=M.shape[1],
        ws_vs_ratio=round(spearmanr(med_ws, med_ratio)[0], 3),
        rs_vs_ratio=round(spearmanr(med_rs, med_ratio)[0], 3),
        catds_vs_ratio=round(spearmanr(med_catds, med_ratio)[0], 3),
        catds_vs_entropy=round(spearmanr(med_catds, med_entropy)[0], 3),
        catds_vs_gini=round(spearmanr(med_catds, med_gini)[0], 3),
    )


if __name__ == "__main__":
    rows = [analyze(name, M, floor, windows, baselines)
            for name, (M, floor, windows, baselines) in load().items()]
    df = pd.DataFrame(rows)
    df.to_csv(f"{HERE}/bosc_catds_comparison.csv", index=False)
    print(df.to_string(index=False))
    ws_lo, ws_hi = df['ws_vs_ratio'].min(), df['ws_vs_ratio'].max()
    rs_lo, rs_hi = df['rs_vs_ratio'].min(), df['rs_vs_ratio'].max()
    print(f"\nWindow score vs ratio family, range across datasets: {ws_lo:.2f}-{ws_hi:.2f}")
    print(f"Ranking score vs ratio family, range across datasets: {rs_lo:.2f}-{rs_hi:.2f}")
    kl = df[df['dataset'] == 'Klaeger'].iloc[0]
    print(f"\nOn Klaeger: CATDS proxy vs ratio = {kl['catds_vs_ratio']}, "
          f"vs entropy = {kl['catds_vs_entropy']}, vs Gini = {kl['catds_vs_gini']}")
