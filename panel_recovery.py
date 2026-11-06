"""
Recoverability from small panels: how much of each measure's full-panel ranking a
random subpanel preserves, for how many compounds (Figure 6).

All measures are compared on the same curated panels (curate_panels.py), the same
random subpanels (curated targets sampled uniformly without replacement, 50 per
size) and the same reporting gate: a compound is scored on a panel only if it has
a target above tau* there. Measures:
  S-score      hit fraction above tau*
  count        number of targets above tau* (negated); ranks every panel exactly as
               the S-score does, so its recovery is identical, although it passes
               weak-addition monotonicity and the S-score does not
  S-score x0   hit fraction above the floor x0 (any detected binding); passes
               weak-addition monotonicity, since every allowed addition is a hit
  entropy x0   baseline-subtracted entropy, beta = x0
  entropy tau* baseline-subtracted entropy, beta = tau*
  Gini         baseline-subtracted Gini, beta = x0
  ratio        R(1), clipping level x0
  candidate    gated smooth-hinge entropy (Methods)
  entropy Ka   published association-constant entropy, p_i ~ 1/Kd_i, with the
               floor-coded entries given zero weight (affinity datasets only)
  partition Ka maximum partition index (Cheng et al. 2010), max_i Ka_i / sum_i Ka_i,
               floor-coded entries given zero weight (affinity datasets only); it
               satisfies weak-addition monotonicity, since an addition at or below
               tau* leaves the maximum unchanged and raises the denominator

For each panel size m:
  coverage       share of the compounds scored on the full panel that are also
                 scored on the subpanel (the same for every measure, by the gate);
  Spearman       between full-panel and subpanel scores, over compounds scored on both;
  pair recovery  share of full-panel pairs separated by at least a quarter of the
                 scored compounds' ranks that are scored on the subpanel and keep
                 their order (ties count one half): coverage and ordering together.

Influence diagnostic: for each compound, the largest change in its score when one
target is left out of the full panel, divided by the spread (SD) of that measure's
scores across compounds; and, for the entropy-like measures, the effective number
of weighted targets (sum w)^2 / sum w^2. Each is correlated with the compound's
recovery error at m = 50 (mean absolute rank change, as a share of the compounds,
over the subpanels on which it is scored).
"""
import script_logging; script_logging.capture(__file__)
import os
os.environ['KMP_DUPLICATE_LIB_OK'] = 'TRUE'
import numpy as np
import pandas as pd
from scipy.stats import spearmanr, rankdata
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
R_PANELS = 50
M_DIAG = 50


def neg_entropy(w):
    tot = w.sum(-1, keepdims=True)
    p = w / np.where(tot <= 0, 1.0, tot)
    return np.where(p > 0, p * np.log2(np.where(p > 0, p, 1.0)), 0.0).sum(-1)


def gini(M, b):
    s = np.sort(np.maximum(M - b, 0.0), -1)
    n = s.shape[-1]
    t = s.sum(-1)
    g = (2 * (np.arange(1, n + 1) * s).sum(-1)) / (n * np.where(t == 0, 1, t)) - (n + 1) / n
    return np.where(t == 0, 0.0, g)


def measures(x0, tau, affinity=True):
    T = tau - x0
    m = {
        "S-score":      lambda M: -(M > tau).mean(-1),
        "count":        lambda M: -(M > tau).sum(-1).astype(float),
        "S-score x0":   lambda M: -(M > x0).mean(-1),
        "entropy x0":   lambda M: neg_entropy(np.maximum(M - x0, 0.0)),
        "entropy tau*": lambda M: neg_entropy(np.maximum(M - tau, 0.0)),
        "Gini":         lambda M: gini(M, x0),
        "ratio":        lambda M: (lambda s: s[..., 0] - np.maximum(s[..., 1], x0))(-np.sort(-M, -1)),
        "candidate":    lambda M: neg_entropy(np.where(M > x0, T * np.logaddexp(0.0, (M - x0) / T), 0.0)),
    }
    if affinity:
        m["entropy Ka"] = lambda M: neg_entropy(np.where(M > x0, 10.0 ** (M - 9.0), 0.0))
        m["partition Ka"] = lambda M: (lambda K: K.max(-1) / np.where(K.sum(-1) > 0, K.sum(-1), 1.0))(
            np.where(M > x0, 10.0 ** (M - 9.0), 0.0))
        # sensitivity (printed, not plotted): floor-coded entries kept at the floor value
        m["entropy Ka, floor kept"] = lambda M: neg_entropy(10.0 ** (M - 9.0))
    return m


def weights(M, x0, tau, name):
    T = tau - x0
    if name == "entropy x0":
        return np.maximum(M - x0, 0.0)
    if name == "entropy tau*":
        return np.maximum(M - tau, 0.0)
    if name == "entropy Ka":
        return np.where(M > x0, 10.0 ** (M - 9.0), 0.0)
    return np.where(M > x0, T * np.logaddexp(0.0, (M - x0) / T), 0.0)


def pair_recovery(full, sub, ok_sub, elig):
    """Share of well-separated full-panel pairs that are scored and correctly ordered."""
    idx = np.where(elig)[0]
    r = rankdata(-full[idx])
    sep = np.abs(r[:, None] - r[None, :]) >= len(idx) / 4
    iu = np.triu(sep, 1)
    a, b = np.where(iu)
    ia, ib = idx[a], idx[b]
    both = ok_sub[ia] & ok_sub[ib]
    sf, ss = np.sign(full[ia] - full[ib]), np.sign(sub[ia] - sub[ib])
    credit = np.where(ss == sf, 1.0, np.where(ss == 0, 0.5, 0.0)) * both
    return credit.mean()


aff = pd.read_csv(f"{HERE}/davis_affinity.csv")
DATA = {
    "Davis":         (aff.pivot(index="Drug_Index", columns="Protein_Index", values="Affinity").fillna(5.0).values, 5.0, 6.0),
    "Klaeger":       (pd.read_csv(f"{HERE}/klaeger_matrix.csv", index_col=0).values, 5.0, 6.0),
    "Anastassiadis": (pd.read_csv(f"{HERE}/anastassiadis_matrix.csv", index_col=0).values, 20.0, 50.0),
    "Metz":          (pd.read_csv(f"{HERE}/metz_matrix.csv", index_col=0).values, 4.0, 6.0),
}
SIZES = [10, 20, 30, 40, 50, 65, 80, 100, 125, 150, 200, 250, 300]

curves = {}
for ds, (M, x0, tau) in DATA.items():
    n, N = M.shape
    meas = measures(x0, tau, affinity=(ds != "Anastassiadis"))
    elig = (M > tau).any(1)
    full = {k: f(M) for k, f in meas.items()}
    sizes = [m for m in SIZES if m < N] + [N]
    rng = np.random.RandomState(42)
    res = {k: {"rho": {}, "pair": {}} for k in meas}
    cover = {}
    err = {k: np.zeros(n) for k in meas}
    err_n = np.zeros(n)
    for m in sizes:
        cov, rho, pair = [], {k: [] for k in meas}, {k: [] for k in meas}
        for _ in range(R_PANELS if m < N else 1):
            Ms = M[:, rng.choice(N, m, replace=False)]
            ok = (Ms > tau).any(1)
            both = ok & elig
            cov.append(both.sum() / elig.sum())
            for k, f in meas.items():
                s = f(Ms)
                rho[k].append(spearmanr(full[k][both], s[both])[0] if both.sum() > 2 else np.nan)
                pair[k].append(pair_recovery(full[k], s, ok, elig))
                if m == M_DIAG:
                    rf = rankdata(-full[k][both]) / both.sum()
                    rs = rankdata(-s[both]) / both.sum()
                    err[k][both] += np.abs(rf - rs)
            if m == M_DIAG:
                err_n += both
        cover[m] = np.mean(cov)
        for k in meas:
            res[k]["rho"][m] = (np.nanmean(rho[k]), np.nanstd(rho[k]), np.mean(np.array(rho[k]) > 0.90))
            res[k]["pair"][m] = (np.mean(pair[k]), np.std(pair[k]))
    curves[ds] = (sizes, cover, res)

    print(f"\n=== {ds}: {n} compounds x {N} curated targets; {elig.sum()} scored on the full panel; {R_PANELS} random subpanels per size ===")
    print("coverage (share of full-panel-scored compounds still scored): "
          + ", ".join(f"{m}: {100 * cover[m]:.0f}%" for m in sizes))
    print(f"{'measure':13s} " + "".join(f"{m:>7d}" for m in sizes) + "   (mean Spearman over compounds scored on both)")
    for k in meas:
        print(f"{k:13s} " + "".join(f"{res[k]['rho'][m][0]:7.3f}" for m in sizes))
    print(f"{'measure':13s} " + "".join(f"{m:>7d}" for m in sizes) + "   (pair recovery)")
    for k in meas:
        print(f"{k:13s} " + "".join(f"{res[k]['pair'][m][0]:7.3f}" for m in sizes))
    for k in meas:
        ps = next((m for m in sizes if res[k]["rho"][m][0] > 0.90), None)
        pp = next((m for m in sizes if res[k]["pair"][m][0] > 0.90), None)
        prev = sizes[sizes.index(ps) - 1] if ps and sizes.index(ps) > 0 else None
        mu, sd, share = res[k]["rho"][ps]
        print(f"  {k:13s} p* (Spearman > 0.90): {ps} (largest smaller size tested: {prev}; at p*: mean {mu:.3f}, "
              f"SD {sd:.3f}, {100 * share:.0f}% of subpanels above 0.90);  smallest size with pair recovery > 0.90: {pp}")
    diff = max(abs(res['S-score']['rho'][m][0] - res['count']['rho'][m][0]) for m in sizes)
    print(f"  S-score vs count: largest difference between their Spearman curves = {diff:.2e}")
    print("  full-panel Spearman with the S-score at x0, over scored compounds: "
          + ", ".join(f"{k} {spearmanr(full[k][elig], full['S-score x0'][elig])[0]:.4f}"
                      for k in meas if k != "S-score x0"))

    # p* at other bars (the 0.90 bar is a chosen descriptive benchmark)
    for bar in (0.80, 0.85, 0.95):
        print(f"  p* at bar {bar:.2f}: " + ", ".join(
            f"{k} {next((m for m in sizes if res[k]['rho'][m][0] > bar), None)}" for k in meas))
    # who drops out under D1: probability that a compound with a active targets keeps one
    # on a uniform random 50-target subpanel (hypergeometric)
    from math import comb
    print("  probability of staying scored at 50 targets, by number of active targets: "
          + ", ".join(f"{a}: {1 - comb(N - a, 50) / comb(N, 50):.2f}" for a in (1, 2, 5, 10)))
    # how close the share-following scores are to the gated entropy
    for k in ("entropy x0", "Gini", "S-score x0"):
        dr = max(abs(res[k]["rho"][m][0] - res["candidate"]["rho"][m][0]) for m in sizes)
        dp = max(abs(res[k]["pair"][m][0] - res["candidate"]["pair"][m][0]) for m in sizes)
        print(f"  {k} vs gated entropy: largest difference in mean Spearman {dr:.3f}, in recoverability {dp:.3f}")

    # influence diagnostic on the full panel
    if M_DIAG in sizes:
        e = {k: np.where(err_n > 0, err[k] / np.maximum(err_n, 1), np.nan) for k in meas}
        print(f"  Influence vs recovery error at m = {M_DIAG} (Spearman over scored compounds):")
        for k, f in meas.items():
            loo = np.zeros(n)
            for j in range(N):
                keep = np.ones(N, bool); keep[j] = False
                loo = np.maximum(loo, np.abs(f(M[:, keep]) - full[k]))
            infl = loo / np.std(full[k][elig])
            ok = elig & ~np.isnan(e[k])
            if np.ptp(infl[ok]) == 0:
                print(f"    {k:13s} max leave-one-out influence: the same for every scored compound (one target)")
                continue
            line = f"    {k:13s} max leave-one-out influence: r = {spearmanr(infl[ok], e[k][ok])[0]:+.2f}"
            if k in ("entropy x0", "entropy tau*", "candidate", "entropy Ka"):
                w = weights(M, x0, tau, k)
                neff = w.sum(1) ** 2 / np.maximum((w ** 2).sum(1), 1e-12)
                line += f";  effective number of targets: r = {spearmanr(neff[ok], e[k][ok])[0]:+.2f}"
            print(line)

# ------------------------------------------------------------------ Figure 6
ORDER = ["candidate", "entropy x0", "entropy tau*", "S-score", "Gini", "ratio", "S-score x0", "entropy Ka",
         "partition Ka"]
# plotted: the entropy and Gini at x0 and the S-score at x0 lie on the gated entropy's curve
# (differences printed above), so they are not drawn
PLOT = ["candidate", "S-score", "entropy tau*", "ratio", "entropy Ka", "partition Ka"]
plt.rcParams.update({"font.size": 13, "axes.titlesize": 13, "axes.labelsize": 13,
                     "xtick.labelsize": 11, "ytick.labelsize": 11})
LABEL = {"candidate": "gated entropy (≈ entropy, Gini, S-score at x₀)", "entropy x0": "entropy, β = x₀", "entropy tau*": "entropy, β = τ*",
         "S-score": "S-score, τ = τ* (= active count)", "Gini": "Gini, β = x₀", "ratio": "ratio R(1)",
         "S-score x0": "S-score, τ = x₀", "entropy Ka": "entropy, published (1/Kd)",
         "partition Ka": "partition index, max (1/Kd)"}
COLOR = dict(zip(ORDER, ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948",
                         "#e34948"]))   # the two association-constant measures share a hue, told apart by line style
MARK = dict(zip(ORDER, ["o", "s", "^", "D", "v", "P", "X", "*", "<"]))
STYLE = {k: ("--" if k == "partition Ka" else "-") for k in ORDER}
fig, axes = plt.subplots(3, 4, figsize=(11, 8.5), sharex="col",
                         gridspec_kw={"height_ratios": [1.2, 3, 3]})
for c, (ds, (sizes, cover, res)) in enumerate(curves.items()):
    axc = axes[0, c]
    axc.plot(sizes, [cover[m] for m in sizes], color="#555555", lw=2, marker="o", ms=4)
    axc.set_ylim(0, 1.01)
    axc.grid(alpha=0.3, lw=0.5)
    axc.set_title(f"{ds} ({sizes[-1]})")
    if c == 0:
        axc.set_ylabel("coverage")
    for r, (key, ylab) in enumerate([("pair", "recoverability"),
                                     ("rho", "Spearman with\nfull-panel ranking")], start=1):
        ax = axes[r, c]
        for k in [k for k in PLOT if k in res]:
            mu = np.array([res[k][key][m][0] for m in sizes])
            ax.plot(sizes, mu, color=COLOR[k], lw=2, ls=STYLE[k], marker=MARK[k], ms=4, label=LABEL[k])
            if key == "rho":   # +/- 1 SD across subpanels
                sd = np.array([res[k][key][m][1] for m in sizes])
                ax.fill_between(sizes, mu - sd, mu + sd, color=COLOR[k], alpha=0.10, lw=0)
        ax.axhline(0.9, color="#888888", lw=1, ls="--")
        ax.set_ylim(0, 1.01)
        ax.grid(alpha=0.3, lw=0.5)
        if c == 0:
            ax.set_ylabel(ylab)
        if r == 2:
            ax.set_xlabel("targets")
h, l = axes[1, 0].get_legend_handles_labels()
fig.legend(h, l, loc="lower center", ncol=3, fontsize=12, frameon=False, bbox_to_anchor=(0.5, -0.09))
plt.tight_layout()
plt.savefig("panel_recovery.png", dpi=150, bbox_inches="tight")
print("\nSaved panel_recovery.png")
