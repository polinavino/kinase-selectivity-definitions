"""
Monotonicity under weak off-target addition (the paper's D3) for the entropy and
the candidate.

As stated in Methods: adding a kinase with affinity in (x0, tau*] must not raise
a compound's selectivity score. This script supports Methods (Propositions 1
and 3, and D4) and Results, "Criteria and the Gated Entropy" (SI Section S4):

1. The condition on the smoothing width T under which the candidate provably
   satisfies D3 (with emphasis baseline beta = floor x0): every gated weight is
   above m = T ln 2, an added kinase at or below tau* has weight at most
   a = T ln(1 + exp((tau* - x0)/T)), and the proof needs a/m <= e, i.e.
   T >= (tau* - x0) / ln(2^e - 1) ~= 0.58 (tau* - x0). With T = tau* - x0 the
   same argument covers every beta in B = [x0, tau*], with m = T ln(1 + exp((x0 - beta)/T)).
2. That the condition is sharp: just below it, a profile with one strong target
   and enough kinases just above the floor violates D3.
3. A worked counterexample for the baseline-subtracted entropy.
4. How often, on the real data, adding one kinase at tau* raises the entropy's
   selectivity score, and that it never raises the candidate's.
5. A random adversarial search for candidate violations with T = tau* - x0.
6. How much moving T to tau* - x0 changes the Metz and Anastassiadis rankings.
7. Two limits of the candidate: its jump at the floor, and that strengthening an
   existing off-target (rather than adding one) can raise its score.
8. A D3 counterexample for the published entropy.
9. D4 (orientation) on reference profiles with m equally strong targets.
10. Proposition 3 (turning a null pair into a weak target never raises -H, P_max
    or the Gini), and the examples that show why its condition is needed.
11. D4 for the published entropy and P_max with floor-coded pairs kept at the floor.
12. D4 for the entropy and Gini over the Anastassiadis baseline sweep, which
    includes baselines below the floor.
"""
import script_logging; script_logging.capture(__file__)
import os
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
T_FACTOR = 1.0 / np.log(2.0 ** np.e - 1.0)   # T >= T_FACTOR * (tau* - x0)


def neg_entropy(w):
    """Selectivity orientation (-H, bits) of each row of nonnegative weights."""
    w = np.atleast_2d(np.asarray(w, float))
    tot = w.sum(1, keepdims=True)
    p = w / np.where(tot <= 0, 1.0, tot)
    return np.where(p > 0, p * np.log2(np.where(p > 0, p, 1.0)), 0.0).sum(1)


def candidate(P, x0, T):
    P = np.atleast_2d(np.asarray(P, float))
    return neg_entropy(np.where(P > x0, T * np.logaddexp(0.0, (P - x0) / T), 0.0))


def entropy(P, beta):
    return neg_entropy(np.maximum(np.atleast_2d(np.asarray(P, float)) - beta, 0.0))


def add_column(P, value):
    P = np.atleast_2d(np.asarray(P, float))
    return np.hstack([P, np.full((P.shape[0], 1), value)])


def worst_rise_near_floor(x0, tau, T, top, n_max=400):
    """Largest rise in the candidate's score from adding a kinase at tau*, over
    profiles with one strong target and n kinases just above the floor."""
    worst, worst_n = 0.0, None
    for n in range(1, n_max + 1):
        prof = [top] + [x0 + 1e-3 * T] * n
        d = (candidate(add_column(prof, tau), x0, T) - candidate(prof, x0, T))[0]
        if d > worst:
            worst, worst_n = d, n
    return worst, worst_n


def load():
    aff = pd.read_csv(f"{HERE}/davis_affinity.csv")
    davis = aff.pivot(index="Drug_Index", columns="Protein_Index", values="Affinity").fillna(5.0).values
    # (matrix, floor x0, D1 threshold tau*, top value for constructed profiles,
    #  entropy baseline grid used in the main analysis, T used before this change)
    return {
        "Davis":         (davis, 5.0, 6.0, 8.0, np.arange(5.0, 6.75, 0.25), 1.0),
        "Klaeger":       (pd.read_csv(f"{HERE}/klaeger_matrix.csv", index_col=0).values,
                          5.0, 6.0, 8.0, np.arange(5.0, 6.75, 0.25), 1.0),
        "Anastassiadis": (pd.read_csv(f"{HERE}/anastassiadis_matrix.csv", index_col=0).values,
                          20.0, 50.0, 90.0, np.arange(0.0, 35.0, 5.0), 10.0),
        "Metz":          (pd.read_csv(f"{HERE}/metz_matrix.csv", index_col=0).values,
                          4.0, 6.0, 8.0, np.arange(4.0, 5.75, 0.25), 1.0),
    }


data = load()

print("1. Condition for the candidate's D3 guarantee: T >= %.4f (tau* - x0)" % T_FACTOR)
print(f"   {'dataset':14s} {'tau*-x0':>8s} {'T now':>6s} {'T before':>9s} {'r=a/m now':>10s} {'r before':>9s}  (need r <= e = {np.e:.3f})")
for name, (M, x0, tau, top, grid, T_old) in data.items():
    d = tau - x0
    r = lambda T: np.logaddexp(0.0, d / T) / np.log(2.0)
    print(f"   {name:14s} {d:8.1f} {d:6.1f} {T_old:9.1f} {r(d):10.3f} {r(T_old):9.3f}")
u = np.linspace(0.0, 1.0, 1001)   # (beta - x0) / (tau* - x0), with T = tau* - x0
r_B = np.logaddexp(0.0, 1.0 - u) / np.logaddexp(0.0, -u)
print(f"   emphasis baseline beta swept over B, T = tau* - x0 (any dataset): "
      f"r = a/m from {r_B.min():.3f} (beta = x0) to {r_B.max():.3f} (beta = tau*)")

print("\n2. Sharpness on profiles with one strong target and n kinases just above the floor")
for name in ("Klaeger", "Metz", "Anastassiadis"):
    _, x0, tau, top, _, T_old = data[name]
    d = tau - x0
    for label, T in (("T = tau*-x0", d), ("T before", T_old),
                     ("T just above bound", 1.01 * T_FACTOR * d), ("T just below bound", 0.97 * T_FACTOR * d)):
        rise, n = worst_rise_near_floor(x0, tau, T, top, n_max=2000)
        verdict = "D3 holds" if rise <= 1e-12 else f"D3 VIOLATED (max rise {rise:+.1e} bits at n = {n})"
        print(f"   {name:14s} {label:20s} T = {T:7.3f}: {verdict}")

print("\n3. Worked counterexample for the entropy (beta = floor = 5.0, tau* = 6.0)")
prof = [7.0] + [5.1] * 50
before, after = entropy(prof, 5.0)[0], entropy(add_column(prof, 6.0), 5.0)[0]
print(f"   profile: one kinase at pKd 7.0 and fifty at 5.1; add one kinase at pKd 6.0")
print(f"   entropy selectivity  {before:+.4f} -> {after:+.4f}  ({'RISES: D3 violated' if after > before else 'does not rise'})")
cb, ca = candidate(prof, 5.0, 1.0)[0], candidate(add_column(prof, 6.0), 5.0, 1.0)[0]
print(f"   candidate (T = 1)    {cb:+.4f} -> {ca:+.4f}  ({'RISES' if ca > cb else 'does not rise'})")
e6, e51 = 6.0 - 5.0, 5.1 - 5.0
g6, g51 = float(np.logaddexp(0.0, 6.0 - 5.0)), float(np.logaddexp(0.0, 5.1 - 5.0))
print(f"   weights: entropy {e6:.2f} (added kinase at 6.0) and {e51:.2f} (each at 5.1), ratio {e6 / e51:.1f}; "
      f"candidate {g6:.4f} and {g51:.4f}, ratio {g6 / g51:.2f}")

print("\n4. Real data: add one kinase at tau* to every compound with an active kinase")
for name, (M, x0, tau, top, grid, T_old) in data.items():
    act = (M > tau).any(1)
    Ma = M[act]
    ent_floor = (entropy(add_column(Ma, tau), x0) > entropy(Ma, x0) + 1e-12)
    ent_any = np.zeros(len(Ma), bool)
    for b in grid:
        ent_any |= (entropy(add_column(Ma, tau), b) > entropy(Ma, b) + 1e-12)
    d = tau - x0
    cand_now = (candidate(add_column(Ma, tau), x0, d) > candidate(Ma, x0, d) + 1e-12)
    cand_old = (candidate(add_column(Ma, tau), x0, T_old) > candidate(Ma, x0, T_old) + 1e-12)
    print(f"   {name:14s} n={act.sum():4d}  entropy rises: {100 * ent_floor.mean():5.1f}% at beta = floor, "
          f"{100 * ent_any.mean():5.1f}% at some beta in the sweep;  candidate rises: "
          f"{100 * cand_now.mean():.1f}% (T = tau*-x0), {100 * cand_old.mean():.1f}% (T before)")
    by_beta = ", ".join(f"{b:g}: {100 * (entropy(add_column(Ma, tau), b) > entropy(Ma, b) + 1e-12).mean():.1f}%"
                        for b in grid)
    print(f"   {'':14s} entropy rises, by baseline: {by_beta}")

print("\n5. Random adversarial search, candidate with T = tau* - x0 (additions anywhere in (x0, tau*])")
rng = np.random.default_rng(0)
for name in ("Klaeger", "Metz", "Anastassiadis"):
    _, x0, tau, top, _, _ = data[name]
    d = tau - x0
    worst = 0.0
    for _ in range(20000):
        n = rng.integers(1, 300)
        prof = np.concatenate([[rng.uniform(tau + 1e-3, top + d)],
                               x0 + d * rng.uniform(1e-4, 1.0, n - 1) ** rng.uniform(1, 8)])
        add = rng.uniform(x0 + 1e-6, tau)
        worst = max(worst, (candidate(add_column(prof, add), x0, d) - candidate(prof, x0, d))[0])
    print(f"   {name:14s} 20000 random profiles: max rise = {worst:+.2e} bits")

print("\n6. Effect of the change of T on the rankings (Spearman, all compounds)")
from scipy.stats import spearmanr
for name, (M, x0, tau, top, grid, T_old) in data.items():
    d = tau - x0
    if T_old == d:
        print(f"   {name:14s} T unchanged ({d:g})")
        continue
    r = spearmanr(candidate(M, x0, T_old), candidate(M, x0, d))[0]
    print(f"   {name:14s} T = {T_old:g} -> {d:g}: Spearman = {r:.3f}")


print("\n7. Limits of the candidate (beta = x0 = 5.0, T = 1)")
one = lambda xs: float(candidate(np.array([xs]), 5.0, 1.0)[0])
print(f"   jump at the floor: one target at 7 with a second at 5.0: {one([7.0, 5.0]):+.4f} bits; "
      f"second at 5.0 + 1e-8: {one([7.0, 5.0 + 1e-8]):+.4f} bits")
base = [7.0] + [5.1] * 50
stronger = [7.0, 6.0] + [5.1] * 49
print(f"   one target at 7 and fifty at 5.1: {one(base):+.5f}; one of the fifty raised to 6.0: {one(stronger):+.5f} "
      f"({'higher' if one(stronger) > one(base) else 'not higher'}: strengthening an existing off-target "
      f"is not covered by D3)")

print("\n8. The published (association-constant) entropy, floor-coded pairs given zero weight")
def ka_score(xs, floor=5.0):
    x = np.asarray(xs, float)
    w = np.where(x > floor, 10.0 ** (x - 9.0), 0.0)
    p = w / w.sum()
    return float((p * np.log2(np.where(p > 0, p, 1.0))).sum())
base = [6.01] + [5.01] * 50
print(f"   one target at 6.01 and fifty at 5.01: {ka_score(base):+.4f} bits; one kinase added at 6.0: "
      f"{ka_score(base + [6.0]):+.4f} bits ({'rises: D3 violated' if ka_score(base + [6.0]) > ka_score(base) else 'does not rise'})")

print("\n9. D4 (orientation) on reference profiles: m targets at a common value a, the rest at the floor,")
print("   panel size N fixed; D4 needs f(v_{m+1}) <= f(v_m) for all m and f(v_1) > f(v_N)")
X0, TAU, NREF = 5.0, 6.0, 30


def gini_w(w):
    s = np.sort(np.asarray(w, float))
    n = s.size
    return 2 * np.sum(np.arange(1, n + 1) * s) / (n * s.sum()) - (n + 1) / n


D4_SCORES = {
    "active count":            lambda x: -float((x > TAU).sum()),
    "S-score, tau = tau*":     lambda x: -float((x > TAU).mean()),
    "S-score, tau = x0":       lambda x: -float((x > X0).mean()),
    "entropy, beta = x0":      lambda x: float(entropy(x, X0)[0]),
    "entropy, beta = tau*":    lambda x: float(entropy(x, TAU)[0]),
    "Gini, beta = x0":         lambda x: gini_w(np.maximum(x - X0, 0.0)),
    "ratio R(1), b = x0":      lambda x: float(np.sort(x)[::-1][0] - max(np.sort(x)[::-1][1], X0)),
    "ratio R(3), b = x0":      lambda x: float(np.sort(x)[::-1][0] - max(np.sort(x)[::-1][3], X0)),
    "gated entropy, T = 1":    lambda x: float(candidate(x, X0, 1.0)[0]),
    "published entropy":       lambda x: float(neg_entropy(np.where(x > X0, 10.0 ** (x - 9.0), 0.0))[0]),
    "P_max":                   lambda x: (lambda K: float(K.max() / K.sum()))(np.where(x > X0, 10.0 ** (x - 9.0), 0.0)),
    "gated entropy, no gate":  lambda x: float(neg_entropy(1.0 * np.logaddexp(0.0, (x - X0) / 1.0))[0]),
    "active count, reversed":  lambda x: float((x > TAU).sum()),
}
for name, f in D4_SCORES.items():
    ok = True
    for a in (6.5, 7.0, 8.0, 9.0):
        v = [f(np.r_[np.full(m, a), np.full(NREF - m, X0)]) for m in range(1, NREF + 1)]
        chain = all(v[i + 1] <= v[i] + 1e-12 for i in range(NREF - 1))
        ok &= chain and v[0] > v[-1] + 1e-12
    print(f"   {name:24s}: {'D4 holds' if ok else 'D4 FAILS'} (a = 6.5, 7, 8, 9; N = {NREF})")

print("\n10. Proposition 3: on a fixed panel, turning a zero weight into a positive weight no larger than")
print("    the smallest positive weight never raises -H, the largest share (P_max) or the Gini;")
print("    random check over 20,000 weight vectors")
rng3 = np.random.default_rng(3)
viol = {"-H": 0, "largest share": 0, "Gini": 0}
for _ in range(20000):
    r = int(rng3.integers(1, 15)); w = rng3.lognormal(0.0, 1.5, r)
    add = rng3.uniform(0.0, w.min()); n = r + int(rng3.integers(1, 10))
    w0 = np.r_[w, np.zeros(n - r)]; w1 = np.r_[w, add, np.zeros(n - r - 1)]
    viol["-H"] += neg_entropy(w1)[0] > neg_entropy(w0)[0] + 1e-12
    viol["largest share"] += w1.max() / w1.sum() > w0.max() / w0.sum() + 1e-12
    viol["Gini"] += gini_w(w1) > gini_w(w0) + 1e-12
print("    violations:", {k: int(v) for k, v in viol.items()})
H2 = lambda w: float(-neg_entropy(np.asarray(w, float))[0])
print(f"    why the weakest-weight condition is needed: weights (10, fifty 1s) have H = {H2([10] + [1] * 50):.6f} bits;")
print(f"    adding another 10, tied with the strongest, gives {H2([10, 10] + [1] * 50):.6f} bits (-H rises)")
print(f"    a concentration score can rank more active targets as more selective: H(1, 1) = {H2([1, 1]):.6f} bits,")
print(f"    H(1000, 1, 1) = {H2([1000, 1, 1]):.6f} bits")
print("    adding one active kinase (value 12) to one target at 9 and fifty at 5.1 (x0 = 5):")
base = np.r_[9.0, np.full(50, 5.1)]
for name in ("gated entropy, T = 1", "entropy, beta = x0", "P_max", "active count"):
    f = D4_SCORES[name]
    print(f"      {name:22s}: {f(base):+.3f} -> {f(np.r_[base, 12.0]):+.3f}")
print("    D4 depends on the floor treatment: with floor weight 1 and hit weight 2,")
print(f"      -H of (2, 1, 1, 1) = {-H2([2, 1, 1, 1]):.6f}; of (2, 2, 1, 1) = {-H2([2, 2, 1, 1]):.6f} (rises: broader is not more even)")
print(f"    a newly strongest target need not raise a concentration score: -H of (1) = {-H2([1]):.6f}; of (1, 2) = {-H2([1, 2]):.6f}")
print("\n11. D4 with floor-coded pairs kept at the floor (association-constant weights 10^(x - 9), floor 10^(x0 - 9))")
for N in (30, 172, 217, 363):
    parts = []
    for a in (6.01, 7.0, 8.0, 9.0):
        prof = lambda m: np.r_[np.full(m, a), np.full(N - m, X0)]
        ent = [float(neg_entropy(10.0 ** (prof(m) - 9.0))[0]) for m in range(1, N + 1)]
        pmx = [(lambda K: float(K.max() / K.sum()))(10.0 ** (prof(m) - 9.0)) for m in range(1, N + 1)]
        e_ok = all(ent[i + 1] <= ent[i] + 1e-12 for i in range(N - 1)) and ent[0] > ent[-1]
        p_ok = all(pmx[i + 1] <= pmx[i] + 1e-15 for i in range(N - 1)) and pmx[0] > pmx[-1]
        parts.append(f"a={a:g}: entropy {'holds' if e_ok else 'FAILS'}, P_max {'holds' if p_ok else 'FAILS'}")
    print(f"   N = {N:3d}: " + "; ".join(parts))
N, a = 363, 7.0
e1, e2 = (float(neg_entropy(10.0 ** (np.r_[np.full(m, a), np.full(N - m, X0)] - 9.0))[0]) for m in (1, 2))
print(f"   example, N = 363, targets at pKd 7: published entropy (floor kept) {e1:+.4f} with one target, {e2:+.4f} with two")

print("\n12. D4 for the entropy and Gini over the Anastassiadis baseline sweep (x0 = 20%, tau* = 50%, N = 291);")
print("    beta < x0 gives floor-coded pairs the positive weight x0 - beta")
XA, NA = 20.0, 291
for beta in (0.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0):
    parts = []
    for name, f in (("entropy", lambda w: float(neg_entropy(w)[0])), ("Gini", gini_w)):
        fails = []
        for a in (55.0, 60.0, 70.0, 80.0, 90.0, 100.0):
            v = [f(np.maximum(np.r_[np.full(m, a), np.full(NA - m, XA)] - beta, 0.0)) for m in range(1, NA + 1)]
            if not (all(v[i + 1] <= v[i] + 1e-12 for i in range(NA - 1)) and v[0] > v[-1] + 1e-12):
                fails.append(f"{a:g}")
        parts.append(f"{name} " + (f"FAILS at a = {', '.join(fails)}%" if fails else "holds"))
    print(f"   beta = {beta:4.0f}%: " + "; ".join(parts) + "  (a tested: 55, 60, 70, 80, 90, 100%)")
g1, g2 = (gini_w(np.r_[np.full(m, 60.0), np.full(NA - m, XA)]) for m in (1, 2))
h1, h2 = (float(neg_entropy(np.r_[np.full(m, 60.0), np.full(NA - m, XA)])[0]) for m in (1, 2))
print(f"   example, beta = 0, targets at 60%: Gini {g1:.5f} with one target, {g2:.5f} with two; "
      f"entropy {h1:+.4f} with one, {h2:+.4f} with two")
