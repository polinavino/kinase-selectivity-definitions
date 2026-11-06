"""
Monotonicity under weak off-target addition (the paper's D3) for the entropy and
the candidate.

As stated in Methods: adding a kinase with affinity in (x0, tau*] must not raise
a compound's selectivity score. This script supports Methods (the proposition)
and Results, "Monotonicity Under Weak-Target Addition" (SI Section S4):

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
