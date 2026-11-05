"""
Distributional consistency (boundary sensitivity) for every measure on all four
datasets: the share of compound pairs whose order reverses as a parameter moves
across B = [x0, tau*], the interval between the assay floor and the activity
threshold. A pair reverses in B if f(x; b) - f(y; b) is positive at some b in B
and negative at another. Pairs are taken among the compounds with at least one
target above tau* (those the reliability gate scores).

Parameters swept over B: the S-score's threshold tau, the entropy's and Gini's
baseline beta, the candidate's emphasis baseline beta (gate fixed at x0) and its
gate x0 (beta pinned to the gate, smoothing width T held at tau* - x0 of the fixed
floor), and the ratio's clipping level b in R_b(k) = x_(1) - max(x_(k+1), b), which
replaces any off-target affinity below b, measured or missing, by b (each k).

Aggregate rule: a measure satisfies D2 only if at most eps of pairs reverse on every
dataset, and, for the ratio, for every k.

The S-score is a step function of tau, so it is evaluated exactly, at every
observed affinity in B and between consecutive ones. The other parameters are
evaluated on a 201-point grid; the same calculation on a 101-point grid is
printed as a convergence check.
"""
import script_logging; script_logging.capture(__file__)
import os
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
EPS = 0.01


def neg_entropy(w):
    tot = w.sum(1, keepdims=True)
    p = w / np.where(tot <= 0, 1.0, tot)
    return np.where(p > 0, p * np.log2(np.where(p > 0, p, 1.0)), 0.0).sum(1)


def s_score(M, t):
    return -(M > t).mean(1)


def entropy(M, b):
    return neg_entropy(np.maximum(M - b, 0.0))


def gini(M, b):
    s = np.sort(np.maximum(M - b, 0.0), 1)
    n = s.shape[1]
    t = s.sum(1)
    g = (2 * (np.arange(1, n + 1) * s).sum(1)) / (n * np.where(t == 0, 1, t)) - (n + 1) / n
    return np.where(t == 0, 0.0, g)


def ratio(M, k, floor):
    s = -np.sort(-M, 1)
    return s[:, 0] - np.maximum(s[:, k], floor)


def candidate(M, x0, T, b):
    return neg_entropy(np.where(M > x0, T * np.logaddexp(0.0, (M - b) / T), 0.0))


def reversals(S):
    """Fraction of compound pairs whose score difference changes sign across the
    rows of S (parameter values x compounds)."""
    G, n = S.shape
    chunk = max(1, int(4e7 // (G * n)))
    count = 0
    for i0 in range(0, n, chunk):
        d = S[:, i0:i0 + chunk, None] - S[:, None, :]
        rev = (d > 0).any(0) & (d < 0).any(0)
        rows = np.arange(i0, min(i0 + chunk, n))[:, None]
        count += int((rev & (np.arange(n)[None, :] > rows)).sum())
    return count / (n * (n - 1) / 2)


aff = pd.read_csv(f"{HERE}/davis_affinity.csv")
DATA = {  # matrix, floor x0, activity threshold tau*
    "Davis":         (aff.pivot(index="Drug_Index", columns="Protein_Index", values="Affinity").fillna(5.0).values,
                      5.0, 6.0),
    "Klaeger":       (pd.read_csv(f"{HERE}/klaeger_matrix.csv", index_col=0).values, 5.0, 6.0),
    "Anastassiadis": (pd.read_csv(f"{HERE}/anastassiadis_matrix.csv", index_col=0).values, 20.0, 50.0),
    "Metz":          (pd.read_csv(f"{HERE}/metz_matrix.csv", index_col=0).values, 4.0, 6.0),
}


def sweeps(M, x0, tau, n_grid):
    B = np.linspace(x0, tau, n_grid)
    T = tau - x0
    v = np.unique(M[(M >= x0) & (M <= tau)])
    tau_exact = np.unique(np.concatenate([[x0, tau], v, (v[:-1] + v[1:]) / 2]))
    rows = {
        "S-score (tau, exact)":      [s_score(M, t) for t in tau_exact],
        "entropy (beta)":            [entropy(M, b) for b in B],
        "Gini (beta)":               [gini(M, b) for b in B],
        "candidate (emphasis beta)": [candidate(M, x0, T, b) for b in B],
        "candidate (gate x0)":       [candidate(M, b, T, b) for b in B],
    }
    rows.update({f"ratio (clipping level, k={k})": [ratio(M, k, b) for b in B] for k in range(1, 6)})
    return rows


print(f"Distributional consistency over B = [x0, tau*]: fraction of scored-compound pairs reversing; "
      f"pass if <= {EPS:.0%}. Continuous parameters on 201 points (101 points in parentheses).")
largest = {}
shift = 0.0
for name, (M, x0, tau) in DATA.items():
    keep = (M > tau).any(1)
    fine = sweeps(M, x0, tau, 201)
    coarse = sweeps(M, x0, tau, 101)
    print(f"\n{name}: B = [{x0:g}, {tau:g}], {keep.sum()} scored compounds")
    for label in fine:
        f = reversals(np.array(fine[label])[:, keep])
        c = reversals(np.array(coarse[label])[:, keep])
        if "exact" not in label:
            shift = max(shift, abs(f - c))
        key = label.split(", k=")[0] + (")" if ", k=" in label else "")
        largest[key] = max(largest.get(key, 0.0), f)
        extra = "" if "exact" in label else f" ({100 * c:5.2f}%)"
        print(f"  {label:30s} {100 * f:6.2f}%{extra}   {'PASS' if f <= EPS else 'fail'}")

print(f"\nLargest change from 101 to 201 grid points: {100 * shift:.2f} percentage points")
print("Largest fraction over the four datasets (and over k for the ratio), per swept parameter:")
for key, v in largest.items():
    print(f"  {key:30s} {100 * v:6.2f}%")
swept = [v for k, v in largest.items() if k != "candidate (emphasis beta)"]
print(f"Only the candidate's emphasis baseline stays within eps on every dataset "
      f"({100 * largest['candidate (emphasis beta)']:.2f}%); every other swept parameter reaches at least "
      f"{100 * min(swept):.2f}% on at least one dataset.")
