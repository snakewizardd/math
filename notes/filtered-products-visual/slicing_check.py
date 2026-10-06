#!/usr/bin/env python3
"""Numeric check of the simplex-slicing lemma (Section 2 of the paper).

rho_r(x)_a = ell * max(0, min(s_{a+1}, r/ell) - max(s_a, (r-1)/ell)),
with s_a = x_0 + ... + x_{a-1}. The lemma claims:

  (i)   each rho_r(x) is a point of the simplex,
  (ii)  at most j of rho_1(x), ..., rho_ell(x) are non-vertices,
  (iii) alpha_* rho_r(x) = rho_r(alpha_* x) for every monotone
        alpha: [j] -> [k], surjections and empty fibres included.

(ii) is what makes the union with >= n restricted factors collapse in
simplicial degrees j < n. This script samples random points, including
points with breaks exactly on slice endpoints, and checks all three.
Needs numpy.
"""

import itertools
import numpy as np

rng = np.random.default_rng(11)
TOL = 1e-12


def rho(x, ell):
    s = np.concatenate([[0.0], np.cumsum(x)])
    s[-1] = 1.0
    out = np.empty((ell, len(x)))
    for r in range(1, ell + 1):
        lo, hi = (r - 1) / ell, r / ell
        out[r - 1] = ell * np.maximum(0.0, np.minimum(s[1:], hi) - np.maximum(s[:-1], lo))
    return out


def push(alpha, x, k):
    y = np.zeros(k + 1)
    for a, xa in enumerate(x):
        y[alpha[a]] += xa
    return y


def is_vertex(p):
    return np.isclose(p.max(), 1.0, atol=1e-9)


def monotone_maps(j, k):
    return list(itertools.combinations_with_replacement(range(k + 1), j + 1))


def sample(j, ell, kind):
    if kind == "random":
        return rng.dirichlet(np.ones(j + 1))
    breaks = np.sort(rng.integers(0, ell + 1, size=j)) / ell
    if kind == "mixed":
        m = rng.integers(0, j + 1)
        breaks[:m] = rng.random(m)
        breaks = np.sort(breaks)
    return np.diff(np.concatenate([[0.0], breaks, [1.0]]))


worst_sum = worst_neg = worst_nat = 0.0
max_ratio = 0.0
checked = nat_checked = 0
for j in range(0, 5):
    for ell in (1, 2, 3, 5, 8, 13):
        maps = {k: monotone_maps(j, k) for k in range(0, 4)}
        for kind in ("random", "endpoints", "mixed"):
            for _ in range(60):
                x = sample(j, ell, kind)
                R = rho(x, ell)
                worst_sum = max(worst_sum, float(np.abs(R.sum(axis=1) - 1).max()))
                worst_neg = min(worst_neg, float(R.min()))
                nonv = sum(not is_vertex(p) for p in R)
                if j > 0:
                    max_ratio = max(max_ratio, nonv / j)
                elif nonv:
                    max_ratio = float("inf")
                checked += 1
                for k, ms in maps.items():
                    for alpha in ms:
                        lhs = np.array([push(alpha, R[r], k) for r in range(ell)])
                        rhs = rho(push(alpha, x, k), ell)
                        worst_nat = max(worst_nat, float(np.abs(lhs - rhs).max()))
                        nat_checked += 1

print(f"points checked: {checked}, (point, alpha) pairs for naturality: {nat_checked}")
print(f"(i)   max |sum - 1| = {worst_sum:.2e}, min coordinate = {worst_neg:.2e}")
print(f"(ii)  max (#non-vertices / j) = {max_ratio:.3f}  (lemma: <= 1)")
print(f"(iii) max |alpha_* rho - rho alpha_*| = {worst_nat:.2e}")
ok = worst_sum < 1e-9 and worst_neg > -TOL and max_ratio <= 1 and worst_nat < 1e-9
print("all claims hold on the sample" if ok else "a claim FAILED on the sample")
