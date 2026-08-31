#!/usr/bin/env python3
"""
critical_point.py

Solves the critical system (rho_m, tau_m) directly, by damped Picard
iteration on the tree of values s(x^j), j = 1..depth, rather than by a
ratio test on series coefficients. This is the reference method used
throughout the repository: its output matches Table 3 of the paper to
5-6 significant figures, and every other script here is checked against
it.

Note (FR) : ce script resout le systeme critique (rho_m, tau_m) par
iteration de Picard directe sur l'arbre des valeurs s(x^j), sans passer
par un test de ratio sur les coefficients de la serie. C'est la methode
de reference du depot : ses resultats concordent avec le Tableau 3 de
l'article a 5-6 chiffres significatifs, et tous les autres scripts sont
verifies contre lui.

Author: Frederic G. Speyser
Run with: python3 critical_point.py
"""
import math


def evaluate_s_tree(m, x, depth=24, niter=80):
    """Approximate s(x^j) for j=1..depth by simultaneous Picard."""
    # xs[j] = x^j, svals[j] ~ s(x^j)
    xs = [0.0] + [x ** j for j in range(1, depth + 1)]
    svals = xs[:]  # start from the atom s(z)~z
    odd = (m % 2 == 1)

    def KC_at(j, svals, xs):
        # z = x^j, need s(z) and s(z^2)=s(x^{2j})
        sz = svals[j]
        j2 = 2 * j
        sz2 = svals[j2] if j2 <= depth else xs[j] ** 2  # s(w)~w for tiny w
        if odd:
            return 0.5 * (sz ** (m - 1) + sz2 ** ((m - 1) / 2))
        return 0.5 * (sz ** (m - 1) + sz * (sz2 ** ((m - 2) / 2)))

    for _ in range(niter):
        news = [0.0]
        for j in range(1, depth + 1):
            expo = 0.0
            for i in range(1, depth // j + 1):
                expo += KC_at(j * i, svals, xs) / i if j * i <= depth else 0.0
            if expo > 20:
                news.append(float("inf"))
            else:
                news.append(xs[j] * math.exp(expo))
        if max(abs(news[j] - svals[j]) for j in range(1, depth + 1)) < 1e-16:
            svals = news
            break
        svals = news
    return svals, xs


def phi_y(m, x, svals):
    """partial Phi / partial y at (x, s(x)), tails frozen."""
    s = svals[1]
    s2 = svals[2]
    odd = (m % 2 == 1)
    if odd:
        dKC = 0.5 * (m - 1) * s ** (m - 2)
        KC = 0.5 * (s ** (m - 1) + s2 ** ((m - 1) / 2))
    else:
        A = 0.5 * (s2 ** ((m - 2) / 2))
        dKC = 0.5 * (m - 1) * s ** (m - 2) + A
        KC = 0.5 * (s ** (m - 1) + s * s2 ** ((m - 2) / 2))
    h = 0.0
    depth = len(svals) - 1
    for i in range(2, depth + 1):
        si = svals[i]
        si2 = svals[2 * i] if 2 * i <= depth else 0.0
        if odd:
            KCi = 0.5 * (si ** (m - 1) + si2 ** ((m - 1) / 2)) if si2 > 0 else 0.5 * si ** (m - 1)
        else:
            KCi = 0.5 * (si ** (m - 1) + si * (si2 ** ((m - 2) / 2))) if si2 > 0 else 0.5 * si ** (m - 1)
        h += KCi / i
    # Phi = x exp(KC+h), Phi_y = x exp(KC+h) * dKC = s * dKC
    return s * dKC, s, KC, h


def find_rho(m, lo=0.50, hi=0.85, steps=36):
    """Bisection on Phi_y(x,s(x)) - 1."""
    def fy(x):
        svals, _ = evaluate_s_tree(m, x)
        if any(not math.isfinite(v) for v in svals):
            return 10.0, float("nan")
        val, s, KC, h = phi_y(m, x, svals)
        if not math.isfinite(val):
            return 10.0, s
        return val - 1.0, s

    # ensure fy(lo)<0<fy(hi) typically
    for _ in range(steps):
        mid = 0.5 * (lo + hi)
        d, s = fy(mid)
        if d < 0:
            lo = mid
        else:
            hi = mid
    x = 0.5 * (lo + hi)
    d, tau = fy(x)
    return x, tau, d


def main():
    print(f"{'m':>3} {'rho':>12} {'tau':>12} {'mu':>10} {'fy-1':>12} {'tau_odd_form':>14}")
    for m in range(5, 13):
        lo, hi = 0.55, min(0.55 + 0.03 * m, 0.88)
        rho, tau, d = find_rho(m, lo=lo, hi=hi)
        form = (2.0 / (m - 1)) ** (1.0 / (m - 1))
        print(f"{m:3d} {rho:12.8f} {tau:12.8f} {1/rho:10.5f} {d:12.3e} {form:14.8f}")


if __name__ == "__main__":
    main()
