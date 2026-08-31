#!/usr/bin/env python3
"""
newton_rho.py

A second, independent ratio-test estimate of rho_m (different truncation
and refinement strategy from compute_rho.py), used as a cross-check
against critical_point.py and against the known values quoted in the
paper for m=5..9. Includes the same n^(-3/2) polynomial correction as
compute_rho.py; without it, this script also showed a systematic bias
of order 10^-2, now removed.

Note (FR) : une seconde estimation, independante, de rho_m par test de
ratio (strategie de troncature et d'affinement differente de
compute_rho.py), utilisee en verification croisee de critical_point.py
et des valeurs deja connues pour m=5..9. Inclut la meme correction
polynomiale en n^(-3/2) que compute_rho.py ; sans elle, ce script
presentait aussi un biais systematique d'ordre 10^-2, desormais retire.

Author: Frederic G. Speyser
Run with: python3 newton_rho.py
"""
import math
import numpy as np
from numpy.polynomial.polynomial import polypow, polymul


def pow_pad(p, e, N):
    if e == 0:
        out = np.zeros(N + 1)
        out[0] = 1.0
        return out
    q = polypow(p, int(e))
    out = np.zeros(N + 1)
    L = min(N + 1, len(q))
    out[:L] = q[:L]
    return out


def series_s(m, N, niter=60):
    s = np.zeros(N + 1)
    s[1] = 1.0
    odd = (m % 2 == 1)
    for it in range(niter):
        sm1 = pow_pad(s, m - 1, N)
        s2 = np.zeros(N + 1)
        for k in range((N // 2) + 1):
            s2[2 * k] = s[k]
        if odd:
            refl = pow_pad(s2, (m - 1) // 2, N)
        else:
            mid = pow_pad(s2, (m - 2) // 2, N)
            prod = polymul(s, mid)
            refl = np.zeros(N + 1)
            L = min(N + 1, len(prod))
            refl[:L] = prod[:L]
        KC = 0.5 * (sm1 + refl)
        E = np.zeros(N + 1)
        max_i = min(N, 40)
        for i in range(1, max_i + 1):
            max_t = N // i
            for t in range(1, max_t + 1):
                E[i * t] += KC[t] / i
        F = np.zeros(N + 1)
        F[0] = 1.0
        for n in range(1, N + 1):
            acc = 0.0
            for k in range(1, n + 1):
                acc += k * E[k] * F[n - k]
            F[n] = acc / n
        news = np.zeros(N + 1)
        news[1:] = F[:N]
        news = np.maximum(news, 0.0)
        err = np.max(np.abs(news - s))
        s = news
        if err < 1e-15:
            break
    return s


def rho_estimates(s, m):
    step = m - 1
    idxs = [n for n in range(1, len(s), step) if s[n] > 0]
    est = []
    for i in range(len(idxs) - 1):
        a, b = idxs[i], idxs[i + 1]
        if s[b] > 0:
            # correct for the n^(-3/2) polynomial factor before extracting rho
            corrected = (s[a] / s[b]) * (a / b) ** 1.5
            est.append(corrected ** (1.0 / (b - a)))
    return idxs, est


def main():
    known = {5: 0.604765, 6: 0.633235, 7: 0.669930, 8: 0.690268, 9: 0.71494}
    print(f"{'m':>3} {'N':>5} {'kmax':>5} {'rho_last':>10} {'rho_med5':>10} {'known':>10} {'bias':>10} {'1/rho':>10}")
    for m in range(5, 13):
        kmax = 28 if m <= 8 else (22 if m <= 10 else 18)
        N = 1 + (m - 1) * kmax
        N = min(N, 260)
        s = series_s(m, N, niter=45)
        idxs, est = rho_estimates(s, m)
        last = est[-1] if est else float("nan")
        med = float(np.median(est[-5:])) if len(est) >= 5 else last
        kn = known.get(m, float("nan"))
        bias = last - kn if kn == kn else float("nan")
        print(f"{m:3d} {N:5d} {len(idxs)-1:5d} {last:10.6f} {med:10.6f} {kn:10.6f} {bias:10.6f} {1/last:10.4f}")


if __name__ == "__main__":
    main()
