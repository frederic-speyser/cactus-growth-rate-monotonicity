#!/usr/bin/env python3
"""
compute_rho.py

High-order series expansion of s_m(x), plus a ratio-test estimate of
rho_m, cross-checked against the closed-form tau_m for odd m (Theorem
3.1 of the paper). The ratio estimate includes the n^(-3/2) polynomial
correction implied by the transfer theorem (Theorem 3.3): omitting it
introduces a systematic bias of order 10^-2, corrected here.

Note (FR) : developpement en serie a haut degre de s_m(x), puis
estimation de rho_m par test de ratio, comparee a la forme close de
tau_m pour m impair (Theoreme 3.1 de l'article). L'estimation par ratio
inclut la correction polynomiale en n^(-3/2) exigee par le theoreme de
transfert (Theoreme 3.3) : l'omettre introduit un biais systematique
d'ordre 10^-2, corrige ici.

Author: Frederic G. Speyser
Run with: python3 compute_rho.py
"""
import math
import numpy as np
from numpy.polynomial.polynomial import polypow, polymul


def series_s(m, N, niter=40):
    """Coefficients of s_m up to degree N (inclusive)."""
    s = np.zeros(N + 1)
    s[1] = 1.0  # start from the atom
    odd = (m % 2 == 1)
    for _ in range(niter):
        # KC = 1/2 (s^{m-1} + reflection)
        def pow_pad(p, e):
            q = polypow(p, e)
            out = np.zeros(N + 1)
            L = min(N + 1, len(q))
            out[:L] = q[:L]
            return out

        sm1 = pow_pad(s, m - 1)
        s2 = np.zeros(N + 1)
        for k in range((N // 2) + 1):
            s2[2 * k] = s[k]
        if odd:
            refl = pow_pad(s2, (m - 1) // 2)
        else:
            mid = pow_pad(s2, (m - 2) // 2)
            refl = np.zeros(N + 1)
            prod = polymul(s, mid)
            L = min(N + 1, len(prod))
            refl[:L] = prod[:L]
        KC = 0.5 * (sm1 + refl)
        if len(KC) < N + 1:
            KC = np.pad(KC, (0, N + 1 - len(KC)))
        KC = KC[: N + 1]
        # E = sum_{i>=1} KC(x^i)/i
        E = np.zeros(N + 1)
        for i in range(1, N + 1):
            max_t = N // i
            if max_t < 1:
                break
            for t in range(1, max_t + 1):
                E[i * t] += KC[t] / i
        # exp(E) via recurrence: F' = F E'
        F = np.zeros(N + 1)
        F[0] = 1.0
        for n in range(1, N + 1):
            F[n] = sum(k * E[k] * F[n - k] for k in range(1, n + 1)) / n
        news = np.zeros(N + 1)
        # s = x * F
        news[1:] = F[:N]
        # keep nonnegative numerical noise cleanup
        news = np.maximum(news, 0.0)
        if np.max(np.abs(news - s)) < 1e-16:
            s = news
            break
        s = news
    return s


def rho_from_ratios(s, m):
    """Estimate rho from [x^{n}] / [x^{n+(m-1)}] ^{1/(m-1)}."""
    step = m - 1
    idxs = [n for n in range(1, len(s), step) if s[n] > 0]
    if len(idxs) < 4:
        return float("nan"), []
    estimates = []
    for a, b in zip(idxs[-8:-1], idxs[-7:]):
        if s[a] > 0 and s[b] > 0:
            # s_n ~ C*rho^-n*n^-3/2, so s_a/s_b ~ rho^(b-a) * (b/a)^(3/2)
            # correcting for the polynomial factor before extracting rho:
            corrected = (s[a] / s[b]) * (a / b) ** 1.5
            estimates.append(corrected ** (1.0 / (b - a)))
    return float(np.median(estimates)), estimates


def main():
    results = []
    W = 0.2784645427610738  # W(1/e)
    print(f"{'m':>4} {'rho_hat':>12} {'tau_odd':>12} {'1-rho':>12} {'log(m/2)/m':>12} {'pred':>12} {'ratio':>10}")
    for m in range(5, 21):
        # need enough blocks: n = 1+(m-1)k, take k up to ~25
        N = 1 + (m - 1) * (22 if m <= 12 else 16)
        N = min(N, 280)
        niter = 50 if m <= 12 else 35
        s = series_s(m, N, niter=niter)
        rho, ests = rho_from_ratios(s, m)
        tau = (2.0 / (m - 1)) ** (1.0 / (m - 1))
        pred = (math.log((m - 1) / 2.0) + 1.0 + W) / (m - 1)
        line = (m, rho, tau, 1 - rho, math.log(m / 2.0) / m, pred, (1 - rho) / pred if pred else float("nan"))
        results.append(line)
        print(f"{m:4d} {rho:12.6f} {tau:12.6f} {1-rho:12.6f} {math.log(m/2)/m:12.6f} {pred:12.6f} {(1-rho)/pred:10.4f}")
        # leading coeffs
        step = m - 1
        coeffs = [(n, s[n]) for n in range(1, min(len(s), 1 + step * 6), step)]
        print("   coeffs", ["%.4g" % c[1] for c in coeffs[:6]])


if __name__ == "__main__":
    main()
