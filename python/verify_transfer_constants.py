#!/usr/bin/env python3
"""
verify_transfer_constants.py

Independent recomputation of Phi_y and Phi_hh (Theorem 3.3 of the paper) by
direct finite differences on the (y, h) reduced critical system, WITHOUT
reusing the closed-form derivation given in the paper's proof.

This is a genuine adversarial cross-check, not a restatement: it works
directly in the y = x^(m-1), h = s(x)/x coordinates used by the transfer
theorem, holds one variable fixed while perturbing the other, and compares
the resulting numerical derivatives against the values tabulated in the
paper's Table 3.

Note (FR) : recalcul independant de Phi_y et Phi_hh (Theoreme 3.3 de
l'article) par differences finies directes sur le systeme critique reduit
(y, h), sans reprendre la derivation en forme close de la preuve. C'est
une verification adversariale genuine, pas une reformulation : le calcul
travaille dans les coordonnees y = x^(m-1), h = s(x)/x du theoreme de
transfert, en fixant une variable pendant qu'on perturbe l'autre, et
compare le resultat aux valeurs du Tableau 3 de l'article.

Method:
  - Phi_h (should equal 1 at criticality, sanity check) and Phi_hh: perturb
    h around h_c = tau_m/rho_m, at fixed x = rho_m (i.e. fixed y = y_c).
  - Phi_y: perturb x around rho_m (i.e. perturb y = x^(m-1) away from y_c),
    holding the FREE h variable fixed at h_c -- this is what the paper's
    exp-schema transfer theorem calls Phi_y, distinct from the trivial
    d/dx derivative at a moving self-consistent point.

Reference: F. G. Speyser, "Strict Monotonicity and a Lambert-W Asymptotic
for Growth Rates of Non-Plane Strict m-Gonal Cacti", Theorem 3.3, Table 3.

Author: Frederic G. Speyser
Run with: python3 verify_transfer_constants.py
"""
import math
from critical_point import evaluate_s_tree, find_rho


def kernel_terms(m, x, depth, h_free=None):
    """i=1 term of the kernel (with h possibly frozen at h_free) plus the
    i>=2 tail, at the point x, using the true self-consistent series values
    for every i>=2 argument."""
    svals, xs = evaluate_s_tree(m, x, depth=depth, niter=150)
    n = m - 1
    odd = (m % 2 == 1)
    y = x ** n

    def kappa(h_val, h_y2):
        if odd:
            return 0.5 * (h_val ** (m - 1) + h_y2 ** ((m - 1) / 2))
        return 0.5 * (h_val ** (m - 1) + h_val * h_y2 ** ((m - 2) / 2))

    h_true = svals[1] / x if x > 0 else 1.0
    idx2 = 2 * n
    h_y2 = svals[idx2] / x ** idx2 if idx2 <= depth else 0.0
    h_used = h_free if h_free is not None else h_true
    i1 = y * kappa(h_used, h_y2)

    tail = 0.0
    i = 2
    while i * n <= depth:
        yi = y ** i
        h_yi = svals[i * n] / x ** (i * n)
        idx2i = 2 * i * n
        h_y2i = svals[idx2i] / x ** idx2i if idx2i <= depth else 0.0
        tail += yi * kappa(h_yi, h_y2i) / i
        i += 1
    return i1 + tail


def verify(m, depth=60):
    rho, tau, _ = find_rho(m, lo=0.5, hi=0.9, steps=60)
    h_c = tau / rho

    def Phi_of_h(h_val):
        return math.exp(kernel_terms(m, rho, depth, h_free=h_val))

    eps_h = 1e-6
    Phi_h = (Phi_of_h(h_c + eps_h) - Phi_of_h(h_c - eps_h)) / (2 * eps_h)
    Phi_hh = (Phi_of_h(h_c + eps_h) - 2 * Phi_of_h(h_c) + Phi_of_h(h_c - eps_h)) / eps_h ** 2

    def Phi_of_x(xv):
        return math.exp(kernel_terms(m, xv, depth, h_free=h_c))

    n = m - 1
    eps_x = rho * 1e-5
    Phi_y_via_x = (Phi_of_x(rho + eps_x) - Phi_of_x(rho - eps_x)) / (2 * eps_x)
    dydx = n * rho ** (n - 1)
    Phi_y = Phi_y_via_x / dydx

    return rho, tau, Phi_h, Phi_hh, Phi_y


TABLE_3 = {
    5: dict(Phi_hh=2.8768, Phi_y=3.6076),
    6: dict(Phi_hh=3.6484, Phi_y=3.4283),
    7: dict(Phi_hh=4.8273, Phi_y=3.0917),
    8: dict(Phi_hh=5.5898, Phi_y=3.0607),
}

if __name__ == "__main__":
    print(f"{'m':>3} {'Phi_h (~=1?)':>13} {'Phi_hh':>10} {'expected':>10} "
          f"{'Phi_y':>10} {'expected':>10}")
    for m in [5, 6, 7, 8]:
        rho, tau, Phi_h, Phi_hh, Phi_y = verify(m)
        exp = TABLE_3[m]
        print(f"{m:3d} {Phi_h:13.6f} {Phi_hh:10.4f} {exp['Phi_hh']:10.4f} "
              f"{Phi_y:10.4f} {exp['Phi_y']:10.4f}")
