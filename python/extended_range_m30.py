#!/usr/bin/env python3
"""
extended_range_m30.py

Extends the critical-data table (rho_m, tau_m, mu_m) and the check of
Theorem 3.3's two-term equivalent (Table 4 of the paper) beyond m=12, up
to m=30, using the already-verified critical_point.py solver.

Also confirms that the residual (mu_m - two-term approximation), rescaled
by m^2/log(m), stays within a narrow band across the whole range -- the
numerical evidence, quoted in the paper, that the stated O((log m)/m^2)
error term is the right order, not merely a safe upper bound.

Note (FR) : etend le tableau des donnees critiques (rho_m, tau_m, mu_m)
et la verification de l'equivalent a deux termes du Theoreme 3.3 (Tableau
4 de l'article) au-dela de m=12, jusqu'a m=30, en utilisant le solveur
deja verifie critical_point.py. Confirme aussi que le residu (mu_m moins
l'approximation a deux termes), remis a l'echelle par m^2/log(m), reste
dans une bande etroite sur toute la plage -- la preuve numerique,
mentionnee dans l'article, que le terme d'erreur annonce O((log m)/m^2)
est le bon ordre, pas juste une borne large et prudente.

Reference: F. G. Speyser, "Strict Monotonicity and a Lambert-W Asymptotic
for Growth Rates of Non-Plane Strict m-Gonal Cacti", Table 3, Table 4,
Figure 3.

Author: Frederic G. Speyser
Run with: python3 extended_range_m30.py
"""
import math
from critical_point import find_rho

OMEGA = 0.2784645427610738  # W(e^{-1}), Lambert function principal value


def lambda_m(m):
    n = m - 1
    return (math.log(n / 2) + 1 + OMEGA) / n


def main():
    m_values = [5, 6, 7, 8, 9, 10, 11, 12, 15, 18, 20, 22, 25, 28, 30]

    print(f"{'m':>4} {'rho_m':>10} {'tau_m':>10} {'mu_m':>10}")
    results = {}
    for m in m_values:
        lo, hi = 0.55, min(0.55 + 0.02 * m, 0.995)
        rho, tau, _ = find_rho(m, lo=lo, hi=hi, steps=80)
        mu = 1 / rho
        results[m] = (rho, tau, mu)
        print(f"{m:4d} {rho:10.6f} {tau:10.6f} {mu:10.6f}")

    print()
    print("Check against the two-term approximation of Theorem 5.3:")
    print(f"{'m':>4} {'1-rho_m':>10} {'lambda_m':>10} {'approx':>10} "
          f"{'residual':>10} {'residual*m^2/log(m)':>20}")
    for m in m_values:
        rho, tau, mu = results[m]
        lam = lambda_m(m)
        approx = 1 + lam - 0.5 * lam ** 2
        residual = mu - approx
        rescaled = residual * m ** 2 / math.log(m)
        print(f"{m:4d} {1 - rho:10.6f} {lam:10.6f} {approx - 1:10.6f} "
              f"{residual:10.6f} {rescaled:20.4f}")


if __name__ == "__main__":
    main()
