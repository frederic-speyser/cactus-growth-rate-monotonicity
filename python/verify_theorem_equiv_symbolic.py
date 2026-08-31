#!/usr/bin/env python3
"""
verify_theorem_equiv_symbolic.py

Independent symbolic verification of the algebraic chain underlying the
proof of Theorem 5.3 (the Lambert-W large-m equivalent), using SymPy.
This checks the identities themselves, not the numerical values: every
step below is an exact symbolic manipulation, distinct from the purely
numerical cross-checks (critical_point.py, verify_transfer_constants.py,
verify_critical_pari.gp) elsewhere in this repository.

Note (FR) : verification symbolique independante de la chaine algebrique
qui sous-tend la preuve du Theoreme 5.3 (l'equivalent de Lambert W pour
m grand), a l'aide de SymPy. Ceci verifie les identites elles-memes, pas
des valeurs numeriques : chaque etape ci-dessous est une manipulation
symbolique exacte, distincte des verifications purement numeriques
(critical_point.py, verify_transfer_constants.py,
verify_critical_pari.gp) presentes ailleurs dans ce depot.

Three identities are checked, matching the three key steps of the proof:

  1. The defining equation v*e^v = e^(-1) has the principal Lambert
     value Omega = W(e^(-1)) as its solution -- the substitution used to
     pass from u = 1 + e^(-u) to u = 1 + W(e^(-1)).

  2. At u = 1 + Omega, the derivative of g(u) = u - 1 - e^(-u) equals
     1 + Omega itself (not merely "nonzero", the precise value used to
     justify the linearisation step in the proof).

  3. The expansion 1 - e^(-Lambda) = Lambda - Lambda^2/2 + O(Lambda^3),
     used in the final "Transfer" step of the proof to go from
     -log(rho_m) = lambda_m + O(m^-2) to the stated equivalent.

Reference: F. G. Speyser, "Strict Monotonicity and a Lambert-W Asymptotic
for Growth Rates of Non-Plane Strict m-Gonal Cacti", Theorem 5.3.

Author: Frederic G. Speyser
Run with: python3 verify_theorem_equiv_symbolic.py
"""
import sympy as sp


def check_1_lambert_definition():
    """v*e^v = e^{-1} has the principal branch solution v = W(e^{-1})."""
    v = sp.symbols('v', real=True)
    solutions = sp.solve(sp.Eq(v * sp.exp(v), sp.exp(-1)), v)
    Omega = sp.LambertW(sp.exp(-1))
    ok = any(sp.simplify(s - Omega) == 0 for s in solutions)
    print("Check 1 -- v*e^v = e^(-1)  =>  v = W(e^(-1))")
    print(f"  sympy.solve gives: {solutions}")
    print(f"  matches Omega = W(e^-1): {ok}")
    print(f"  numerical value: Omega = {sp.N(Omega, 20)}")
    return ok


def check_2_derivative_at_fixed_point():
    """g(u) = u - 1 - e^{-u}; g'(1+Omega) should equal 1+Omega exactly,
    using the defining property Omega*e^Omega = e^{-1}, i.e.
    e^{-Omega} = Omega*e."""
    u = sp.symbols('u', real=True)
    Omega = sp.symbols('Omega', positive=True)
    g = u - 1 - sp.exp(-u)
    gprime = sp.diff(g, u)  # = 1 + e^{-u}... wait, d/du(-e^{-u}) = e^{-u}
    gprime_at = gprime.subs(u, 1 + Omega)
    # gprime_at = 1 + e^{-(1+Omega)} = 1 + e^{-1}*e^{-Omega}
    # Using the defining relation Omega*e^Omega = e^{-1}  =>  e^{-Omega} = Omega*e
    # so e^{-1}*e^{-Omega} = e^{-1}*Omega*e = Omega
    e_neg_Omega_substitute = Omega * sp.E  # from the defining relation
    gprime_at_substituted = gprime_at.subs(sp.exp(-Omega), e_neg_Omega_substitute)
    simplified = sp.simplify(gprime_at_substituted - (1 + Omega))
    print()
    print("Check 2 -- g(u) = u - 1 - e^(-u); g'(1+Omega) = 1+Omega")
    print(f"  g'(u) = {gprime}")
    print(f"  g'(1+Omega) before substitution: {gprime_at}")
    print(f"  after using Omega*e^Omega = e^(-1): {gprime_at_substituted}")
    print(f"  difference from (1+Omega), simplified: {simplified}")
    ok = simplified == 0
    # Cross-check numerically to high precision as well, independently
    Omega_num = sp.LambertW(sp.exp(-1))
    numeric_check = sp.N(1 + sp.exp(-(1 + Omega_num)) - (1 + Omega_num), 30)
    print(f"  independent high-precision numerical check: {numeric_check}")
    return ok and abs(numeric_check) < sp.Float('1e-25')


def check_3_expansion():
    """1 - e^{-Lambda} = Lambda - Lambda^2/2 + O(Lambda^3)."""
    Lambda = sp.symbols('Lambda')
    expr = 1 - sp.exp(-Lambda)
    series = sp.series(expr, Lambda, 0, 4).removeO()
    expected = Lambda - Lambda**2 / 2
    diff = sp.simplify(series - expected - Lambda**3 / 6)
    print()
    print("Check 3 -- 1 - e^(-Lambda) = Lambda - Lambda^2/2 + O(Lambda^3)")
    print(f"  sympy series (through cubic term): {series}")
    print(f"  matches Lambda - Lambda^2/2 + Lambda^3/6 exactly: {diff == 0}")
    return diff == 0


def check_4_full_chain_numeric():
    """End-to-end numerical sanity check of the full substitution chain,
    at high precision, independent of the symbolic steps above: solve
    u = 1 + e^{-u} directly by Newton's method and compare to
    1 + W(e^{-1})."""
    print()
    print("Check 4 -- end-to-end: numerically solve u=1+e^-u directly,")
    print("           compare to 1+W(e^-1), at 40-digit precision")
    u = sp.Symbol('u')
    u_num = sp.nsolve(sp.Eq(u, 1 + sp.exp(-u)), u, 1.3, prec=40)
    closed_form = sp.N(1 + sp.LambertW(sp.exp(-1)), 40)
    diff = sp.N(u_num - closed_form, 40)
    print(f"  Newton solution of u=1+e^-u:  {u_num}")
    print(f"  1 + W(e^-1) (closed form):    {closed_form}")
    print(f"  difference:                   {diff}")
    return abs(diff) < sp.Float('1e-35')


if __name__ == "__main__":
    results = [
        check_1_lambert_definition(),
        check_2_derivative_at_fixed_point(),
        check_3_expansion(),
        check_4_full_chain_numeric(),
    ]
    print()
    print("=" * 60)
    if all(results):
        print("All four checks PASSED.")
    else:
        print("SOME CHECKS FAILED -- see output above.")
    print("=" * 60)
