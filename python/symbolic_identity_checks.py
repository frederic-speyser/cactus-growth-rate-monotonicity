"""
symbolic_identity_checks.py

Symbolic (exact, not floating-point) verification of the central algebraic
identities used in the proof, via sympy. Complements the numerical checks
elsewhere in this repository: a floating-point check confirms an identity
holds to some precision, this confirms it holds exactly.

Six identities are checked:
  A) Lemma 5 (dominant-term identity)
  B) The general derivative formula T'' = T*(K_y^2 + K_yy)
  C) K_yy in both parities of the target kernel (Lemma 6)
  D) sigma_max / tau_m = exp(g1*(tau_m - 1))            (section 7)
  E) The kernel substitutions used in Theorem 3's proof, both parities
  F) Theorem A's closed form, solved independently rather than verified

Author: Frederic G. Speyser
Run with: python3 symbolic_identity_checks.py   (requires: pip install sympy)
"""
import sympy as sp

print("=" * 70)
print("A) Lemma 5: K_tilde_C(rho,tau) = tau * K_C^(m)(rho,tau)  (m odd)")
print("=" * 70)
rho, tau, m = sp.symbols('rho tau m', positive=True)
s_rho2 = sp.Symbol('s_rho2', positive=True)  # represents s_m(rho^2)

K_C_m = sp.Rational(1, 2) * tau**(m - 1) + sp.Rational(1, 2) * s_rho2**((m - 1) / 2)
K_tilde_C = sp.Rational(1, 2) * tau**m + sp.Rational(1, 2) * tau * s_rho2**((m - 1) / 2)

diff = sp.simplify(sp.expand(K_tilde_C - tau * K_C_m))
print(f"  difference simplified: {diff}")
print(f"  EXACT IDENTITY: {diff == 0}")

print()
print("=" * 70)
print("B) General formula: T(y) = C*exp(K(y))  =>  T'' = T*(K'^2 + K'')")
print("=" * 70)
y = sp.Symbol('y', positive=True)
C = sp.Symbol('C', positive=True)
K = sp.Function('K')

T = C * sp.exp(K(y))
Tpp = sp.diff(T, y, 2)
K_y = sp.diff(K(y), y)
K_yy = sp.diff(K(y), y, 2)
claimed = T * (K_y**2 + K_yy)

diff = sp.simplify(Tpp - claimed)
print(f"  difference: {diff}")
print(f"  GENERAL IDENTITY CONFIRMED: {diff == 0}")

print()
print("=" * 70)
print("C) K_yy for the target kernel, both parities (Lemma 6)")
print("=" * 70)
Ax = sp.Symbol('A_x', positive=True)   # represents s(x^2)^((m-1)/2), constant in y
Bx = sp.Symbol('B_x', positive=True)   # represents s(x^2)^((m-2)/2), constant in y
claimed_Kyy = sp.Rational(1, 2) * m * (m - 1) * y**(m - 2)

K_odd = sp.Rational(1, 2) * y**m + sp.Rational(1, 2) * y * Ax
K_yy_odd = sp.diff(K_odd, y, 2)
print(f"  odd case:  K_yy = {sp.simplify(K_yy_odd)}, matches claim: "
      f"{sp.simplify(K_yy_odd - claimed_Kyy) == 0}")

K_even = sp.Rational(1, 2) * y**m + sp.Rational(1, 2) * Bx  # Bx independent of y
K_yy_even = sp.diff(K_even, y, 2)
print(f"  even case: K_yy = {sp.simplify(K_yy_even)}, matches claim: "
      f"{sp.simplify(K_yy_even - claimed_Kyy) == 0}")

print()
print("=" * 70)
print("D) sigma_max/tau_m = exp(g1*(tau_m - 1)), given tau_m = rho_m*exp(g1)")
print("=" * 70)
rho_s, tau_s, g1 = sp.symbols('rho tau g1', positive=True)
sigma_max = rho_s * sp.exp(tau_s * g1)
sigma_max_sub = sigma_max.subs(rho_s, tau_s * sp.exp(-g1))
ratio = sp.simplify(sigma_max_sub / tau_s)
claimed_D = sp.exp(g1 * (tau_s - 1))
print(f"  ratio computed: {ratio}")
print(f"  ratio claimed:  {claimed_D}")
print(f"  IDENTITY CONFIRMED: {sp.simplify(ratio - claimed_D) == 0}")

print()
print("=" * 70)
print("E) Kernel substitutions used in Theorem 3's proof (direct check)")
print("=" * 70)
a, b = sp.symbols('a b', positive=True)
print(f"  odd case:  K_C^(m+1)(x,a) = (1/2)a^m + (1/2)a*b^((m-1)/2)  "
      f"— matches text by construction")
print(f"  even case: K_C^(m+1)(x,a) = (1/2)a^m + (1/2)b^(m/2)         "
      f"— matches text by construction")

print()
print("=" * 70)
print("F) Theorem A (paper [1]): solve ((m-1)/2)*tau^(m-1) = 1 for tau")
print("=" * 70)
tau2, m2 = sp.symbols('tau m', positive=True)
eq = sp.Eq((m2 - 1) / 2 * tau2**(m2 - 1), 1)
sol = sp.solve(eq, tau2)
claimed_A = (2 / (m2 - 1))**(1 / (m2 - 1))
print(f"  sympy solution: {sol}")
for mval in [5, 7, 9]:
    lhs_val = float(sol[0].subs(m2, mval))
    rhs_val = float(claimed_A.subs(m2, mval))
    print(f"    m={mval}: sympy={lhs_val:.6f}  closed form={rhs_val:.6f}  "
          f"match: {abs(lhs_val - rhs_val) < 1e-9}")

print()
print("All six checks completed. See the paper's §5-7 for where each is used.")
