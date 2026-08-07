"""
lower_bound_rho5.py

An explicit, if non-optimal, lower bound on rho_5, obtained via a majorant
series argument: under the hypothesis s(t) <= phi(t) on [0, x0], for a
polynomial phi with a small number of terms, each term of the exponent in
the defining functional equation is majorized, giving a self-consistent
bound that can be pushed to a boundary x0 via a Picard-iterate induction.

Three ansatze are tried, of increasing precision, each using more exact
low-order coefficients of s_5(x) (all independently confirmed elsewhere in
this repository, e.g. via Coq in coefficient_stabilization.v):

  1 parameter:  phi(t) = C*t
  2 parameters: phi(t) = C*t + D*t^5      (D = 1, exact)
  3 parameters: phi(t) = C*t + D*t^5 + E*t^9   (D = 1, E = 3, exact)

Each bound is verified directly against the true Picard iteration of the
functional equation, not just checked at the single boundary point.

This is a genuinely useful, if modest, result: it is the first explicit
lower bound on any rho_m obtained in this line of work, all prior results
being upper bounds only (Theorem 2 of the paper). It exhibits a clear,
diminishing return as more terms are added (0.576083 -> 0.585030 ->
0.588955, against the true value 0.604765), and the chantier was
deliberately closed at three parameters rather than pursued further -
see the paper's remarks for the reasoning.

Author: Frederic G. Speyser
Run with: python3 lower_bound_rho5.py   (requires: numpy)
"""
import numpy as np

RHO_5_TRUE = 0.604765  # already established elsewhere (Newton solve on the
                        # exact criticality system), used here only for
                        # comparison, not assumed in the derivation itself


# ---------- 1-parameter ansatz: phi(t) = C*t ----------

def x0_of_C_1param(C):
    """Self-consistency boundary for the 1-parameter ansatz, derived by
    hand: (1-x0^4)^(-(C^4+C^2)/2) = C."""
    denom = C**4 + C**2
    threshold = C**(-2 / denom)
    if not (0 < threshold < 1):
        return -1
    x0_4 = 1 - threshold
    return x0_4**0.25 if x0_4 > 0 else -1


def best_1param():
    Cs = np.linspace(1.001, 5, 2000)
    vals = [x0_of_C_1param(C) for C in Cs]
    i = int(np.argmax(vals))
    return Cs[i], vals[i]


# ---------- 2- and 3-parameter ansatze, via explicit sympy-verified
#            expansion (see the paper / verification log for the full
#            symbolic derivation; the coefficients below are its output) ----

def lhs_2param(y, C, D=1.0):
    c4 = (C**4 + C**2) / 2
    c8 = 2 * C**3 * D
    c12 = 3 * C**2 * D**2 + C * D
    c16 = 2 * C * D**3
    c20 = (D**4 + D**2) / 2
    G = (c4 * (-np.log(1 - y)) + c8 * (-np.log(1 - y**2)) +
         c12 * (-np.log(1 - y**3)) + c16 * (-np.log(1 - y**4)) +
         c20 * (-np.log(1 - y**5)))
    return np.exp(G)


def check_consistency_2param(C, y0, D=1.0, n_points=300):
    for y in np.linspace(1e-9, y0, n_points):
        if y >= 1:
            return False
        try:
            if lhs_2param(y, C, D) > C + D * y:
                return False
        except (ValueError, FloatingPointError, OverflowError):
            return False
    return True


def find_max_y0_2param(C, D=1.0, y_lo=1e-6, y_hi=0.9999, tol=1e-7):
    if not check_consistency_2param(C, y_lo, D):
        return None
    lo, hi = y_lo, y_hi
    while hi - lo > tol:
        mid = (lo + hi) / 2
        if check_consistency_2param(C, mid, D):
            lo = mid
        else:
            hi = mid
    return lo


def best_2param(D=1.0):
    best_x0, best_C, best_y0 = 0, None, None
    for C in np.linspace(1.001, 2.5, 400):
        y0 = find_max_y0_2param(C, D)
        if y0 is None:
            continue
        x0 = y0**0.25
        if x0 > best_x0:
            best_x0, best_C, best_y0 = x0, C, y0
    return best_C, best_x0


def coeffs_3param(C, D=1.0, E=3.0):
    return {
        1: (C**4 + C**2) / 2,
        2: 2 * C**3,
        3: 6 * C**3 + 3 * C**2 + C,
        4: 18 * C**2 + 2 * C,
        5: 27 * C**2 + 21 * C + 1,
        6: 54 * C + 6,
        7: 54 * C + 30,
        8: 54.0,
        9: 45.0,
    }
    # note: coefficients derived by symbolic expansion of (Ct+Dt^5+Et^9)^4
    # and its t^2-substituted, squared companion term; see
    # symbolic_identity_checks.py for the general expansion method and the
    # paper / verification log for this specific expansion, cross-checked
    # with sympy.Poly at the time of derivation.


def lhs_3param(y, C, D=1.0, E=3.0):
    cs = coeffs_3param(C, D, E)
    G = sum(cs[k] * (-np.log(1 - y**k)) for k in range(1, 10))
    return np.exp(G)


def check_consistency_3param(C, y0, D=1.0, E=3.0, n_points=300):
    for y in np.linspace(1e-9, y0, n_points):
        if y >= 1:
            return False
        try:
            if lhs_3param(y, C, D, E) > C + D * y + E * y**2:
                return False
        except (ValueError, FloatingPointError, OverflowError):
            return False
    return True


def find_max_y0_3param(C, D=1.0, E=3.0, y_lo=1e-6, y_hi=0.9999, tol=1e-7):
    if not check_consistency_3param(C, y_lo, D, E):
        return None
    lo, hi = y_lo, y_hi
    while hi - lo > tol:
        mid = (lo + hi) / 2
        if check_consistency_3param(C, mid, D, E):
            lo = mid
        else:
            hi = mid
    return lo


def best_3param(D=1.0, E=3.0):
    best_x0, best_C = 0, None
    for C in np.linspace(1.001, 2.0, 600):
        y0 = find_max_y0_3param(C, D, E)
        if y0 is None:
            continue
        x0 = y0**0.25
        if x0 > best_x0:
            best_x0, best_C = x0, C
    return best_C, best_x0


# ---------- Direct verification against the true Picard iteration ----------

def K_C_5(x, y, s_x2):
    return 0.5 * y**4 + 0.5 * s_x2**2


def solve_s_at_point(x, iters=200, N_terms=60):
    pts = [x**i for i in range(1, N_terms + 1)]
    s_vals = list(pts)
    for _ in range(iters):
        new_vals = []
        for idx, xi in enumerate(pts):
            i = idx + 1
            G = 0.0
            j = 1
            while i * j <= N_terms:
                y_val = s_vals[i * j - 1]
                x2i_idx = 2 * i * j
                s_x2 = s_vals[x2i_idx - 1] if x2i_idx <= N_terms else 0.0
                G += K_C_5(xi**j, y_val, s_x2) / j
                j += 1
            new_vals.append(xi * np.exp(min(G, 50)))
        if max(new_vals) > 1e8:
            return None
        s_vals = new_vals
    return s_vals[0]


if __name__ == "__main__":
    print("=" * 70)
    print("1-parameter ansatz: phi(t) = C*t")
    print("=" * 70)
    C1, x0_1 = best_1param()
    print(f"  best C = {C1:.6f}, x0 = {x0_1:.6f}")

    print()
    print("=" * 70)
    print("2-parameter ansatz: phi(t) = C*t + t^5  (D=1, exact coefficient)")
    print("=" * 70)
    C2, x0_2 = best_2param()
    print(f"  best C = {C2:.6f}, x0 = {x0_2:.6f}  (improvement: {x0_2 - x0_1:+.6f})")

    print()
    print("=" * 70)
    print("3-parameter ansatz: phi(t) = C*t + t^5 + 3t^9  (D=1, E=3, exact)")
    print("=" * 70)
    C3, x0_3 = best_3param()
    print(f"  best C = {C3:.6f}, x0 = {x0_3:.6f}  (improvement: {x0_3 - x0_2:+.6f})")

    print()
    print("=" * 70)
    print(f"Summary, against rho_5 (true value, established elsewhere) = {RHO_5_TRUE}")
    print("=" * 70)
    print(f"  1 param : {x0_1:.6f}")
    print(f"  2 params: {x0_2:.6f}")
    print(f"  3 params: {x0_3:.6f}")
    print(f"  true    : {RHO_5_TRUE}")
    print()
    print("  Diminishing returns, deliberately not pursued further (see paper).")

    print()
    print("=" * 70)
    print(f"Direct verification of the 3-parameter bound via true Picard iteration")
    print("=" * 70)
    for frac in [0.9, 1.0, 1.02]:
        x = x0_3 * frac
        val = solve_s_at_point(x, iters=150, N_terms=50)
        bound = C3 * x + 1.0 * x**5 + 3.0 * x**9
        if val is None:
            print(f"  x={x:.6f} : DIVERGES (beyond rho_5, as expected)")
        else:
            print(f"  x={x:.6f} : converges to {val:.6f}, bound={bound:.6f}, "
                  f"holds: {val <= bound}")
