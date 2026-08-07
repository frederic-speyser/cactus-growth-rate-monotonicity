"""
mgonal_cactus_series.py

Computes the rooted and unrooted enumeration series for strict m-gonal cactus
graphs (free / non-plane case), m = 5, 6, 7, 8, from the defining functional
equations of the split-decomposition grammar (rooted fixed-point equation and
dissymmetry theorem), by exact rational formal power series arithmetic.

Reference: F. G. Speyser, "Enumeration of Non-Plane m-Gonal Cactus Graphs via
Split-Decomposition", Sections 5.1-5.3.

Author: Frederic G. Speyser
Run: python3 mgonal_cactus_series.py
"""
from fractions import Fraction as F

N = 26  # truncation order (>= highest term needed, m=8 unrooted needs x^22)


# ---------- formal power series utilities (lists of Fraction, index = degree) ----------

def zero():
    return [F(0)] * (N + 1)


def mul(a, b):
    c = zero()
    for i, ai in enumerate(a):
        if ai == 0:
            continue
        maxj = N - i
        if maxj < 0:
            continue
        for j, bj in enumerate(b[:maxj + 1]):
            if bj == 0:
                continue
            c[i + j] += ai * bj
    return c


def add(a, b):
    return [x + y for x, y in zip(a, b)]


def sub(a, b):
    return [x - y for x, y in zip(a, b)]


def scale(a, k):
    return [x * k for x in a]


def stretch(a, r):
    """coefficients of a(x^r)"""
    c = zero()
    for n, an in enumerate(a):
        if n * r <= N:
            c[n * r] = an
    return c


def shift_by_x(a):
    """multiply by x"""
    c = zero()
    for n in range(N):
        c[n + 1] = a[n]
    return c


def power_int(a, k):
    r = [F(0)] * (N + 1)
    r[0] = F(1)
    base = a
    while k > 0:
        if k & 1:
            r = mul(r, base)
        base = mul(base, base)
        k >>= 1
    return r


def exp_series(u):
    """exp(u) for a series u with u[0] == 0, via v' = u' v (Cauchy product recurrence)."""
    assert u[0] == 0
    v = zero()
    v[0] = F(1)
    # (n) v_n = sum_{k=1}^{n} k*u_k*v_{n-k}
    for n in range(1, N + 1):
        s = F(0)
        for k in range(1, n + 1):
            if u[k] != 0:
                s += k * u[k] * v[n - k]
        v[n] = s / n
    return v


# ---------- combinatorial specification (Section 5.1-5.3) ----------

def K_C(s, m):
    """K_C(x) from the rooted-grammar equation, for given rooted series s(x)."""
    s2 = stretch(s, 2)
    if m % 2 == 1:
        term1 = power_int(s, m - 1)
        term2 = power_int(s2, (m - 1) // 2)
        return scale(add(term1, term2), F(1, 2))
    else:
        term1 = power_int(s, m - 1)
        term2 = mul(s, power_int(s2, (m - 2) // 2))
        return scale(add(term1, term2), F(1, 2))


def sum_i_KC_xi_over_i(s, m):
    """Sigma_{i>=1} K_C(x^i)/i, truncated: K_C(x^i) has min degree i*(m-1)."""
    total = zero()
    i = 1
    while i * (m - 1) <= N:
        s_xi = stretch(s, i)
        kc_i = K_C(s_xi, m)          # K_C evaluated with s(x) -> s(x^i)
        kc_i_stretched_correctly = stretch(kc_i, 1)  # already in terms of x (built from s(x^i)), no further stretch needed
        total = add(total, scale(kc_i_stretched_correctly, F(1, i)))
        i += 1
    return total


def solve_s(m, iters=None):
    """Fixed point s = x + x*(exp(Sigma K_C(x^i)/i) - 1)."""
    if iters is None:
        iters = N + 2
    s = zero()
    s[1] = F(1)  # start from s = x
    for _ in range(iters):
        E = exp_series(sum_i_KC_xi_over_i(s, m))
        Eminus1 = list(E)
        Eminus1[0] -= 1
        s_new = zero()
        s_new[1] += 1  # the leading "x"
        s_new = add(s_new, shift_by_x(Eminus1))
        s = s_new
    return s


def phi_totient(n):
    result = n
    p = 2
    nn = n
    while p * p <= nn:
        if nn % p == 0:
            while nn % p == 0:
                nn //= p
            result -= result // p
        p += 1
    if nn > 1:
        result -= result // nn
    return result


def divisors(n):
    return [d for d in range(1, n + 1) if n % d == 0]


def Z_Dm(s, m):
    """Z_{D_m}(s(x), s(x^2), ..., s(x^m)) — dihedral cycle index, per the general
    formula stated in Section 5.3."""
    p = {i: stretch(s, i) for i in range(1, m + 1)}
    total = zero()
    for d in divisors(m):
        term = power_int(p[d], m // d)
        total = add(total, scale(term, F(phi_totient(d), 2 * m)))
    if m % 2 == 1:
        extra = mul(p[1], power_int(p[2], (m - 1) // 2))
        total = add(total, scale(extra, F(1, 2)))
    else:
        extra1 = mul(power_int(p[1], 2), power_int(p[2], (m - 2) // 2))
        extra2 = power_int(p[2], m // 2)
        total = add(total, scale(add(extra1, extra2), F(1, 4)))
    return total


def solve_G(m):
    s = solve_s(m)
    KC = K_C(s, m)
    E = exp_series(sum_i_KC_xi_over_i(s, m))
    Eminus1 = list(E)
    Eminus1[0] -= 1
    S_X = shift_by_x(Eminus1)                 # S_X(x) = x*(E(x)-1)
    S_C = sub(sub(E, [F(1)] + [F(0)] * N), KC)  # S_C(x) = E(x) - 1 - K_C(x)
    T_S = shift_by_x(S_C)                      # T_S(x) = x * S_C(x)
    T_SCm = mul(KC, S_X)                        # T_{S-Cm}(x) = K_C(x) * S_X(x)
    T_Cm = Z_Dm(s, m)
    G = add(sub(T_Cm, T_SCm), T_S)
    return s, G


def series_str(coeffs, maxdeg=None):
    parts = []
    for n, c in enumerate(coeffs):
        if maxdeg is not None and n > maxdeg:
            break
        if c == 0:
            continue
        assert c.denominator == 1, f"non-integer coefficient at x^{n}: {c}"
        ci = c.numerator
        if n == 0:
            parts.append(f"{ci}")
        elif n == 1:
            parts.append(f"{ci}x" if ci != 1 else "x")
        else:
            parts.append(f"{ci}x^{n}" if ci != 1 else f"x^{n}")
    return " + ".join(parts)


# ---------- published values, Section 6.1-6.2, for consistency display ----------

published_rooted = {
    5: "x + x^5 + 3x^9 + 13x^13 + 62x^17 + 333x^21",
    6: "x + x^6 + 4x^11 + 22x^16 + 140x^21",
    7: "x + x^7 + 4x^13 + 25x^19",
    8: "x + x^8 + 5x^15",
}
published_unrooted = {
    5: "x^5 + x^9 + 3x^13 + 8x^17 + 31x^21",
    6: "x^6 + x^11 + 4x^16 + 13x^21",
    7: "x^7 + x^13 + 4x^19",
    8: "x^8 + x^15 + 5x^22",
}

if __name__ == "__main__":
    print(f"Truncation order N = {N}\n")
    all_ok = True
    for m in [5, 6, 7, 8]:
        s, G = solve_G(m)
        maxdeg_rooted = max(int(t.split('^')[1]) if '^' in t else 1
                             for t in published_rooted[m].split(' + '))
        maxdeg_unrooted = max(int(t.split('^')[1]) if '^' in t else 1
                               for t in published_unrooted[m].split(' + '))
        got_rooted = series_str(s, maxdeg_rooted)
        got_unrooted = series_str(G, maxdeg_unrooted)
        ok_r = got_rooted.replace(" ", "") == published_rooted[m].replace(" ", "")
        ok_u = got_unrooted.replace(" ", "") == published_unrooted[m].replace(" ", "")
        all_ok = all_ok and ok_r and ok_u
        print(f"m = {m}")
        print(f"  rooted   (published): {published_rooted[m]}")
        print(f"  rooted   (computed) : {got_rooted}   {'OK' if ok_r else 'MISMATCH'}")
        print(f"  unrooted (published): {published_unrooted[m]}")
        print(f"  unrooted (computed) : {got_unrooted}   {'OK' if ok_u else 'MISMATCH'}")
        print()
    print("ALL SERIES MATCH" if all_ok else "DISCREPANCY FOUND")

    # ---------- OEIS submission data: first 25 terms of each sequence ----------
    # Reproduces exactly the DATA fields listed in oeis-submissions.html.
    print("\nOEIS submission data (25 terms per sequence)\n")
    N_OEIS_TERMS = 25
    N_saved = N
    for m in [5, 6, 7, 8]:
        N = max(N_saved, 180)  # ensure enough truncation order to reach 25 nonzero terms, even for m=8
        s, G = solve_G(m)
        rooted_terms = [c.numerator for n, c in enumerate(s) if c != 0][:N_OEIS_TERMS]
        unrooted_terms = [c.numerator for n, c in enumerate(G) if c != 0][:N_OEIS_TERMS]
        print(f"m = {m}")
        print(f"  rooted   (offset 0): {', '.join(map(str, rooted_terms))}")
        print(f"  unrooted (offset 1): {', '.join(map(str, unrooted_terms))}")
        print()
    N = N_saved
