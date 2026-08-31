# Strict Monotonicity and a Lambert-W Asymptotic for Growth Rates of Non-Plane Strict m-Gonal Cacti - verification code

This repository accompanies the paper that proves the growth-rate monotonicity conjecture originally raised for strict non-plane m-gonal cacti, and goes substantially further: it establishes a precise large-m asymptotic for the radius of convergence itself, with a closed-form rate involving the Lambert W function.

A series of working papers accompanied the consolidation and preparation of this result. The working paper itself has since been consolidated into a single, self-contained work, and this repository has been updated to match it exactly: every script here verifies a claim that appears, under the same theorem/table numbers, in the current preprint, titled "Strict Monotonicity and a Lambert-W Asymptotic for Growth Rates of Non-Plane Strict m-Gonal Cacti."


See the **[COMPANION PAGE](https://frederic-speyser.github.io/cactus-growth-rate-monotonicity/)** for an overview of the results and a verification matrix linking each claim to the script that checks it (includes a summary in French).

Specifically, the paper proves:

- **Strict monotonicity.** The radius of convergence rho_m of the rooted enumeration series is strictly increasing in m - the exponential growth rate 1/rho_m is strictly decreasing. The proof combines a purely combinatorial lower bound, a non-circular monotonicity result for the sequence of rooted generating functions (via controlled Picard iterates of the defining functional equation), and a differential criterion from a convexity property of the underlying fixed-point operator.
- **A Lambert-W asymptotic.** As m tends to infinity,

      1 - rho_m = lambda_m - (1/2)lambda_m^2 + O((log m)/m^2),

  where lambda_m = (log((m-1)/2) + 1 + W(e^(-1))) / (m-1) and W is the Lambert function. Thus the exponential growth rate decreases to 1 at rate (log m)/m.

**Companion enumeration paper:** Speyser, "Enumeration and Asymptotic Analysis of Strict Non-Plane m-Gonal Cactus Graphs via Split-Decomposition" - code and data at [non-plane-mgonal-cacti](https://github.com/frederic-speyser/non-plane-mgonal-cacti). The present paper reuses that paper's rooted/unrooted series (Tables 1-2) as its starting point; Tables 3-4 and Figure 3 of the present paper are new and verified here.

## What this repository contains

### `paper/`

The LaTeX source and compiled PDF of the preprint itself.

### `python/`

| File | What it verifies |
|---|---|
| `critical_point.py` | Solves the critical system (rho_m, tau_m) directly, by damped Picard iteration on the (x, s(x^2), s(x^4), ...) tree. This is the accurate reference method: its output matches the paper's Table 3 to 5-6 significant figures, and is used as the trusted base for every other script below. |
| `compute_rho.py` | Independent estimate of rho_m via a ratio test on high-order series coefficients, corrected to include the n^(-3/2) polynomial factor from the transfer theorem (the original, uncorrected version had a systematic bias of order 10^-2; the corrected version agrees with `critical_point.py` to 10^-4). |
| `newton_rho.py` | A second, independent ratio-test estimate of rho_m, same n^(-3/2) correction applied, agreeing with `critical_point.py` to 10^-5. |
| `verify_transfer_constants.py` | Independent recomputation of Phi_y and Phi_hh (Theorem 3.3, Table 3) by direct finite differences on the reduced (y, h) critical system - a genuinely different route from the closed-form derivation in the paper's proof, not a restatement of it. Matches Table 3 to within 1-3%. |
| `extended_range_m30.py` | Extends Table 3/Table 4 and the check of Theorem 5.3's two-term equivalent from m=12 up to m=30, confirming numerically that the residual, rescaled by m^2/log(m), stabilises in [4.4, 5.2] across the whole range - the evidence, quoted in the paper, that the stated O((log m)/m^2) error term is the right order. |
| `verify_theorem_equiv_symbolic.py` | Independent symbolic (SymPy) verification of the algebraic chain in the proof of Theorem 5.3: that v*e^v=e^(-1) has solution v=W(e^(-1)), that the derivative of u-1-e^(-u) at 1+Omega equals 1+Omega exactly, and the Taylor expansion of 1-e^(-Lambda) used in the final transfer step. Checks the identities themselves, not numerical values. |

### `pari/`

| File | What it verifies |
|---|---|
| `verify_critical_pari.gp` | Independent recomputation of rho_m (Table 3) using native PARI/GP truncated power series arithmetic, via an Euler-transform recurrence for the block-indexed rooted series - a different language, a different underlying arithmetic engine, and a different indexing convention (block count rather than vertex count) from every Python script above. Matches Table 3 to within 3x10^-5. |

### `coq/`

Four Coq developments, each compiling without `Admitted`:

- **`convexity_lemma.v`** - Lemma 4.6 of the paper: a general analytic fact about strictly increasing derivatives, proved via the mean value theorem.
- **`operator_monotonicity.v`** - the operator-monotonicity step underlying Theorem 4.4, in full, including the genuine infinite series (using the Coquelicot library), not merely its elementary building blocks.
- **`involution_counting.v`** - an abstract combinatorial counting principle (under an involution on a finite set, the number of orbits is at least half its size), generalised to any type with decidable equality and instantiated concretely for m=5.
- **`coefficient_stabilization.v`** - a finite-order, fully computational confirmation (m=5, truncated at degree 13) that the Picard iteration defining the rooted series stabilises, matching the coefficients independently published and OEIS-verified for that case.

**On the role of this code.** The mathematical proof is complete in the paper itself, on paper, independently of any code. The Coq developments verify specific lemmas and an abstract principle underlying the argument - algebraic identities, numerical bounds, one counting principle - not the full inductive proof of Theorem 4.7 (strict monotonicity) in its entirety, nor the paper's other results. This is a supplementary layer of verification, in the same spirit as the independent numerical cross-checks in `python/`.

All four files were compiled and verified against Coq 8.18.0 with the Coquelicot library (3.4.1), with no custom axioms introduced beyond the classical axioms already standard throughout Coq's own `Reals` library (two of the four files use no axioms at all).

## Requirements

- Python 3, with `numpy` and `sympy`.
- PARI/GP (any recent version) for the script in `pari/`.
- Coq 8.18 or later to compile the files in `coq/`; the Coquelicot library is additionally required for `operator_monotonicity.v` only.

## Usage

```bash
cd python/
python3 critical_point.py
python3 compute_rho.py
python3 newton_rho.py
python3 verify_transfer_constants.py
python3 extended_range_m30.py
python3 verify_theorem_equiv_symbolic.py   # requires: pip install sympy
```

```bash
cd pari/
gp -q verify_critical_pari.gp   # requires PARI/GP
```

```bash
cd coq/
coqc convexity_lemma.v
coqc involution_counting.v
coqc coefficient_stabilization.v
coqc operator_monotonicity.v   # requires Coquelicot
```

## Citation

If you use this code, please cite the preprint:

Frédéric G. Speyser, "Strict Monotonicity and a Lambert-W Asymptotic for Growth Rates of Non-Plane Strict m-Gonal Cacti," preprint, 2026. DOI: https://doi.org/10.5281/zenodo.22206801

A citable archive of this repository (the code accompanying the preprint above) is available via Zenodo. DOI: https://doi.org/10.5281/zenodo.21840266

## Author

Frédéric G. Speyser - Independent Researcher, Association Sciences & Cooperation, France
ORCID: 0000-0002-1767-5325

## License

MIT (see `LICENSE`).
