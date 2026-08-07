# Verification code for "Strict Monotonicity of the Growth Rate for Non-Plane Strict *m*-Gonal Cactus Graphs"


## Rationale

In [2], I enumerate strict non-plane *m*-gonal cactus graphs for each fixed *m* ≥ 5, and conjecture (Conjecture 1) that the exponential growth rate 1/ρ*m* of this family is strictly decreasing in *m*, on strong numerical evidence but without a proof.

This repository accompanies the paper that proves it. The proof combines a purely combinatorial lower bound on the growth rate, a non-circular monotonicity result for the sequence of rooted generating functions obtained by controlling the Picard iterates of the defining functional equation, and a differential criterion obtained from a convexity property of the underlying fixed-point operator.

Three pieces of this proof have been independently verified by machine, using the Coq proof assistant: the convexity lemma, the operator monotonicity lemma (with its true infinite series, not a finite approximation), and the abstract combinatorial counting principle behind the paper's chain-cactus argument. One of these three — the counting principle — carries no axioms at all; the other two rely only on the classical axioms already standard throughout Coq's own Reals library, nothing introduced by this work itself.

**On the role of this code.** The mathematical proof is complete in the paper itself, on paper, independently of any code. The material in this repository consists of tests performed afterward, using independent methods — exact symbolic computation, and mechanized checking in the Coq proof assistant — to check specific claims made in the proof: algebraic identities, numerical bounds, and one abstract combinatorial principle. This code does not replace or complete the proof; it is a supplementary layer of verification, in the same spirit as the independent PARI/GP implementation and structural tests already used to cross-check the enumerative data of [2].

The preprint of [2] is available on Zenodo (DOI [10.5281/zenodo.21461100](https://doi.org/10.5281/zenodo.21461100), a link that will always point to its latest version, even after future updates), and the accompanying code on GitHub: [non-plane-mgonal-cacti](https://github.com/frederic-speyser/non-plane-mgonal-cacti). That paper has already been submitted to the Electronic Journal of Combinatorics. This paper — the one this repository provides verification code for — is still in preparation for submission to the same journal; its Zenodo archive will be linked in the References section below once deposited.

## Related repositories

- [**non-plane-mgonal-cacti**](https://github.com/frederic-speyser/non-plane-mgonal-cacti) — code and data for [2], the enumeration paper this work builds on and proves the main conjecture of.
- [**cactus-split-decomp-omega**](https://github.com/frederic-speyser/cactus-split-decomp-omega) — an independent, exploratory extension toward a mixed set Ω of admissible cycle lengths, rather than a single fixed *m*. Numerical and observational only, no theorem; unrelated to the proof in this repository, which is specific to the singleton case.

## What this repository contains

- **`coq/`** — four Coq developments, each compiling without `Admitted`:
  - `convexity_lemma.v` — Lemma 1 (§4): a general analytic fact about strictly increasing derivatives, proved via the mean value theorem.
  - `operator_monotonicity.v` — Lemma 2 (§4) in full, including the genuine infinite series (using the Coquelicot library), not merely its elementary building blocks.
  - `involution_counting.v` — the abstract counting principle behind Theorem 2 (§5): under an involution on a finite set, the number of orbits is at least half its size. Generalized to any type with decidable equality, and instantiated concretely for the case corresponding to *m* = 5.
  - `coefficient_stabilization.v` — a finite-order, fully computational confirmation (for *m* = 5, truncated at degree 13) that the Picard iteration defining the rooted series stabilizes, matching the coefficients already published in [2].
- **`python/`** — the scripts supporting the numerical and symbolic claims made in the paper, including the exact-rational series solver adapted from [2], the symbolic (sympy) verification of the algebraic identities of §5–7, and the derivation of an explicit lower bound on ρ₅.

## What this repository does not contain

No mechanized proof of Theorem 3 (§6) itself, nor of the paper's main theorem as a whole. The Coq developments above verify specific lemmas and an abstract principle underlying the argument, not the full inductive proof or its combinatorial content.

No formalization of the general convergence theorem for Picard iterates on the true functional equation (as opposed to a finite coefficient prefix): this rests on citing Pivoteau, Salvy and Soria's [5] extension of Joyal's implicit species theorem, whose hypotheses have been checked to apply to the present equation, but which is not reproved here.

No formalization connecting the abstract involution-counting principle to actual cactus graphs and their split-decomposition trees: this would require formalizing graph theory and split-decomposition itself, a substantially larger undertaking not attempted here.

## Requirements

- **Coq** 8.18 or later to compile the files in `coq/`; the **Coquelicot** library is additionally required for `operator_monotonicity.v` only.
- **Python** 3, with `sympy` and `numpy`, to run the scripts in `python/`.

## References

[1] Speyser, F. G. *Strict Monotonicity of the Growth Rate for Non-Plane Strict m-Gonal Cactus Graphs.* Preprint, currently being prepared for submission to the Electronic Journal of Combinatorics, 2026. Zenodo archive: to be added here once deposited.

[2] Speyser, F. G. *Enumeration and Asymptotic Analysis of Strict Non-Plane m-Gonal Cactus Graphs via Split-Decomposition.* Submitted to the Electronic Journal of Combinatorics, 2026. Preprint: [10.5281/zenodo.21461100](https://doi.org/10.5281/zenodo.21461100).

[3] Bahrani, M., Lumbroso, J. *Enumerations, Forbidden Subgraph Characterizations, and the Split-Decomposition.* Electronic Journal of Combinatorics 25(4), #P4.47, 2018. DOI: [10.37236/6431](https://doi.org/10.37236/6431).

[4] Bahrani, M., Lumbroso, J. *Split-Decomposition Trees with Prime Nodes: Enumeration and Random Generation of Cactus Graphs.* Proceedings of ANALCO 2018, pp. 143–157. DOI: [10.1137/1.9781611975062.13](https://doi.org/10.1137/1.9781611975062.13).

[5] Pivoteau, C., Salvy, B., Soria, M. *Algorithms for Combinatorial Structures: Well-Founded Systems and Newton Iterations.* Journal of Combinatorial Theory, Series A 119(8), 2012, pp. 1711–1773. DOI: [10.1016/j.jcta.2012.05.007](https://doi.org/10.1016/j.jcta.2012.05.007).

## Citation

If you use this code, please cite the paper listed at the top of the References section above. A citable archive of this repository will be added here via Zenodo once deposited.

## Author

Frédéric G. Speyser - Independent Researcher, France - ORCID: 0000-0002-1767-5325

## License

MIT (see `LICENSE`).
