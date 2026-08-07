(*
  convexity_lemma.v

  Lemma 1 (Section 4 of the paper): a general analytic fact, independent
  of any combinatorial content. Let g be differentiable on an interval,
  with g' strictly increasing, g(sigma) = 0 and g(tau) < 0 for some
  sigma < tau in that interval. Then g'(sigma) < 0.

  Proved via the mean value theorem (MVT_cor2 from Coq's standard
  Ranalysis library), by contradiction: if g'(sigma) >= 0, strict
  monotonicity of g' forces g(tau) > g(sigma) = 0, contradicting the
  hypothesis.

  Compiles with no axioms beyond those already standard in Coq's Reals
  library (classical logic, functional extensionality) - see the
  Print Assumptions output at the end of the file.

  Author: Frederic G. Speyser
  Compile with: coqc convexity_lemma.v   (Coq 8.18 or later)
*)

Require Import Reals.
Require Import Ranalysis.
Require Import Lra.
Open Scope R_scope.

Theorem convexity_lemma :
  forall (g g' : R -> R) (sigma tau : R),
  sigma < tau ->
  (forall x, sigma <= x <= tau -> derivable_pt_lim g x (g' x)) ->
  (forall x y, sigma <= x -> x < y -> y <= tau -> g' x < g' y) ->
  g sigma = 0 ->
  g tau < 0 ->
  g' sigma < 0.
Proof.
  intros g g' sigma tau Hlt Hderiv Hmono Hgsigma Hgtau.
  destruct (Rlt_dec (g' sigma) 0) as [Hcase | Hcase].
  - exact Hcase.
  - exfalso.
    apply Rnot_lt_le in Hcase.
    destruct (MVT_cor2 g g' sigma tau Hlt Hderiv) as [c [Heq Hrange]].
    destruct Hrange as [Hc1 Hc2].
    assert (Hgtau_eq : g tau = g' c * (tau - sigma)).
    { rewrite Hgsigma in Heq. lra. }
    assert (Hgc_neg : g' c < 0).
    { assert (Htau_sigma_pos : tau - sigma > 0) by lra.
      nra. }
    assert (Hgc_gt : g' sigma < g' c).
    { apply Hmono; lra. }
    lra.
Qed.

Print Assumptions convexity_lemma.
