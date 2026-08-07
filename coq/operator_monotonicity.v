(*
  operator_monotonicity.v

  Lemma 2 (Section 4 of the paper) in full, including its true infinite
  series - not merely its elementary building blocks. Given two families
  of nonnegative reals (u_i) and (v_i) with u_i <= v_i for every i, the
  operator built from u is bounded above, pointwise, by the operator
  built from v, where the operator is rho * exp(K(y, w_2) + sum_i K(w_i,
  w_2i)/i) for an abstract kernel K assumed only to be nonnegative and
  nondecreasing in each argument separately - exactly the properties the
  paper's informal proof actually invokes, nothing specific to the
  cactus kernel.

  The infinite series is handled via the Coquelicot library's Series_le,
  comparing two series termwise under a convergence hypothesis on the
  larger one.

  Compiles with no axioms beyond those already standard in Coq's Reals
  library (classical logic, functional extensionality) - see the
  Print Assumptions output at the end of the file. No axiom specific to
  Coquelicot is introduced.

  Author: Frederic G. Speyser
  Compile with: coqc operator_monotonicity.v
    (Coq 8.18 or later, with the Coquelicot library installed)
*)

Require Import Reals.
Require Import Coquelicot.Coquelicot.
Require Import Lra.
Require Import Lia.
Open Scope R_scope.

Section Lemma2FullSeries.

(* Noyau abstrait K, avec les proprietes de monotonie effectivement utilisees
   dans la demonstration informelle de l'article (composition de monotonies
   elementaires : positivite et croissance en chacun des deux arguments). *)
Variable K : R -> R -> R.
Hypothesis K_nonneg : forall a b, 0 <= a -> 0 <= b -> 0 <= K a b.
Hypothesis K_mono1 : forall a a' b, 0 <= a -> a <= a' -> 0 <= b -> K a b <= K a' b.
Hypothesis K_mono2 : forall a b b', 0 <= a -> 0 <= b -> b <= b' -> K a b <= K a b'.

Variable u v : nat -> R.
Hypothesis u_nonneg : forall n, 0 <= u n.
Hypothesis u_le_v : forall n, u n <= v n.

Lemma v_nonneg : forall n, 0 <= v n.
Proof. intro n. apply Rle_trans with (u n). apply u_nonneg. apply u_le_v. Qed.

(* Terme de la queue de la section 4/7 : K(w_i, w_2i)/i, indices decales de 2
   (n=0 correspond a i=2, etc.) pour tomber dans l'indexation naturelle de
   Coquelicot (series sur nat a partir de 0). *)
Definition term (w : nat -> R) (n : nat) : R :=
  K (w (n+2)%nat) (w (2*(n+2))%nat) / INR (n+2).

Lemma INR_n2_pos : forall n, 0 < INR (n+2).
Proof. intro n. apply lt_0_INR. lia. Qed.

Lemma term_nonneg_gen : forall (w : nat -> R), (forall n, 0 <= w n) -> forall n, 0 <= term w n.
Proof.
  intros w Hw n. unfold term.
  apply Rdiv_le_0_compat.
  - apply K_nonneg; apply Hw.
  - apply INR_n2_pos.
Qed.

Lemma term_le : forall n, term u n <= term v n.
Proof.
  intro n. unfold term.
  unfold Rdiv.
  apply Rmult_le_compat_r.
  - left. apply Rinv_0_lt_compat. apply INR_n2_pos.
  - apply Rle_trans with (K (v (n+2)%nat) (u (2*(n+2))%nat)).
    + apply K_mono1; [apply u_nonneg | apply u_le_v | apply u_nonneg].
    + apply K_mono2; [apply v_nonneg | apply u_nonneg | apply u_le_v].
Qed.

Hypothesis v_series_ex : ex_series (term v).

Theorem tail_series_le : Series (term u) <= Series (term v).
Proof.
  apply Series_le.
  - intro n. split.
    + apply term_nonneg_gen. exact u_nonneg.
    + apply term_le.
  - exact v_series_ex.
Qed.

(* Assemblage final : l'operateur complet, avec le terme K(y,w_2) et
   l'exponentielle, exactement comme dans l'enonce du Lemme 2. *)
Variable y rho : R.
Hypothesis y_nonneg : 0 <= y.
Hypothesis rho_pos : 0 < rho.

Definition full_operator (w : nat -> R) : R :=
  rho * exp (K y (w 2%nat) + Series (term w)).

Theorem operator_monotone : full_operator u <= full_operator v.
Proof.
  unfold full_operator.
  apply Rmult_le_compat_l.
  - lra.
  - assert (Hstep1 : K y (u 2%nat) <= K y (v 2%nat)).
    { apply K_mono2; [exact y_nonneg | apply u_nonneg | apply u_le_v]. }
    assert (Hsum_le : K y (u 2%nat) + Series (term u) <= K y (v 2%nat) + Series (term v)).
    { apply Rplus_le_compat; [exact Hstep1 | exact tail_series_le]. }
    destruct (Rle_lt_or_eq_dec _ _ Hsum_le) as [Hlt | Heq].
    + left. apply exp_increasing. exact Hlt.
    + right. rewrite Heq. reflexivity.
Qed.

End Lemma2FullSeries.

Print Assumptions operator_monotone.
