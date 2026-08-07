(*
  involution_counting.v

  The abstract combinatorial principle behind Theorem 2 (Section 5 of the
  paper): under an involution f on a finite set (f(f(x)) = x for every x),
  the number of orbits {x, f(x)} is at least half the size of the set.
  This is the exact mechanism behind the "at most 2-to-1" argument used
  for the path-reversal symmetry of chain cacti, stripped of any graph-
  theoretic content.

  Generalized to any type with decidable equality (not just nat), then
  instantiated concretely: chain configurations are represented as
  sequences over a finite alphabet (all_words), path-reversal as list
  reversal (rev, whose involutivity is already in Coq's standard
  library), and the concrete case s=2, L=3 corresponds to m=5, k=5
  blocks - giving at least 2^(k-3) = 4 pairwise distinct configurations
  out of 8, matching Theorem 2's formula exactly.

  Scope, stated plainly: this formalizes the abstract counting mechanism,
  not that path-reversal is genuinely an involution on actual cactus
  split-decomposition trees - establishing that would require formalizing
  graph theory and split-decomposition itself, not attempted here.

  Compiles with no axioms at all ("Closed under the global context"): a
  purely constructive proof, not even the classical axioms used in the
  two Reals-based files in this repository.

  Author: Frederic G. Speyser
  Compile with: coqc involution_counting.v   (Coq 8.18 or later)
*)

Require Import List.
Require Import PeanoNat.
Require Import Lia.
Import ListNotations.

Section InvolutionCountingGeneric.

(* Generalization of the counting principle to an arbitrary type with
   decidable equality, rather than specifically nat. Same proof structure
   as the nat-only version, adapted. *)

Variable A : Type.
Variable eqA_dec : forall x y : A, {x = y} + {x <> y}.
Variable elems : list A.
Hypothesis elems_NoDup : NoDup elems.

Variable f : A -> A.
Hypothesis Hbound : forall x, In x elems -> In (f x) elems.
Hypothesis Hinvol : forall x, In x elems -> f (f x) = x.

(* An arbitrary total order on A, via its index in elems, used to define
   "representative" as before (i <= f i becomes index(x) <= index(f x)). *)
Fixpoint index_of (x : A) (l : list A) : nat :=
  match l with
  | [] => 0
  | y :: l' => if eqA_dec x y then 0 else S (index_of x l')
  end.

Definition is_rep (x : A) : bool := Nat.leb (index_of x elems) (index_of (f x) elems).

Definition reps : list A := filter is_rep elems.
Definition nonreps : list A := filter (fun x => negb (is_rep x)) elems.

Lemma partition_length : length reps + length nonreps = length elems.
Proof.
  unfold reps, nonreps. apply filter_length.
Qed.

Lemma index_injective_on_list :
  forall (l : list A) (x y : A), In x l -> In y l -> index_of x l = index_of y l -> x = y.
Proof.
  induction l as [|a l' IH]; intros x y Hx Hy Heq.
  - destruct Hx.
  - simpl in Heq. destruct (eqA_dec x a) as [Hxa|Hxa]; destruct (eqA_dec y a) as [Hya|Hya].
    + subst. reflexivity.
    + discriminate.
    + discriminate.
    + injection Heq as Heq'.
      destruct Hx as [Hx|Hx]; [congruence|].
      destruct Hy as [Hy|Hy]; [congruence|].
      exact (IH x y Hx Hy Heq').
Qed.

Lemma index_injective_on_elems :
  forall x y, In x elems -> In y elems -> index_of x elems = index_of y elems -> x = y.
Proof. apply index_injective_on_list. Qed.

Lemma f_injective_on_elems :
  forall x y, In x elems -> In y elems -> f x = f y -> x = y.
Proof.
  intros x y Hx Hy Heq.
  assert (H1 : f (f x) = x) by (apply Hinvol; exact Hx).
  assert (H2 : f (f y) = y) by (apply Hinvol; exact Hy).
  rewrite Heq in H1. rewrite H1 in H2. exact H2.
Qed.

Lemma nonrep_image_is_rep :
  forall x, In x elems -> is_rep x = false -> is_rep (f x) = true.
Proof.
  intros x Hin Hnrep.
  unfold is_rep in Hnrep. apply Nat.leb_gt in Hnrep.
  assert (Hfx_in : In (f x) elems) by (apply Hbound; exact Hin).
  assert (Hffx : f (f x) = x) by (apply Hinvol; exact Hin).
  unfold is_rep. apply Nat.leb_le.
  rewrite Hffx. lia.
Qed.

Lemma nonreps_NoDup : NoDup nonreps.
Proof. unfold nonreps. apply NoDup_filter. exact elems_NoDup. Qed.

Lemma nonreps_elements_in_elems : forall x, In x nonreps -> In x elems.
Proof.
  intros x Hx. unfold nonreps in Hx. apply filter_In in Hx. destruct Hx as [Hx _]. exact Hx.
Qed.

Lemma map_injective_preserves_NoDup :
  forall (l : list A), NoDup l -> (forall x, In x l -> In x elems) ->
  NoDup (map f l).
Proof.
  induction l as [|a l IH]; intros Hnd Hbnd.
  - simpl. constructor.
  - simpl. inversion Hnd as [|a0 l0 Ha Hl0]; subst.
    constructor.
    + intro Hin. apply in_map_iff in Hin. destruct Hin as [x [Heqx Hxin]].
      assert (Hx_el : In x elems) by (apply Hbnd; right; exact Hxin).
      assert (Ha_el : In a elems) by (apply Hbnd; left; reflexivity).
      assert (Hai : a = x) by (apply f_injective_on_elems; [exact Ha_el|exact Hx_el|symmetry; exact Heqx]).
      apply Ha. rewrite Hai. exact Hxin.
    + apply IH.
      * exact Hl0.
      * intros y Hy. apply Hbnd. right. exact Hy.
Qed.

Lemma map_f_nonreps_NoDup : NoDup (map f nonreps).
Proof.
  apply map_injective_preserves_NoDup.
  - exact nonreps_NoDup.
  - exact nonreps_elements_in_elems.
Qed.

Lemma map_f_nonreps_incl_reps : incl (map f nonreps) reps.
Proof.
  intros y Hy. apply in_map_iff in Hy. destruct Hy as [x [Heq Hx_in]].
  unfold nonreps in Hx_in. apply filter_In in Hx_in. destruct Hx_in as [Hx_el Hx_nrep].
  apply Bool.negb_true_iff in Hx_nrep.
  assert (Hrep_fx := nonrep_image_is_rep x Hx_el Hx_nrep).
  rewrite Heq in Hrep_fx.
  unfold reps. apply filter_In. split.
  - rewrite <- Heq. apply Hbound. exact Hx_el.
  - exact Hrep_fx.
Qed.

Lemma nonreps_length_le_reps_length : length nonreps <= length reps.
Proof.
  rewrite <- (map_length f nonreps).
  apply NoDup_incl_length.
  - exact map_f_nonreps_NoDup.
  - exact map_f_nonreps_incl_reps.
Qed.

Theorem at_least_half_orbits_generic : 2 * length reps >= length elems.
Proof.
  pose proof partition_length as Hpart.
  pose proof nonreps_length_le_reps_length as Hle.
  lia.
Qed.

End InvolutionCountingGeneric.

Print Assumptions at_least_half_orbits_generic.

(* Concrete instantiation: configurations of length L over an alphabet of
   size s, representing k-2 = L interior blocks with s local classes each. *)

Fixpoint all_words (s L : nat) : list (list nat) :=
  match L with
  | O => [[]]
  | S L' =>
    flat_map (fun prefix => map (fun d => d :: prefix) (seq 0 s)) (all_words s L')
  end.

Eval vm_compute in (all_words 2 3).
Eval vm_compute in (length (all_words 2 3)).

Definition eqlist_dec : forall x y : list nat, {x=y}+{x<>y} := list_eq_dec Nat.eq_dec.

Lemma all_words_charac :
  forall s L w, In w (all_words s L) <-> (length w = L /\ Forall (fun d => d < s) w).
Proof.
  intros s L. induction L as [|L' IH]; intro w; simpl.
  - split.
    + intros [Heq|[]]. subst. split; [reflexivity | constructor].
    + intros [Hlen Hforall]. left. destruct w; [reflexivity | simpl in Hlen; discriminate].
  - split.
    + intro Hin. apply in_flat_map in Hin. destruct Hin as [prefix [Hpre Hmap]].
      apply in_map_iff in Hmap. destruct Hmap as [d [Heqw Hd_in]]. subst w.
      apply in_seq in Hd_in. destruct Hd_in as [_ Hd_lt].
      apply IH in Hpre. destruct Hpre as [Hlen Hforall].
      split.
      * simpl. rewrite Hlen. reflexivity.
      * constructor; [lia | exact Hforall].
    + intros [Hlen Hforall].
      destruct w as [|d prefix]; [simpl in Hlen; discriminate|].
      simpl in Hlen. injection Hlen as Hlen'.
      assert (Hd_lt : d < s) by (inversion Hforall; assumption).
      assert (Hforall' : Forall (fun x => x < s) prefix) by (inversion Hforall; assumption).
      apply in_flat_map. exists prefix. split.
      * apply IH. split; [exact Hlen' | exact Hforall'].
      * apply in_map_iff. exists d. split; [reflexivity | apply in_seq; lia].
Qed.

Lemma rev_closed_on_all_words :
  forall s L w, In w (all_words s L) -> In (rev w) (all_words s L).
Proof.
  intros s L w Hin. apply all_words_charac in Hin. destruct Hin as [Hlen Hforall].
  apply all_words_charac. split.
  - rewrite rev_length. exact Hlen.
  - apply Forall_forall. intros d Hd_in_rev.
    assert (Hd_in_w : In d w) by (apply in_rev; exact Hd_in_rev).
    rewrite Forall_forall in Hforall.
    exact (Hforall d Hd_in_w).
Qed.

Lemma all_words_2_3_NoDup : NoDup (all_words 2 3).
Proof.
  assert (Heq : nodup eqlist_dec (all_words 2 3) = all_words 2 3)
    by (vm_compute; reflexivity).
  rewrite <- Heq. apply NoDup_nodup.
Qed.

(* Final assembly: concrete instantiation for s=2, L=3 (k-2=3, i.e. k=5
   blocks, corresponding to m=5 since ceil((5-1)/2)=2 local classes). *)

Theorem chain_counting_m5_k5 :
  2 * length (reps (list nat) eqlist_dec (all_words 2 3) (@rev nat))
  >= length (all_words 2 3).
Proof.
  apply at_least_half_orbits_generic.
  - exact all_words_2_3_NoDup.
  - exact (rev_closed_on_all_words 2 3).
  - intros w _. apply rev_involutive.
Qed.

Print Assumptions chain_counting_m5_k5.
Eval vm_compute in (length (all_words 2 3)).
