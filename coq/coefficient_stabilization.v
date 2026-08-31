(*
  coefficient_stabilization.v

  A finite-order, fully computational confirmation, for m = 5 truncated
  at degree 13, that the Picard iteration defining the rooted generating
  function s_5(x) stabilizes coefficient by coefficient after a bounded
  number of iterations - and that the stabilized result matches the
  coefficients already published in [1] exactly: 1, 1, 3, 13 at degrees
  1, 5, 9, 13.

  Implemented as exact rational-number arithmetic (QArith), with each
  operation explicitly reduced via Qred - an earlier version without
  this reduction overflowed the stack by the second iteration, as
  unreduced fractions grew without bound (a real implementation issue,
  not a mathematical one).

  Scope, stated plainly: this confirms the stabilization phenomenon for
  one specific truncation order and one specific m. It does not
  establish that this phenomenon holds at every order and for every m -
  that more general fact is established, in the paper, by citing
  Pivoteau, Salvy and Soria's extension of Joyal's implicit species
  theorem, not reproved here.

  Compiles with no axioms at all ("Closed under the global context") for
  both theorems below: purely decidable computation.

  Note (FR) : confirmation entierement calculatoire, pour m=5 tronque au
  degre 13, que l'iteration de Picard definissant s_5(x) se stabilise
  coefficient par coefficient apres un nombre borne d'iterations, et que
  le resultat concorde exactement avec les coefficients deja publies
  (1, 1, 3, 13 aux degres 1, 5, 9, 13). Arithmetique rationnelle exacte
  (QArith). Compile sans aucun axiome.

  Author: Frederic G. Speyser
  Compile with: coqc coefficient_stabilization.v   (Coq 8.18 or later)
*)

Require Import QArith.
Require Import List.
Require Import Lia.
Import ListNotations.
Open Scope Q_scope.

Definition mul_trunc (N : nat) (a b : list Q) : list Q :=
  map (fun n =>
    Qred (fold_left (fun acc k => acc + (nth k a 0) * (nth (n-k) b 0)) (seq 0 (n+1)) 0)
  ) (seq 0 (N+1)).

Definition stretch (N r : nat) (a : list Q) : list Q :=
  map (fun n => if andb (Nat.eqb (n mod r) 0) (Nat.leb (n/r) N)
                then nth (n/r) a 0 else 0)
      (seq 0 (N+1)).

Definition add_series (N : nat) (a b : list Q) : list Q :=
  map (fun n => Qred (nth n a 0 + nth n b 0)) (seq 0 (N+1)).

Definition scale_series (N : nat) (c : Q) (a : list Q) : list Q :=
  map (fun n => Qred (c * nth n a 0)) (seq 0 (N+1)).

Fixpoint pow_series (N : nat) (a : list Q) (k : nat) : list Q :=
  match k with
  | O => map (fun n => if Nat.eqb n 0 then 1 else 0) (seq 0 (N+1))
  | S k' => mul_trunc N a (pow_series N a k')
  end.

Fixpoint exp_series (N : nat) (u : list Q) : list Q :=
  let fix build (n : nat) : list Q :=
    match n with
    | O => [1]
    | S n' =>
      let prev := build n' in
      let s := fold_left (fun a k =>
        a + (Z.of_nat k # 1) * (nth k u 0) * (nth (n - k) prev 0)
      ) (seq 1 n) 0 in
      prev ++ [Qred (s / (Z.of_nat n # 1))]
    end
  in build N.

Definition K_C5 (N : nat) (y s_x2 : list Q) : list Q :=
  add_series N
    (scale_series N (1#2) (pow_series N y 4))
    (scale_series N (1#2) (pow_series N s_x2 2)).

Definition zero_series (N : nat) : list Q := map (fun _ => 0) (seq 0 (N+1)).

Fixpoint full_exponent' (N : nat) (s : list Q) (imax : nat) : list Q :=
  match imax with
  | O => zero_series N
  | S i' =>
    let i := imax in
    let s_xi := stretch N i s in
    let s_x2i := stretch N (2*i) s in
    let term := scale_series N (1 # Pos.of_nat i) (K_C5 N s_xi s_x2i) in
    add_series N term (full_exponent' N s i')
  end.

Definition Phi5 (N : nat) (s : list Q) : list Q :=
  let x := map (fun n => if Nat.eqb n 1 then 1 else 0) (seq 0 (N+1)) in
  mul_trunc N x (exp_series N (full_exponent' N s N)).

Fixpoint picard (N : nat) (k : nat) : list Q :=
  match k with
  | O => map (fun n => if Nat.eqb n 1 then 1 else 0) (seq 0 (N+1))
  | S k' => Phi5 N (picard N k')
  end.

Eval vm_compute in (picard 13 1).
Eval vm_compute in (picard 13 2).
Eval vm_compute in (picard 13 3).
Eval vm_compute in (picard 13 4).
Eval vm_compute in (picard 13 7).
Eval vm_compute in (picard 13 13).

(* Stated as a genuine theorem, not merely a displayed computation:
   the iteration has stabilized from k=3 through k=13, coefficient by
   coefficient, for the truncation at order 13. *)
Theorem picard_stabilizes_by_3 :
  picard 13 3 = picard 13 13.
Proof. vm_compute. reflexivity. Qed.

Theorem picard_matches_known_coefficients :
  picard 13 13 = [0;1;0;0;0;1;0;0;0;3;0;0;0;13].
Proof. vm_compute. reflexivity. Qed.

Print Assumptions picard_stabilizes_by_3.
Print Assumptions picard_matches_known_coefficients.
