/*
  verify_critical_pari.gp

  Independent recomputation of rho_m (Table 3 of the paper) by native
  PARI/GP truncated power series arithmetic, using an Euler-transform
  recurrence for the block-indexed rooted series a_k = s_{m,1+(m-1)k} -
  a genuinely different route from the two Python methods in this
  repository (critical_point.py: pointwise Picard on the (x, s(x^2),
  s(x^4), ...) tree; compute_rho.py / newton_rho.py: coefficient ratio
  test on a vertex-indexed series computed via numpy polynomial
  arithmetic). Here the series is built natively in PARI's Ser type,
  block-indexed rather than vertex-indexed, and rho_m is extracted via
  the corresponding n^(-3/2)-corrected ratio test, giving a cross-check
  in both a different language and a different underlying arithmetic
  engine, and a different indexing convention entirely.

  Note (FR) : recalcul independant de rho_m (Tableau 3 de l'article) par
  arithmetique de series tronquees native de PARI/GP, via une recurrence
  de type transformee d'Euler pour la serie enracinee indexee par nombre
  de blocs a_k = s_{m,1+(m-1)k} - une methode genuinement differente des
  deux approches Python de ce depot. La serie est construite nativement
  dans le type Ser de PARI, indexee par blocs plutot que par degre, et
  rho_m est extrait par le test de ratio corrige correspondant, donnant
  une verification croisee a la fois dans un autre langage, un autre
  moteur arithmetique, et une convention d'indexation entierement
  differente.

  Reference: F. G. Speyser, "Strict Monotonicity and a Lambert-W Asymptotic
  for Growth Rates of Non-Plane Strict m-Gonal Cacti", Table 3.

  Author: Frederic G. Speyser
  Run with: gp -q verify_critical_pari.gp
*/

EulerT(v) = {Vec(exp(x*Ser(dirmul(v, vector(#v, n, 1/n))))-1, -#v)};

\\ Block-indexed rooted series: a[k+1] = number of rooted m-gonal cacti
\\ with k blocks, k = 0, 1, ..., N (a[1] = 1, the trivial single vertex).
series_blocks(mm, N) = {
  my(v = []);
  for(n = 1, N,
    my(g = 1 + x*Ser(v));
    if(mm % 2 == 1,
      v = EulerT(Vec((g^(mm-1) + subst(g^((mm-1)\2), x, x^2)) / 2)),
      v = EulerT(Vec((g^(mm-1) + g*subst(g^((mm-2)\2), x, x^2)) / 2))
    );
  );
  concat([1], v);
};

\\ rho_m estimate via the n^(-3/2)-corrected ratio test, applied to the
\\ last two available block-indexed coefficients, translated back to
\\ vertex count n = 1 + (m-1)*k.
rho_from_blocks(a, mm) = {
  my(K = #a - 1, ak = a[K], akp1 = a[K+1]);
  my(nk = 1+(mm-1)*(K-1), nkp1 = 1+(mm-1)*K);
  my(corrected = (ak/akp1) * (nk/nkp1)^1.5);
  corrected^(1.0/(mm-1));
};

known5 = 0.604765; known6 = 0.633235; known7 = 0.669930; known8 = 0.690268;

print("  m    rho_m (PARI, block-indexed)      known (Table 3)        diff");
a5 = series_blocks(5, 30); r5 = rho_from_blocks(a5, 5); printf("%3d   %20.10f   %20.6f   %10.6f\n", 5, r5, known5, r5-known5);
a6 = series_blocks(6, 30); r6 = rho_from_blocks(a6, 6); printf("%3d   %20.10f   %20.6f   %10.6f\n", 6, r6, known6, r6-known6);
a7 = series_blocks(7, 30); r7 = rho_from_blocks(a7, 7); printf("%3d   %20.10f   %20.6f   %10.6f\n", 7, r7, known7, r7-known7);
a8 = series_blocks(8, 30); r8 = rho_from_blocks(a8, 8); printf("%3d   %20.10f   %20.6f   %10.6f\n", 8, r8, known8, r8-known8);
quit;
