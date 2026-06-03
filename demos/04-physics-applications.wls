#!/usr/bin/env wolframscript
(* ::Package:: *)

(* 04-physics-applications.wls -- the algebras at work in physics

   Three vignettes:
     1. SU(3) flavour and the eightfold way (quarks, mesons, baryons),
     2. SU(2) and angular momentum (spin j, dimension 2j+1, Casimir 2 j(j+1)),
     3. a numeric check that a built representation really obeys the algebra.

   Run (no front end needed):
     wolframscript -file demos/04-physics-applications.wls
     /Applications/Wolfram.app/Contents/MacOS/wolframscript -file demos/04-physics-applications.wls
*)

scriptDir = If[$InputFileName === "", Directory[], DirectoryName[$InputFileName]];
root = ParentDirectory[scriptDir];
PacletDirectoryLoad[root];
Needs["ClassicalLieAlgebra`"];

section[s_] := (Print[]; Print["== ", s, " =="]);
p[lbl_, x_] := Print["  ", lbl, " = ", x];
check[lbl_, cond_] := Print[If[TrueQ[cond], "  [ok]   ", "  [FAIL] "], lbl];
bracket[a_, b_] := a . b - b . a;
dim[g_, w_] := RepresentationDimension[Irrep[g, w]];
(* fixed-width text columns -- PaddedForm/Grid need a front end to typeset, and
   ToString on a fraction builds it vertically; InputForm keeps it on one line *)
toStr[s_String] := s;
toStr[x_] := ToString[x, InputForm];
pad[x_, n_] := StringPadRight[toStr[x], n];
row5[a_, b_, c_, d_, e_] := Print["   ", pad[a, 7], pad[b, 9], pad[c, 7], pad[d, 10], toStr[e]];

chevalleyHolds[g_, w_] := Module[
   {m = RepresentationMatrices[Irrep[g, w]], a = CartanMatrix[g], r = Rank[g], z},
   z = 0 m["Cartan"][[1]];
   (And @@ Flatten@Table[
       bracket[m["Raising"][[i]], m["Lowering"][[j]]] == If[i == j, m["Cartan"][[i]], z],
       {i, r}, {j, r}]) &&
   (And @@ Flatten@Table[
       bracket[m["Cartan"][[i]], m["Raising"][[j]]] == a[[j, i]] m["Raising"][[j]],
       {i, r}, {j, r}])];

(* =====================================================================
   1. SU(3) flavour: the eightfold way
   ===================================================================== *)
section["SU(3) flavour and the eightfold way"];

Print["Gell-Mann and Ne'eman organized the hadrons into SU(3) flavour multiplets."];
Print["Each multiplet is an irrep, labelled by its Dynkin labels:"];
Print[];
p["quark triplet      3   = (1,0), dim", dim[SU[3], {1, 0}]];
p["antiquark triplet  3bar= (0,1), dim", dim[SU[3], {0, 1}]];
p["meson/baryon octet 8   = (1,1), dim", dim[SU[3], {1, 1}]];
p["baryon decuplet    10  = (3,0), dim", dim[SU[3], {3, 0}]];
Print[];
Print["The octet's weights are the eight states' (isospin, hypercharge) quantum"];
Print["numbers; the doubly-degenerate central weight is the pi0/eta (or Sigma0/Lambda) pair:"];
p["octet WeightSystem", WeightSystem[Irrep[SU[3], {1, 1}]]];

Print[];
Print["The hadron multiplets follow from tensoring quarks. The dimensions match"];
Print["the classic decompositions:"];
Print["  mesons    3 (x) 3bar = 8 (+) 1     : ", dim[SU[3], {1, 0}] dim[SU[3], {0, 1}],
   " = ", dim[SU[3], {1, 1}], " + ", dim[SU[3], {0, 0}]];
Print["  baryons   3 (x) 3 (x) 3 = 10 (+) 8 (+) 8 (+) 1 : ",
   dim[SU[3], {1, 0}]^3, " = ", dim[SU[3], {3, 0}], " + ", dim[SU[3], {1, 1}],
   " + ", dim[SU[3], {1, 1}], " + ", dim[SU[3], {0, 0}]];

(* =====================================================================
   2. SU(2) and angular momentum
   ===================================================================== *)
section["SU(2) and angular momentum"];

Print["The su(2) irrep with Dynkin label {2j} is the spin-j multiplet: dimension"];
Print["2j+1, and quadratic Casimir 2 j(j+1) in this package's normalization."];
Print[];
row5["j", "Dynkin", "dim", "Casimir", "2 j(j+1)"];
row5["----", "------", "---", "-------", "--------"];
(With[{j = #/2, n = #},
   row5[j, "{" <> ToString[n] <> "}", dim[SU[2], {n}],
      CasimirEigenvalue[Irrep[SU[2], {n}]], 2 j (j + 1)]] &) /@ {1, 2, 3, 4};

(* =====================================================================
   3. A built representation really obeys the algebra
   ===================================================================== *)
section["Numeric check: the octet matrices form the su(3) algebra"];

Print["RepresentationMatrices returns explicit numbers. We can verify directly"];
Print["that they satisfy the defining Chevalley relations of su(3):"];
Print["  [E_i, F_j] = delta_ij H_i,   [H_i, E_j] = A_ji E_j."];
m = RepresentationMatrices[Irrep[SU[3], {1, 1}]];
p["[E1,F1] - H1 == 0", bracket[m["Raising"][[1]], m["Lowering"][[1]]] - m["Cartan"][[1]] == 0 m["Cartan"][[1]]];
p["[E1,F2]    == 0", bracket[m["Raising"][[1]], m["Lowering"][[2]]] == 0 m["Cartan"][[1]]];
p["all su(3) Chevalley relations hold", chevalleyHolds[SU[3], {1, 1}]];
Print[];
Print["The same engine and the same check work for spinors and for sp:"];
p["so(5) spinor (0,1) obeys its algebra", chevalleyHolds[SO[5], {0, 1}]];
p["sp(4) defining (1,0) obeys its algebra", chevalleyHolds[Sp[4], {1, 0}]];

(* =====================================================================
   self-checks
   ===================================================================== *)
section["self-checks"];
check["3 (x) 3bar = 8 (+) 1 by dimension",
   dim[SU[3], {1, 0}] dim[SU[3], {0, 1}] == dim[SU[3], {1, 1}] + dim[SU[3], {0, 0}]];
check["3 (x) 3 (x) 3 = 10 (+) 8 (+) 8 (+) 1 by dimension",
   dim[SU[3], {1, 0}]^3 == dim[SU[3], {3, 0}] + 2 dim[SU[3], {1, 1}] + dim[SU[3], {0, 0}]];
check["spin-j dimension is 2j+1", And @@ (dim[SU[2], {#}] == # + 1 & /@ {1, 2, 3, 4})];
check["spin-j Casimir is 2 j(j+1)",
   And @@ (CasimirEigenvalue[Irrep[SU[2], {#}]] == With[{j = #/2}, 2 j (j + 1)] & /@ {1, 2, 3, 4})];
check["octet matrices satisfy su(3) Chevalley relations", chevalleyHolds[SU[3], {1, 1}]];
check["sp(4) defining satisfies its Chevalley relations", chevalleyHolds[Sp[4], {1, 0}]];
Print[];
Print["done."];
