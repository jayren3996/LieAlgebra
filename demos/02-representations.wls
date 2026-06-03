#!/usr/bin/env wolframscript
(* ::Package:: *)

(* 02-representations.wls -- the representation engine

   Irrep[g, w] builds the irreducible representation of g with highest weight w
   (Dynkin labels). From it you can read off the dimension, the weight system,
   the Casimir eigenvalue, and the explicit Chevalley generator matrices -- for
   every classical type, including the orthogonal/symplectic spinors that do not
   live in tensor powers of the defining representation.

   Run (no front end needed):
     wolframscript -file demos/02-representations.wls
     /Applications/Wolfram.app/Contents/MacOS/wolframscript -file demos/02-representations.wls
*)

scriptDir = If[$InputFileName === "", Directory[], DirectoryName[$InputFileName]];
root = ParentDirectory[scriptDir];
PacletDirectoryLoad[root];
Needs["ClassicalLieAlgebra`"];

section[s_] := (Print[]; Print["== ", s, " =="]);
p[lbl_, x_] := Print["  ", lbl, " = ", x];
check[lbl_, cond_] := Print[If[TrueQ[cond], "  [ok]   ", "  [FAIL] "], lbl];
bracket[a_, b_] := a . b - b . a;
pad[x_, n_] := StringPadRight[ToString[x], n];

(* Do the built generators satisfy the Chevalley relations of g?
   [E_i, F_j] = delta_ij H_i   and   [H_i, E_j] = A_ji E_j  (Cartan matrix A). *)
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
   1. One irrep, end to end: the su(3) adjoint (octet)
   ===================================================================== *)
section["The su(3) adjoint, (1,1) -- the 'octet'"];

ir = Irrep[SU[3], {1, 1}];
p["Irrep", ir];
p["HighestWeight (Dynkin labels)", HighestWeight[ir]];
p["RepresentationDimension (Weyl formula)", RepresentationDimension[ir]];
p["CasimirEigenvalue", CasimirEigenvalue[ir]];
Print[];
Print["WeightSystem gives <|weight -> multiplicity|>, weights as Dynkin labels."];
Print["The zero weight is doubly degenerate, which is why the octet has a"];
Print["two-dimensional center:"];
ws = WeightSystem[ir];
p["WeightSystem", ws];
p["sum of multiplicities", Total[Values[ws]]];
p["multiplicity of weight {0,0}", ws[{0, 0}]];

Print[];
Print["RepresentationMatrices gives the Chevalley generators as matrices in the"];
Print["irrep, in an orthonormal basis (E_i = ConjugateTranspose[F_i], H_i diagonal):"];
m = RepresentationMatrices[ir];
p["keys", Keys[m]];
p["one H, E, F per simple root", Length /@ Values[m]];
p["each matrix is", Dimensions[m["Cartan"][[1]]]];
p["H1 diagonal", Diagonal[m["Cartan"][[1]]]];
p["E1 == ConjugateTranspose[F1]", m["Raising"][[1]] == ConjugateTranspose[m["Lowering"][[1]]]];
p["Chevalley relations hold", chevalleyHolds[SU[3], {1, 1}]];

(* =====================================================================
   2. A dimension table across the families
   ===================================================================== *)
section["Dimensions and Casimirs across su / so / sp"];

rows = {
   {"su(2)  spin-1/2", SU[2], {1}},
   {"su(2)  spin-1  ", SU[2], {2}},
   {"su(3)  fund 3  ", SU[3], {1, 0}},
   {"su(3)  adj 8   ", SU[3], {1, 1}},
   {"su(3)  10      ", SU[3], {3, 0}},
   {"so(5)  vector 5", SO[5], {1, 0}},
   {"so(5)  spinor 4", SO[5], {0, 1}},
   {"sp(4)  defng 4 ", Sp[4], {1, 0}},
   {"so(7)  spinor 8", SO[7], {0, 0, 1}}
   };
Print["  ", pad["label", 18], pad["dim", 6], "Casimir"];
Print["  ", pad["-----", 18], pad["---", 6], "-------"];
(With[{ir2 = Irrep[#[[2]], #[[3]]]},
    Print["  ", pad[#[[1]], 18], pad[RepresentationDimension[ir2], 6],
       CasimirEigenvalue[ir2]]] &) /@ rows;

(* =====================================================================
   3. The spinor representations
   ===================================================================== *)
section["A spinor: so(5), highest weight (0,1)"];

Print["The so(5) spinor is 4-dimensional and does NOT appear in any tensor power"];
Print["of the 5-dimensional vector representation. The abstract highest-weight"];
Print["(Shapovalov) construction reaches it anyway:"];
sp5 = Irrep[SO[5], {0, 1}];
p["dimension", RepresentationDimension[sp5]];
ms = RepresentationMatrices[sp5];
p["generator matrix size", Dimensions[ms["Cartan"][[1]]]];
p["Chevalley relations hold", chevalleyHolds[SO[5], {0, 1}]];
p["so(7) spinor relations hold", chevalleyHolds[SO[7], {0, 0, 1}]];

(* =====================================================================
   self-checks
   ===================================================================== *)
section["self-checks"];
check["octet dimension == 8", RepresentationDimension[ir] == 8];
check["dimension == sum of weight multiplicities", RepresentationDimension[ir] == Total[Values[ws]]];
check["octet Casimir == 6", CasimirEigenvalue[ir] == 6];
check["zero weight has multiplicity 2", ws[{0, 0}] == 2];
check["orthonormal basis: E1 == F1^dagger", m["Raising"][[1]] == ConjugateTranspose[m["Lowering"][[1]]]];
check["octet satisfies su(3) Chevalley relations", chevalleyHolds[SU[3], {1, 1}]];
check["so(5) spinor dimension == 4", RepresentationDimension[sp5] == 4];
check["so(5) spinor satisfies Chevalley relations", chevalleyHolds[SO[5], {0, 1}]];
Print[];
Print["done."];
