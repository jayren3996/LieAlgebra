#!/usr/bin/env wolframscript
(* ::Package:: *)

(* 01-getting-started.wls -- a tour of ClassicalLieAlgebra

   Construct the classical algebras, read off their root-system data, and look
   at the defining generators in the three standard bases.

   Run it with a command-line kernel (no front end needed):

     wolframscript -file demos/01-getting-started.wls

   On macOS the system `wolframscript` may be broken; use the app binary:

     /Applications/Wolfram.app/Contents/MacOS/wolframscript -file demos/01-getting-started.wls
*)

scriptDir = If[$InputFileName === "", Directory[], DirectoryName[$InputFileName]];
root = ParentDirectory[scriptDir];
PacletDirectoryLoad[root];
Needs["ClassicalLieAlgebra`"];

section[s_] := (Print[]; Print["== ", s, " =="]);
p[lbl_, x_] := Print["  ", lbl, " = ", x];
(* print a matrix one row per line -- MatrixForm needs a front end to typeset *)
pm[lbl_, m_] := (Print["  ", lbl, " ="]; Scan[Print["    ", #] &, m]);
bracket[a_, b_] := a . b - b . a;

check[lbl_, cond_] := Print[If[TrueQ[cond], "  [ok]   ", "  [FAIL] "], lbl];

(* =====================================================================
   1. The four families and how to name an algebra
   ===================================================================== *)
section["Constructing algebras"];

Print["The classical simple Lie algebras come in four families. Each has a"];
Print["canonical form LieAlgebra[type, rank], and a familiar-name shorthand:"];
Print[];
p["SU[3]  (special unitary)   ", SU[3]];
p["SO[5]  (odd orthogonal, B) ", SO[5]];
p["SO[6]  (even orthogonal, D)", SO[6]];
p["Sp[4]  (symplectic, C)     ", Sp[4]];
Print[];
Print["The shorthands normalize to the canonical form, so SU[3] and"];
p["LieAlgebra[\"A\", 2]", LieAlgebra["A", 2]];
p["SU[3] === LieAlgebra[\"A\", 2]", SU[3] === LieAlgebra["A", 2]];

(* =====================================================================
   2. Root-system data
   ===================================================================== *)
section["Root-system data (SO[5] = B_2)"];

g = SO[5];
p["Rank", Rank[g]];
p["Dimension", LieAlgebraDimension[g]];
p["CartanMatrix", CartanMatrix[g]];
p["SimpleRoots (Euclidean)", SimpleRoots[g]];
p["PositiveRoots (Euclidean)", PositiveRoots[g]];
p["FundamentalWeights", FundamentalWeights[g]];
Print[];
Print["The dimension is (# positive roots) x 2 + rank:"];
p["2 * #posRoots + rank", 2 Length[PositiveRoots[g]] + Rank[g]];

(* =====================================================================
   3. Generators in three bases
   ===================================================================== *)
section["Generators of su(3) in three bases"];

Print["Generators[g] returns the defining-representation generators."];
Print["For su(3) these are the eight generators T_a = lambda_a/2 (3x3):"];
gens = Generators[SU[3]];
p["count", Length[gens]];
p["each is", Dimensions[First[gens]]];
pm["the first generator, lambda_1/2", First[gens]];

Print[];
Print["CartanWeyl[g] and Chevalley[g] return an association keyed by"];
Print["\"Cartan\", \"Raising\", \"Lowering\":"];
cw = CartanWeyl[SU[3]];
ch = Chevalley[SU[3]];
p["CartanWeyl keys", Keys[cw]];
p["counts (Cartan/Raising/Lowering)", Length /@ Values[cw]];
Print[];
Print["The Chevalley basis is the convenient one for building representations."];
pm["H1 (Chevalley)", ch["Cartan"][[1]]];
pm["E1 (Chevalley, a raising operator)", ch["Raising"][[1]]];

Print[];
Print["These really are generators of the same algebra: [E1, F1] = H1."];
p["[E1,F1] == H1", bracket[ch["Raising"][[1]], ch["Lowering"][[1]]] == ch["Cartan"][[1]]];

Print[];
Print["For so/sp, BasisTransform[g] conjugates the antisymmetric realization"];
Print["into the diagonal one. For SO[5] it is a 5x5 matrix:"];
p["Dimensions[BasisTransform[SO[5]]]", Dimensions[BasisTransform[SO[5]]]];

(* =====================================================================
   self-checks
   ===================================================================== *)
section["self-checks"];
check["SU[3] canonicalizes to A_2", SU[3] === LieAlgebra["A", 2]];
check["Sp[4] canonicalizes to C_2", Sp[4] === LieAlgebra["C", 2]];
check["Rank[SO[5]] == 2", Rank[SO[5]] == 2];
check["dim su(3) == 8", LieAlgebraDimension[SU[3]] == 8];
check["8 Gell-Mann generators", Length[Generators[SU[3]]] == 8];
check["B_2 Cartan matrix", CartanMatrix[SO[5]] == {{2, -2}, {-1, 2}}];
check["Chevalley [E1,F1] == H1", bracket[ch["Raising"][[1]], ch["Lowering"][[1]]] == ch["Cartan"][[1]]];
Print[];
Print["done."];
