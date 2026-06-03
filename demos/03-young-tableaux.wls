#!/usr/bin/env wolframscript
(* ::Package:: *)

(* 03-young-tableaux.wls -- many-body wavefunctions with Young tableaux

   This is the SU(3) computation from docs/walkthrough.md, made executable. An
   SU(3) irrep corresponds to a Young diagram of at most two rows; its states are
   the standard tableaux, read as symmetrized many-body wavefunctions. The toolkit
   builds them, takes inner products, normalizes, and orthogonalizes degenerate
   weight spaces -- exactly the Gram-Schmidt step you do in quantum mechanics.

   Run (no front end needed):
     wolframscript -file demos/03-young-tableaux.wls
     /Applications/Wolfram.app/Contents/MacOS/wolframscript -file demos/03-young-tableaux.wls
*)

scriptDir = If[$InputFileName === "", Directory[], DirectoryName[$InputFileName]];
root = ParentDirectory[scriptDir];
PacletDirectoryLoad[root];
Needs["ClassicalLieAlgebra`"];

section[s_] := (Print[]; Print["== ", s, " =="]);
p[lbl_, x_] := Print["  ", lbl, " = ", x];
check[lbl_, cond_] := Print[If[TrueQ[cond], "  [ok]   ", "  [FAIL] "], lbl];
T = TensorTableau;

(* =====================================================================
   1. The highest-weight tableau and the Young symmetrizer
   ===================================================================== *)
section["The Young symmetrizer"];

Print["The su(3) (1,1) irrep has a highest-weight tableau of shape [2,1], filled"];
Print["with 1s along the first row and 2s along the second:  {{1,1},{2}}."];
Print[];
Print["Tableau[rows] is the shape (its Young symmetrizer); Psi[i1,i2,...] is a"];
Print["basis tensor |i1 i2 ...>. TableauPermute applies the symmetrizer to a tensor."];
Print["For the column-standard shape {{1,2},{3}} acting on |1 1 2>:"];
ct = Tableau[{{1, 2}, {3}}];
sym = TableauPermute[ct, Psi[1, 1, 2]];
p["TableauPermute[Tableau[{{1,2},{3}}], Psi[1,1,2]]", sym];
Print["  (antisymmetric in the column {1,3}, symmetric in the row {1,2}.)"];

(* =====================================================================
   2. Tensor tableaux as compact many-body states
   ===================================================================== *)
section["Tensor tableaux and their wavefunctions"];

Print["TensorTableau[rows] is the symmetrized state in compact form; ToTensor"];
Print["expands it into the Psi basis, and TableauForm displays it as a grid"];
Print["(the grid renders best in a notebook front end)."];
psi = ToTensor[T[{{1, 2}, {2}}]];
p["ToTensor[TensorTableau[{{1,2},{2}}]]", psi];
Print["  TableauForm[TensorTableau[{{1,2},{2}}]] renders the same state as a"];
Print["  2-row Young grid in a notebook front end."];

Print[];
Print["Columns are antisymmetric, so a repeated entry in a column annihilates the"];
Print["state. The tableau {{1,2},{1}} has 1 twice in its first column:"];
p["ToTensor[TensorTableau[{{1,2},{1}}]]", ToTensor[T[{{1, 2}, {1}}]]];

(* =====================================================================
   3. Inner products and normalization
   ===================================================================== *)
section["Inner products and normalization"];

Print["TableauDot is the Hermitian inner product; TableauNormalization scales to"];
Print["unit norm. Take |a> = {{1,2},{2}} and normalize it:"];
a = T[{{1, 2}, {2}}];
na = TableauNormalization[a];
p["raw <a|a>", TableauDot[a, a]];
p["normalized |a>", na];
p["<a|a> after normalization", TableauDot[na, na]];
Print[];
Print["TensorDot / TensorNorm do the same at the level of raw Psi tensors:"];
p["TensorNorm[ToTensor[a]]", TensorNorm[ToTensor[a]]];

(* =====================================================================
   4. Orthogonalizing a degenerate weight space
   ===================================================================== *)
section["Orthogonalizing the degenerate (0,0) weight space"];

Print["In the (1,1) irrep the weight (0,0) is shared by two independent states."];
Print["They are not orthogonal, so -- as with Gram-Schmidt in QM -- we keep the"];
Print["first and orthogonalize the second against it. TableauOrthogonalization"];
Print["does exactly  {t1,  t2 - <t2,t1>/<t1,t1> t1}:"];
t1 = 2 T[{{1, 2}, {3}}] - T[{{1, 3}, {2}}];
t2 = T[{{1, 3}, {2}}] + T[{{1, 2}, {3}}];
{c, d} = TableauOrthogonalization[t1, t2];
p["kept state  c", c];
p["orthogonalized d", d];
p["<c|d>", TableauDot[c, d]];

(* =====================================================================
   5. A transition matrix element
   ===================================================================== *)
section["A transition matrix element"];

Print["Lowering normalized |a> with F2 reaches the superposition"];
Print["  |b> = ({{1,2},{3}} + {{1,3},{2}}) / Sqrt[6]   (see the walkthrough)."];
Print["Its overlap with a (0,0) basis state is the transition matrix element:"];
b = T[{{1, 2}, {3}}]/Sqrt[6] + T[{{1, 3}, {2}}]/Sqrt[6];
nd = TableauNormalization[2 T[{{1, 2}, {3}}] - T[{{1, 3}, {2}}]];
p["<b|b>", TableauDot[b, b]];
p["<d|b>  (transition element)", TableauDot[nd, b]];

(* =====================================================================
   self-checks
   ===================================================================== *)
section["self-checks"];
check["symmetrizer: TableauPermute matches the hand result",
   sym === 2 Psi[1, 1, 2] - Psi[1, 2, 1] - Psi[2, 1, 1]];
check["repeated column entry annihilates the state", ToTensor[T[{{1, 2}, {1}}]] === 0];
check["normalization gives unit norm", TableauDot[na, na] === 1];
check["TensorNorm agrees with TableauDot", TensorNorm[ToTensor[a]]^2 === TableauDot[a, a]];
check["orthogonalized states are orthogonal", TableauDot[c, d] === 0];
Print[];
Print["done."];
