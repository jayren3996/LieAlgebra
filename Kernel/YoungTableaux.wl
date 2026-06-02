BeginPackage["ClassicalLieAlgebra`YoungTableaux`"];
Needs["ClassicalLieAlgebra`"];
Begin["`Private`"];

(* ---- Permutation helpers (private) ---- *)
CycleDotAtom[c1_Cycles, c2_Cycles] := PermutationProduct[c2, c1];
CycleDotAtom[a_*c1_Cycles, c2_Cycles] := a * PermutationProduct[c2, c1];
CycleDotAtom[c1_Cycles, b_*c2_Cycles] := b * PermutationProduct[c2, c1];
CycleDotAtom[a_*c1_Cycles, b_*c2_Cycles] := a*b * PermutationProduct[c2, c1];

CycleDot[a_, b_] := Distribute[CycleDotAtom[a, b]];
CycleProduct[a__] := Fold[CycleDot, {a}];

(* ---- Young tableau helpers (private) ---- *)
TableauTranspose[t_List] := Module[
  {l = Length[t], li = Length /@ t, p = 1, out = {}, i, j},
  For[i = 1, i <= l, i++,
    j = l - i + 1;
    out = Join[out, Transpose @ t[[1 ;; j, p ;; li[[j]]]]];
    p = li[[j]] + 1;
  ];
  out
];
ShowTableau[t_List] := Grid[Map[Item[#, Frame -> True]&, #]& /@ t];
TableauTranspose[t_Tableau] := Tableau @ TableauTranspose @ t[[1]];

(* TableauForm for bare Tableau (grid display of rows) *)
TableauForm[t_Tableau] := ShowTableau @ t[[1]];

(* ---- Young symmetrizer (private) ---- *)
iParity[p_] := (-1)^(Total[Length /@ p[[1]]] + Length[p[[1]]]);
iPermutationGroupElements[l_] := GroupElements @ PermutationGroup[Cycles[{#}]& /@ {l[[1 ;; 2]], l}];
iPermuteElements[l_, 1] := Total @ iPermutationGroupElements[l];
iPermuteElements[l_, 2] := (Dot[iParity /@ #, #]&) @ iPermutationGroupElements[l];
YoungSymmetrizer[t_List] := Module[
  {p = Cycles[{{}}], q = Cycles[{{}}], tt = TableauTranspose[t], l, i},
  For[i = 1, i <= Length[t], i++,
    If[Length[t[[i]]] < 2, Break[]];
    p = CycleDot[iPermuteElements[t[[i]], 1], p];
  ];
  l = Length[tt];
  For[i = 1, i <= l, i++,
    If[Length[tt[[i]]] < 2, Break[]];
    q = CycleDot[iPermuteElements[tt[[i]], 2], q];
  ];
  CycleDot[p, q]
];
Symmetrizer[t_Tableau] := YoungSymmetrizer @ t[[1]];

(* ---- Tensor permutation (private atoms; public TableauPermute) ---- *)
TensorPermuteAtom[c_Cycles, v_Psi] := Permute[v, c];
TensorPermuteAtom[a_*c_Cycles, v_Psi] := a * Permute[v, c];
TensorPermuteAtom[c_Cycles, b_*v_Psi] := b * Permute[v, c];
TensorPermuteAtom[a_*c_Cycles, b_*v_Psi] := a*b * Permute[v, c];
iTensorPermute[a_, b_Psi] := Distribute @ TensorPermuteAtom[a, b];
TableauPermute[t_Tableau, v_Psi] := iTensorPermute[Symmetrizer[t], v];

(* ---- Tensor norm (private iSquarePsi / iTensorNorm; public TensorNorm) ---- *)
iSquarePsi[p_Psi] := 1;
iSquarePsi[a_*p_Psi] := Abs[a]^2;
iTensorNorm[p_Psi] := 1;
iTensorNorm[a_*p_Psi] := Abs[a];
iTensorNorm[p_] := Sqrt @ Total[iSquarePsi /@ List @@ p];

(* Public TensorNorm Expands first so scalars distribute over sums *)
TensorNorm[expr_] := iTensorNorm[Expand[expr]];

(* ---- Tensor inner product (public TensorDot) ---- *)
iDotPsi[p1_Psi, p2_Psi] := If[SameQ[p1, p2], 1, 0];
iDotPsi[a_*p1_Psi, p2_Psi] := If[SameQ[p1, p2], Conjugate[a], 0];
iDotPsi[p1_Psi, b_*p2_Psi] := If[SameQ[p1, p2], b, 0];
iDotPsi[a_*p1_Psi, b_*p2_Psi] := If[SameQ[p1, p2], Conjugate[a]*b, 0];
TensorDot[p1_, p2_] := Distribute @ iDotPsi[p1, p2];

(* ---- TensorTableau -> tensor conversion (private helpers) ---- *)
ListToTensor[t_List] := Module[
  {l = Length[t], li = Length /@ t, i, p = 1, tab = {}},
  For[i = 1, i <= l, i++,
    AppendTo[tab, Range[p, p + li[[i]] - 1]];
    p += li[[i]];
  ];
  TableauPermute[Tableau[tab], Psi @@ Flatten[t]]
];

(* Private inner ToTensor rules that handle atoms *)
iToTensor[t_TensorTableau] := ListToTensor @ t[[1]];
iToTensor[a_*t_TensorTableau] := a * ListToTensor @ t[[1]];
iToTensor[t1_ + t2_] := iToTensor[t1] + iToTensor[t2];

(* Public ToTensor Expands first so scalar*sum distributes *)
ToTensor[expr_] := iToTensor[Expand[expr]];

(* ---- TableauForm for TensorTableau (private inner rules; public entry point) ---- *)
iTableauForm[t_TensorTableau] := ShowTableau @ t[[1]];
iTableauForm[a_*t_TensorTableau] := a * ShowTableau @ t[[1]];
iTableauForm[t1_ + t2_] := iTableauForm[t1] + iTableauForm[t2];

TableauForm[expr_TensorTableau] := iTableauForm[expr];
TableauForm[expr_] := iTableauForm[Expand[expr]];

(* ---- Public TableauDot ---- *)
TableauDot[t1_, t2_] := TensorDot[ToTensor[t1], ToTensor[t2]];

(* ---- Normalization & Orthogonalization ---- *)
TableauNormalization[t_] := Module[{e = Expand[t]},
  e / iTensorNorm[Expand @ ToTensor[e]]
];

TableauOrthogonalization[t1_, t2_] := Module[
  {v1 = Expand @ ToTensor[t1], v2 = Expand @ ToTensor[t2], ov},
  ov = TensorDot[v1, v2] / iTensorNorm[v1]^2;
  {t1, Expand[t2 - ov*t1]}
];

End[];
EndPackage[];
