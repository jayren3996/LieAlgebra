BeginPackage["ClassicalLieAlgebra`SpecialUnitary`"];
Needs["ClassicalLieAlgebra`"];
Needs["ClassicalLieAlgebra`Common`"];
Begin["`Private`"];

(* ---- Standard (Hermitian traceless) basis ---- *)
(* Type-1: symmetric off-diagonal pair *)
suT[m_, 1, i_, j_] := Module[{mat = ConstantArray[0, {m, m}], a, b},
  {a, b} = If[i < j, {i, j}, {j, i}];
  mat[[a, b]] += 1/2;
  mat[[b, a]] += 1/2;
  mat
];
(* Type-2: antisymmetric off-diagonal pair *)
suT[m_, 2, i_, j_] := Module[{mat = ConstantArray[0, {m, m}], a, b},
  {a, b} = If[i < j, {i, j}, {j, i}];
  mat[[a, b]] -= I/2;
  mat[[b, a]] += I/2;
  mat
];
(* Type-3: diagonal generator *)
suT[m_, 3, k_] := Module[{d = ConstantArray[0, m]},
  d[[1 ;; k - 1]] = 1/Sqrt[2 k (k - 1)];
  d[[k]] = -Sqrt[(k - 1)/(2 k)];
  DiagonalMatrix[d]
];

suStandard[m_] := Module[{list = {}, i, j},
  For[j = 2, j <= m, j++,
    For[i = 1, i < j, i++,
      AppendTo[list, suT[m, 1, i, j]];
      AppendTo[list, suT[m, 2, i, j]];
    ];
    AppendTo[list, suT[m, 3, j]];
  ];
  list
];

(* ---- Cartan-Weyl basis ---- *)
suCWH[m_, i_] := Module[{d = ConstantArray[0, m], jj},
  jj = m + 1 - i;
  d[[1 ;; jj - 1]] = 1/Sqrt[jj (jj - 1)];
  d[[jj]] = -Sqrt[(jj - 1)/jj];
  DiagonalMatrix[d]
];

suCWE[m_, i_, j_] := Module[{mat = ConstantArray[0, {m, m}], a, b},
  {a, b} = If[i < j, {i, j}, {j, i}];
  mat[[a, b]] = 1;
  mat
];

suCWF[m_, i_, j_] := Transpose[suCWE[m, i, j]];

suCartanWeyl[m_] := Module[{hl, el = {}, fl, i, j},
  hl = Table[suCWH[m, i], {i, m - 1}];
  For[j = 2, j <= m, j++,
    For[i = 1, i < j, i++,
      AppendTo[el, suCWE[m, i, j]];
    ];
  ];
  fl = Transpose /@ el;
  {hl, el, fl}
];

(* ---- Chevalley basis ---- *)
suChevH[m_, i_] := Module[{d = ConstantArray[0, m]},
  d[[i]] = 1;
  d[[i + 1]] = -1;
  DiagonalMatrix[d]
];

suChevE[m_, i_] := suCWE[m, i, i + 1];
suChevF[m_, i_] := suCWF[m, i, i + 1];

suChevalley[m_] := Module[{hl, el, fl},
  hl = Table[suChevH[m, i], {i, m - 1}];
  el = Table[suChevE[m, i], {i, m - 1}];
  fl = Transpose /@ el;
  {hl, el, fl}
];

(* ---- Association helper ---- *)
toAssoc[{h_, e_, f_}] := <|"Cartan" -> h, "Raising" -> e, "Lowering" -> f|>;

(* ---- Generators dispatch for type A (su(n)) ---- *)
Generators[HoldPattern[LieAlgebra["A", r_]]] := suStandard[r + 1];
Generators[HoldPattern[LieAlgebra["A", r_]], "Standard"] := suStandard[r + 1];
Generators[HoldPattern[LieAlgebra["A", r_]], "Standard", OptionsPattern[]] := suStandard[r + 1];
Generators[HoldPattern[LieAlgebra["A", r_]], "CartanWeyl", OptionsPattern[]] := toAssoc[suCartanWeyl[r + 1]];
Generators[HoldPattern[LieAlgebra["A", r_]], "Chevalley", OptionsPattern[]] := toAssoc[suChevalley[r + 1]];

End[];
EndPackage[];
