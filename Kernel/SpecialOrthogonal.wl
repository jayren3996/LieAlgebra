BeginPackage["ClassicalLieAlgebra`SpecialOrthogonal`"];
Needs["ClassicalLieAlgebra`"];
Needs["ClassicalLieAlgebra`Common`"];
Begin["`Private`"];

(* ---- Standard (antisymmetric imaginary) basis ---- *)
soStandardOne[m_, i_, j_] := Module[{mat = ConstantArray[0, {m, m}], a, b},
  {a, b} = If[i < j, {i, j}, {j, i}];
  mat[[a, b]] -= I;
  mat[[b, a]] += I;
  mat
];

soStandard[m_] := Module[{list = {}, i, j},
  For[j = 2, j <= m, j++,
    For[i = 1, i < j, i++,
      AppendTo[list, soStandardOne[m, i, j]];
    ];
  ];
  list
];

(* ---- Cartan-Weyl basis (antisymmetric realization) ---- *)
soCWH[m_, i_] := Module[{mat = ConstantArray[0, {m, m}]},
  mat[[2 i - 1, 2 i    ]] = -I;
  mat[[2 i,     2 i - 1]] = +I;
  mat
];

soCWE[{m_, 1}, i_, j_] := Module[{mat = ConstantArray[0, {m, m}], a, b},
  {a, b} = If[i < j, {i, j}, {j, i}];
  mat[[2 a,     2 b - 1]] = +1/2;
  mat[[2 a - 1, 2 b - 1]] = -I/2;
  mat[[2 a,     2 b    ]] = -I/2;
  mat[[2 a - 1, 2 b    ]] = -1/2;
  -I*mat + I*Transpose[mat]
];

soCWE[{m_, 2}, i_, j_] := Module[{mat = ConstantArray[0, {m, m}], a, b},
  {a, b} = If[i < j, {i, j}, {j, i}];
  mat[[2 a,     2 b - 1]] = +1/2;
  mat[[2 a - 1, 2 b - 1]] = -I/2;
  mat[[2 a,     2 b    ]] = +I/2;
  mat[[2 a - 1, 2 b    ]] = +1/2;
  -I*mat + I*Transpose[mat]
];

soCWE[{m_, 3}, i_] := Module[{mat = ConstantArray[0, {m, m}]},
  mat[[2 i,     m]] = +1/Sqrt[2];
  mat[[2 i - 1, m]] = -I/Sqrt[2];
  -I*mat + I*Transpose[mat]
];

soCartanWeyl[m_] := Module[{l, p, hh, ee = {}, i, j},
  {l, p} = QuotientRemainder[m, 2];
  hh = Table[soCWH[m, i], {i, l}];
  For[j = 2, j <= l, j++,
    For[i = 1, i < j, i++,
      AppendTo[ee, soCWE[{m, 1}, i, j]];
      AppendTo[ee, soCWE[{m, 2}, i, j]];
    ];
  ];
  If[p == 1,
    For[i = 1, i <= l, i++,
      AppendTo[ee, soCWE[{m, 3}, i]];
    ];
  ];
  {hh, ee, Conjugate @* Transpose /@ ee}
];

(* ---- Chevalley basis in antisymmetric realization (SOCh) ---- *)
soChH[m_?OddQ, i_] := Module[{mat = ConstantArray[0, {m, m}], l = (m - 1)/2},
  If[i == l,
    mat[[2 l - 1, 2 l]] = -2 I;
    mat[[2 l, 2 l - 1]] = +2 I;
    ,
    mat[[2 i - 1, 2 i    ]] = -I;
    mat[[2 i,     2 i - 1]] = +I;
    mat[[2 i + 1, 2 i + 2]] = +I;
    mat[[2 i + 2, 2 i + 1]] = -I;
  ];
  mat
];

soChH[m_?EvenQ, i_] := Module[{mat = ConstantArray[0, {m, m}], l = m/2},
  If[i == l,
    mat[[2 l - 1, 2 l    ]] = -I;
    mat[[2 l,     2 l - 1]] = +I;
    mat[[2 l - 3, 2 l - 2]] = -I;
    mat[[2 l - 2, 2 l - 3]] = +I;
    ,
    mat[[2 i - 1, 2 i    ]] = -I;
    mat[[2 i,     2 i - 1]] = +I;
    mat[[2 i + 1, 2 i + 2]] = +I;
    mat[[2 i + 2, 2 i + 1]] = -I;
  ];
  mat
];

soChE[m_?OddQ, i_] := Module[{mat = ConstantArray[0, {m, m}], l = (m - 1)/2},
  If[i == l,
    mat[[2 i,     m]] = +1;
    mat[[2 i - 1, m]] = -I;
    ,
    mat[[2 i,     2 i + 1]] = +1/2;
    mat[[2 i - 1, 2 i + 1]] = -I/2;
    mat[[2 i,     2 i + 2]] = -I/2;
    mat[[2 i - 1, 2 i + 2]] = -1/2;
  ];
  -I*mat + I*Transpose[mat]
];

soChE[m_?EvenQ, i_] := Module[{mat = ConstantArray[0, {m, m}], l = m/2},
  If[i == l,
    mat[[2 i - 2, 2 i - 1]] = +1/2;
    mat[[2 i - 3, 2 i - 1]] = -I/2;
    mat[[2 i - 2, 2 i    ]] = +I/2;
    mat[[2 i - 3, 2 i    ]] = +1/2;
    ,
    mat[[2 i,     2 i + 1]] =  1/2;
    mat[[2 i - 1, 2 i + 1]] = -I/2;
    mat[[2 i,     2 i + 2]] = -I/2;
    mat[[2 i - 1, 2 i + 2]] = -1/2;
  ];
  -I*mat + I*Transpose[mat]
];

soChevalleyAntisym[m_] := Module[{l, hh, ee},
  l = Quotient[m, 2];
  hh = Table[soChH[m, i], {i, l}];
  ee = Table[soChE[m, i], {i, l}];
  {hh, ee, Conjugate @* Transpose /@ ee}
];

(* ---- Change-of-basis matrix (antisymmetric -> diagonal) ----
   BasisTransform is, by definition, the change of basis that turns the antisymmetric
   Chevalley realization into the diagonal one. We compute it directly as the intertwiner
   B with  B . X_antisym . Inverse[B] == X_diagonal  for every Chevalley generator X
   (unique up to scale by Schur, since the defining rep is irreducible), then normalize B
   to be unitary. Solving  B . A_k == D_k . B  is the linear system
   (I (x) A_k^T - D_k (x) I) . vec(B) == 0  (row-major vec), stacked over all generators. *)
soBasis[m_] := Module[{anti, diag, sys, b0, c},
  anti = soChevalleyAntisym[m];
  diag = soChevalleyDiagonal[m];
  sys = Join @@ MapThread[
     KroneckerProduct[IdentityMatrix[m], Transpose[#1]] - KroneckerProduct[#2, IdentityMatrix[m]] &,
     {Join @@ anti, Join @@ diag}];
  b0 = ArrayReshape[First[NullSpace[sys]], {m, m}];
  c = (b0 . ConjugateTranspose[b0])[[1, 1]];
  b0 / Sqrt[c]
];

(* ---- Chevalley basis in diagonal realization (SOn) ---- *)
soH[m_?OddQ, i_] := Module[{mat = ConstantArray[0, {m, m}], l = (m - 1)/2},
  If[i == l,
    mat[[i,       i      ]] = +2;
    mat[[i + 2,   i + 2  ]] = -2;
    ,
    mat[[i,       i      ]] = +1;
    mat[[i + 1,   i + 1  ]] = -1;
    mat[[m - i,   m - i  ]] = +1;
    mat[[m - i + 1, m - i + 1]] = -1;
  ];
  mat
];

soH[m_?EvenQ, i_] := Module[{mat = ConstantArray[0, {m, m}], l = m/2},
  If[i == l,
    mat[[i - 1, i - 1]] = +1;
    mat[[i,     i    ]] = +1;
    mat[[i + 1, i + 1]] = -1;
    mat[[i + 2, i + 2]] = -1;
    ,
    mat[[i,       i      ]] = +1;
    mat[[i + 1,   i + 1  ]] = -1;
    mat[[m - i,   m - i  ]] = +1;
    mat[[m - i + 1, m - i + 1]] = -1;
  ];
  mat
];

soE[m_?OddQ, i_] := Module[{mat = ConstantArray[0, {m, m}], l = (m - 1)/2},
  If[i == l,
    mat[[i,     i + 1]] = +Sqrt[2];
    mat[[i + 1, i + 2]] = +Sqrt[2];
    ,
    mat[[i,     i + 1  ]] = +1;
    mat[[m - i, m - i + 1]] = +1;
  ];
  mat
];

soE[m_?EvenQ, i_] := Module[{mat = ConstantArray[0, {m, m}], l = m/2},
  If[i == l,
    mat[[i - 1, i + 1]] = +1;
    mat[[i,     i + 2]] = +1;
    ,
    mat[[i,     i + 1  ]] = +1;
    mat[[m - i, m - i + 1]] = +1;
  ];
  mat
];

soChevalleyDiagonal[m_] := Module[{l, hh, ee},
  l = Quotient[m, 2];
  hh = Table[soH[m, i], {i, l}];
  ee = Table[soE[m, i], {i, l}];
  {hh, ee, Transpose /@ ee}
];

(* ---- Association helper (local copy matches SpecialUnitary) ---- *)
toAssoc[{h_, e_, f_}] := <|"Cartan" -> h, "Raising" -> e, "Lowering" -> f|>;

(* ---- Matrix dimension from LieAlgebra head ---- *)
soMatDim[g_LieAlgebra] := Switch[g[[1]], "B", 2 g[[2]] + 1, "D", 2 g[[2]]];

(* ---- Option reader ---- *)
realizationOf[opts___] := OptionValue[Generators, {opts}, "Realization"];

(* ---- Generators dispatch for types B and D (so(n)) ---- *)
Generators[g_LieAlgebra /; MatchQ[g[[1]], "B" | "D"]] :=
  soStandard[soMatDim[g]];

Generators[g_LieAlgebra /; MatchQ[g[[1]], "B" | "D"], "Standard"] :=
  soStandard[soMatDim[g]];

Generators[g_LieAlgebra /; MatchQ[g[[1]], "B" | "D"], "Standard", OptionsPattern[]] :=
  soStandard[soMatDim[g]];

Generators[g_LieAlgebra /; MatchQ[g[[1]], "B" | "D"], "CartanWeyl", OptionsPattern[]] :=
  toAssoc[soCartanWeyl[soMatDim[g]]];

Generators[g_LieAlgebra /; MatchQ[g[[1]], "B" | "D"], "Chevalley", opts : OptionsPattern[]] :=
  Module[{rz = realizationOf[opts]},
    Switch[rz,
      "Diagonal",      toAssoc[soChevalleyDiagonal[soMatDim[g]]],
      "Antisymmetric", toAssoc[soChevalleyAntisym[soMatDim[g]]],
      _,               Message[Generators::badrealization, rz]; $Failed
    ]
  ];

BasisTransform[g_LieAlgebra /; MatchQ[g[[1]], "B" | "D"]] := soBasis[soMatDim[g]];

End[];
EndPackage[];
