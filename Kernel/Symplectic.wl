BeginPackage["ClassicalLieAlgebra`Symplectic`"];
Needs["ClassicalLieAlgebra`"];
Needs["ClassicalLieAlgebra`Common`"];
Begin["`Private`"];

(* ---- Standard basis (spStandard): anti-Hermitian compact form ---- *)
(* Cross-pair generators for 1 <= a < b <= l (matrix indices 2a-1,2a,2b-1,2b) *)
spStandardCross[m_, a_, b_] := Module[{mat1, mat2, mat3, mat4},
  (* type 1: real antisymmetric *)
  mat1 = ConstantArray[0, {m, m}];
  mat1[[2 a - 1, 2 b - 1]] = +1;
  mat1[[2 a,     2 b    ]] = +1;
  mat1[[2 b - 1, 2 a - 1]] = -1;
  mat1[[2 b,     2 a    ]] = -1;
  (* type 2: real mixed-antisymmetric *)
  mat2 = ConstantArray[0, {m, m}];
  mat2[[2 a - 1, 2 b    ]] = +1;
  mat2[[2 a,     2 b - 1]] = -1;
  mat2[[2 b - 1, 2 a    ]] = +1;
  mat2[[2 b,     2 a - 1]] = -1;
  (* type 3: imaginary symmetric *)
  mat3 = ConstantArray[0, {m, m}];
  mat3[[2 a - 1, 2 b    ]] = +I;
  mat3[[2 a,     2 b - 1]] = +I;
  mat3[[2 b - 1, 2 a    ]] = +I;
  mat3[[2 b,     2 a - 1]] = +I;
  (* type 4: imaginary antisymmetric *)
  mat4 = ConstantArray[0, {m, m}];
  mat4[[2 a - 1, 2 b - 1]] = -I;
  mat4[[2 a,     2 b    ]] = +I;
  mat4[[2 b - 1, 2 a - 1]] = -I;
  mat4[[2 b,     2 a    ]] = +I;
  {mat1, mat2, mat3, mat4}
];

(* Within-pair generators for site i *)
spStandardSite[m_, i_] := Module[{mat5, mat6, mat7},
  (* type 5: imaginary diagonal *)
  mat5 = ConstantArray[0, {m, m}];
  mat5[[2 i - 1, 2 i - 1]] = +I;
  mat5[[2 i,     2 i    ]] = -I;
  (* type 6: real antisymmetric *)
  mat6 = ConstantArray[0, {m, m}];
  mat6[[2 i - 1, 2 i]] = +1;
  mat6[[2 i,     2 i - 1]] = -1;
  (* type 7: imaginary symmetric *)
  mat7 = ConstantArray[0, {m, m}];
  mat7[[2 i - 1, 2 i]] = +I;
  mat7[[2 i,     2 i - 1]] = +I;
  {mat5, mat6, mat7}
];

spStandard[m_] := Module[{l = m/2, list = {}, a, b, i},
  For[a = 1, a <= l, a++,
    For[b = a + 1, b <= l, b++,
      list = Join[list, spStandardCross[m, a, b]];
    ];
  ];
  For[i = 1, i <= l, i++,
    list = Join[list, spStandardSite[m, i]];
  ];
  list
];

(* ---- Cartan-Weyl basis (spCartanWeyl) ---- *)
spCWH[m_, i_] := Module[{mat = ConstantArray[0, {m, m}]},
  mat[[2 i - 1, 2 i - 1]] = 1/Sqrt[2];
  mat[[2 i,     2 i    ]] = -1/Sqrt[2];
  mat
];

spCWE[m_, type_, i_, j_] := Module[{mat = ConstantArray[0, {m, m}], a, b},
  {a, b} = If[i < j, {i, j}, {j, i}];
  Switch[type,
    1,
      mat[[2 a - 1, 2 b - 1]] = +1/Sqrt[2];
      mat[[2 b,     2 a    ]] = -1/Sqrt[2];,
    2,
      mat[[2 a - 1, 2 b]] = 1/Sqrt[2];
      mat[[2 b - 1, 2 a]] = 1/Sqrt[2];
  ];
  mat
];

spCWElong[m_, i_] := Module[{mat = ConstantArray[0, {m, m}]},
  mat[[2 i - 1, 2 i]] = 1;
  mat
];

spCartanWeyl[m_] := Module[{l = m/2, hh, ee = {}, i, j},
  hh = Table[spCWH[m, i], {i, l}];
  For[j = 2, j <= l, j++,
    For[i = 1, i < j, i++,
      AppendTo[ee, spCWE[m, 1, i, j]];
      AppendTo[ee, spCWE[m, 2, i, j]];
    ];
  ];
  For[i = 1, i <= l, i++,
    AppendTo[ee, spCWElong[m, i]];
  ];
  {hh, ee, Conjugate @* Transpose /@ ee}
];

(* ---- Change-of-basis matrix (antisymmetric -> diagonal) — a permutation ---- *)
spBasis[m_] := Module[{l = m/2, b = ConstantArray[0, {m, m}], i},
  For[i = 1, i <= l, i++,
    b[[i,     2 i - 1]] = 1;
    b[[l + i, m + 2 - 2 i]] = 1;
  ];
  b
];

(* ---- Chevalley basis in antisymmetric realization (SpCh -> spChevalleyAntisym) ---- *)
spChH[m_, i_] := Module[{mat = ConstantArray[0, {m, m}]},
  If[i == m/2,
    mat[[2 i - 1, 2 i - 1]] = +1;
    mat[[2 i,     2 i    ]] = -1;
    ,
    mat[[2 i - 1, 2 i - 1]] = +1;
    mat[[2 i,     2 i    ]] = -1;
    mat[[2 i + 1, 2 i + 1]] = -1;
    mat[[2 i + 2, 2 i + 2]] = +1;
  ];
  mat
];

spChE[m_, i_] := Module[{mat = ConstantArray[0, {m, m}]},
  If[i == m/2,
    mat[[2 i - 1, 2 i]] = 1;
    ,
    mat[[2 i - 1, 2 i + 1]] = +1;
    mat[[2 i + 2, 2 i    ]] = -1;
  ];
  mat
];

spChevalleyAntisym[m_] := Module[{l = m/2, hh, ee},
  hh = Table[spChH[m, i], {i, l}];
  ee = Table[spChE[m, i], {i, l}];
  {hh, ee, Conjugate @* Transpose /@ ee}
];

(* ---- Chevalley basis in diagonal realization (Spn -> spChevalleyDiagonal) ---- *)
spH[m_, i_] := Module[{mat = ConstantArray[0, {m, m}]},
  If[i == m/2,
    mat[[i,     i    ]] = 1;
    mat[[i + 1, i + 1]] = -1;
    ,
    mat[[i,         i        ]] = 1;
    mat[[i + 1,     i + 1    ]] = -1;
    mat[[m - i,     m - i    ]] = 1;
    mat[[m - i + 1, m - i + 1]] = -1;
  ];
  mat
];

spE[m_, i_] := Module[{mat = ConstantArray[0, {m, m}]},
  If[i == m/2,
    mat[[i, i + 1]] = 1;
    ,
    mat[[i,     i + 1    ]] = +1;
    mat[[m - i, m - i + 1]] = -1;
  ];
  mat
];

spChevalleyDiagonal[m_] := Module[{l = m/2, hh, ee},
  hh = Table[spH[m, i], {i, l}];
  ee = Table[spE[m, i], {i, l}];
  {hh, ee, Conjugate @* Transpose /@ ee}
];

(* ---- Association helper ---- *)
toAssoc[{h_, e_, f_}] := <|"Cartan" -> h, "Raising" -> e, "Lowering" -> f|>;

(* ---- Matrix dimension from LieAlgebra head ---- *)
spMatDim[g_LieAlgebra] := 2 g[[2]];

(* ---- Option reader ---- *)
realizationOf[opts___] := OptionValue[Generators, {opts}, "Realization"];

(* ---- Generators dispatch for type C (sp(2n)) ---- *)
Generators[g_LieAlgebra /; MatchQ[g[[1]], "C"]] :=
  spStandard[spMatDim[g]];

Generators[g_LieAlgebra /; MatchQ[g[[1]], "C"], "Standard"] :=
  spStandard[spMatDim[g]];

Generators[g_LieAlgebra /; MatchQ[g[[1]], "C"], "Standard", OptionsPattern[]] :=
  spStandard[spMatDim[g]];

Generators[g_LieAlgebra /; MatchQ[g[[1]], "C"], "CartanWeyl", OptionsPattern[]] :=
  toAssoc[spCartanWeyl[spMatDim[g]]];

Generators[g_LieAlgebra /; MatchQ[g[[1]], "C"], "Chevalley", opts : OptionsPattern[]] :=
  Module[{rz = realizationOf[opts]},
    Switch[rz,
      "Diagonal",      toAssoc[spChevalleyDiagonal[spMatDim[g]]],
      "Antisymmetric", toAssoc[spChevalleyAntisym[spMatDim[g]]],
      _,               Message[Generators::badrealization, rz]; $Failed
    ]
  ];

BasisTransform[g_LieAlgebra /; MatchQ[g[[1]], "C"]] := spBasis[spMatDim[g]];

End[];
EndPackage[];
