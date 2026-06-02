BeginPackage["ClassicalLieAlgebra`Representations`"];
Needs["ClassicalLieAlgebra`"];
Needs["ClassicalLieAlgebra`Common`"];
Needs["ClassicalLieAlgebra`Weights`"];
Begin["`Private`"];

validWeightQ[g_, w_] := ListQ[w] && Length[w] === Rank[g] && AllTrue[w, IntegerQ[#] && # >= 0 &];

(* fire only on INVALID input; valid Irrep[g,w] stays inert. Guard on g_LieAlgebra so the
   literal-LieAlgebra-pattern validation problem (Phase 1) does not bite. *)
Irrep[g_LieAlgebra, w_] /; ! validWeightQ[g, w] := (Message[Irrep::badweight, w, g, Rank[g]]; $Failed);

HighestWeight[Irrep[g_LieAlgebra, w_]] := w;

RepresentationDimension[Irrep[g_LieAlgebra, w_]] := ClassicalLieAlgebra`Weights`Private`weylDim[g, w];

CasimirEigenvalue[Irrep[g_LieAlgebra, w_]] := ClassicalLieAlgebra`Weights`Private`casimir[g, w];

WeightSystem[Irrep[g_LieAlgebra, w_]] := ClassicalLieAlgebra`Weights`Private`weightSystemDynkin[g, w];

(* ── Task 7: Shapovalov contravariant form on lowering-words ───────────────── *)

(* Weight of a word in Dynkin-label coordinates.
   w_λ has weight λ; applying f_i shifts the weight by -α_i, which in
   Dynkin coordinates subtracts the i-th row of the Cartan matrix. *)
wordWeight[g_, lam_, word_List] :=
  lam - Total[CartanMatrix[g][[#]] & /@ word, 1];

(* applyF[i, vec]: prepend root index i to every word key (= act with f_i). *)
applyF[i_Integer, vec_Association] :=
  KeyMap[Prepend[#, i] &, vec];

(* applyH[g, lam, i, word]: h_i acts as the i-th Dynkin label of the word's weight. *)
applyH[g_, lam_, i_Integer, word_List] :=
  Module[{val = wordWeight[g, lam, word][[i]]},
    If[val === 0, <||>, <|word -> val|>]
  ];

(* applyE[g, lam, i, word]: single-word action of e_i.
   Uses e_i f_j = f_j e_i + δ_{ij} h_i and e_i v_λ = 0. *)
applyE[g_, lam_, i_Integer, word_List] :=
  If[word === {},
    <||>,  (* e_i annihilates v_λ *)
    Module[{j = First[word], rest = Rest[word]},
      mergeVecs[
        applyF[j, applyE[g, lam, i, rest]],
        If[i === j, applyH[g, lam, i, rest], <||>]
      ]
    ]
  ];

(* Extend applyE linearly over a vector (Association word -> coeff). *)
applyE[g_, lam_, i_Integer, vec_Association] :=
  mergeVecs @@ KeyValueMap[Function[{word, coef},
    If[coef === 0, <||>, Map[Times[coef, #] &, applyE[g, lam, i, word]]]
  ], vec];

(* Merge a sequence of word-vectors by adding coefficients; drop zeros. *)
mergeVecs[vecs__Association] :=
  DeleteCases[Merge[{vecs}, Total], 0];
mergeVecs[] := <||>;

(* shapWord[word, vec]: ⟨word, vec⟩ where word is a single lowering-word.
   ⟨{}, v⟩ = coefficient of {} in v.
   ⟨{i, rest}, v⟩ = ⟨rest, e_i v⟩  (contravariant: ⟨f_i u, w⟩ = ⟨u, e_i w⟩). *)
shapWord[g_, lam_, {}, vec_Association] :=
  Lookup[vec, Key[{}], 0];
shapWord[g_, lam_, word_List, vec_Association] :=
  shapWord[g, lam, Rest[word], applyE[g, lam, First[word], vec]];

(* shapovalov[g, lam, vec1, vec2]: bilinear contravariant form ⟨vec1, vec2⟩. *)
shapovalov[g_, lam_, vec1_Association, vec2_Association] :=
  Total[KeyValueMap[Function[{word, coef},
    coef * shapWord[g, lam, word, vec2]
  ], vec1]];

End[];
EndPackage[];
