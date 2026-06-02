Needs["ClassicalLieAlgebra`"];
VerificationTest[ HighestWeight[Irrep[SU[3], {1, 1}]], {1, 1}, TestID->"hw" ];
VerificationTest[ Head[Irrep[SU[3], {1, 1}]], Irrep, TestID->"irrep-inert" ];
VerificationTest[ Irrep[SU[3], {1, 1, 0}], $Failed, {Irrep::badweight}, TestID->"wrong-length" ];
VerificationTest[ Irrep[SU[3], {1, -1}], $Failed, {Irrep::badweight}, TestID->"negative-label" ];
VerificationTest[ Irrep[SO[5], {0, 1}], Irrep[LieAlgebra["B", 2], {0, 1}], TestID->"sugar-canonicalizes" ];
VerificationTest[ RepresentationDimension[Irrep[SU[3], {1, 0}]], 3, TestID->"dim-su3-fund" ];
VerificationTest[ RepresentationDimension[Irrep[SU[3], {1, 1}]], 8, TestID->"dim-su3-adjoint" ];
VerificationTest[ RepresentationDimension[Irrep[SO[5], {1, 0}]], 5, TestID->"dim-so5-vector" ];
VerificationTest[ RepresentationDimension[Irrep[SO[5], {0, 1}]], 4, TestID->"dim-so5-spinor" ];
VerificationTest[ RepresentationDimension[Irrep[Sp[4], {1, 0}]], 4, TestID->"dim-sp4-defining" ];
VerificationTest[ RepresentationDimension[Irrep[Sp[4], {0, 1}]], 5, TestID->"dim-sp4-omega2" ];
VerificationTest[ RepresentationDimension[Irrep[SU[4], {0, 0, 0}]], 1, TestID->"dim-trivial" ];
VerificationTest[ RepresentationDimension[Irrep[SO[7], {0, 0, 1}]], 8, TestID->"dim-so7-spinor" ];
VerificationTest[ CasimirEigenvalue[Irrep[SU[2], {1}]], 3/2, TestID->"cas-su2-fund" ];
VerificationTest[ CasimirEigenvalue[Irrep[SU[3], {1, 1}]], 6, TestID->"cas-su3-adjoint-2hv" ];
VerificationTest[ CasimirEigenvalue[Irrep[SU[3], {0, 0}]], 0, TestID->"cas-trivial" ];
VerificationTest[ CasimirEigenvalue[Irrep[Sp[4], {2, 0}]], 6, TestID->"cas-sp4-adjoint-2hv" ];
VerificationTest[ Total[Values[WeightSystem[Irrep[SU[3], {1, 0}]]]], 3, TestID->"ws-su3-fund-total" ];
VerificationTest[ Max[Values[WeightSystem[Irrep[SU[3], {1, 0}]]]], 1, TestID->"ws-su3-fund-mult1" ];
VerificationTest[ Total[Values[WeightSystem[Irrep[SU[3], {1, 1}]]]], 8, TestID->"ws-su3-adj-total" ];
VerificationTest[ WeightSystem[Irrep[SU[3], {1, 1}]][{0, 0}], 2, TestID->"ws-su3-adj-zeroweight-mult2" ];
VerificationTest[ Total[Values[WeightSystem[Irrep[SO[5], {1, 0}]]]], RepresentationDimension[Irrep[SO[5], {1, 0}]], TestID->"ws-so5-sum-eq-dim" ];
VerificationTest[ Total[Values[WeightSystem[Irrep[Sp[4], {0, 1}]]]], 5, TestID->"ws-sp4-sum" ];
VerificationTest[ Total[Values[WeightSystem[Irrep[SU[3], {2, 0}]]]], 6, TestID->"ws-su3-6plet" ];

(* Task 7: Shapovalov contravariant form *)
VerificationTest[ ClassicalLieAlgebra`Representations`Private`shapovalov[SU[2], {2}, <|{} -> 1|>, <|{} -> 1|>], 1, TestID->"shap-hw-norm" ];
VerificationTest[ ClassicalLieAlgebra`Representations`Private`shapovalov[SU[2], {2}, <|{1} -> 1|>, <|{1} -> 1|>], 2, TestID->"shap-f-norm" ];
VerificationTest[ ClassicalLieAlgebra`Representations`Private`shapovalov[SU[2], {2}, <|{1, 1} -> 1|>, <|{1, 1} -> 1|>], 4, TestID->"shap-ff-norm" ];
(* su(3) fundamental {1,0}: <f_1 v, f_1 v> = lambda_1 = 1 *)
VerificationTest[ ClassicalLieAlgebra`Representations`Private`shapovalov[SU[3], {1, 0}, <|{1} -> 1|>, <|{1} -> 1|>], 1, TestID->"shap-su3-f1" ];
(* orthogonality of different weights: <f_1 v, f_2 v> = 0 for su(3) *)
VerificationTest[ ClassicalLieAlgebra`Representations`Private`shapovalov[SU[3], {1, 0}, <|{1} -> 1|>, <|{2} -> 1|>], 0, TestID->"shap-orthog-weights" ];

(* Task 8: RepresentationMatrices via Shapovalov highest-weight construction *)
bracket[a_, b_] := a.b - b.a;
VerificationTest[ Module[{m = RepresentationMatrices[Irrep[SU[2], {2}]]},
   bracket[m["Raising"][[1]], m["Lowering"][[1]]] == m["Cartan"][[1]]], True, TestID->"su2-spin1-EF=H" ];
VerificationTest[ Sort[Diagonal[RepresentationMatrices[Irrep[SU[2], {2}]]["Cartan"][[1]]]], {-2, 0, 2}, TestID->"su2-spin1-H-spectrum" ];
VerificationTest[ Length[RepresentationMatrices[Irrep[SU[2], {2}]]["Cartan"][[1]]], 3, TestID->"su2-spin1-size3" ];
VerificationTest[ Module[{m = RepresentationMatrices[Irrep[SU[3], {1, 0}]]},
   And @@ Table[bracket[m["Raising"][[i]], m["Lowering"][[i]]] == m["Cartan"][[i]], {i, 2}]], True, TestID->"su3-fund-EF=H" ];
VerificationTest[ Length[RepresentationMatrices[Irrep[SU[3], {1, 0}]]["Cartan"][[1]]], 3, TestID->"su3-fund-size3" ];
VerificationTest[ Module[{m = RepresentationMatrices[Irrep[SU[3], {1, 0}]]},
   bracket[m["Raising"][[1]], m["Lowering"][[2]]] == 0 IdentityMatrix[3]], True, TestID->"su3-fund-EF-offdiag-0" ];

(* Task 9: full Chevalley relations on degenerate, spinor, and sp cases *)
chevChecks[g_, w_] := Module[{m = RepresentationMatrices[Irrep[g, w]], A = CartanMatrix[g], r = Rank[g], d, z},
   d = Length[m["Cartan"][[1]]]; z = 0 IdentityMatrix[d];
   (And @@ Flatten[Table[bracket[m["Raising"][[i]], m["Lowering"][[j]]] == If[i == j, m["Cartan"][[i]], z], {i, r}, {j, r}]]) &&
   (And @@ Flatten[Table[bracket[m["Cartan"][[i]], m["Raising"][[j]]] == A[[j, i]] m["Raising"][[j]], {i, r}, {j, r}]])];

VerificationTest[ Length[RepresentationMatrices[Irrep[SU[3], {1, 1}]]["Cartan"][[1]]], 8, TestID->"su3-adjoint-dim8" ];
VerificationTest[ chevChecks[SU[3], {1, 1}], True, TestID->"su3-adjoint-chevalley-relations" ];
VerificationTest[ Length[RepresentationMatrices[Irrep[SO[5], {0, 1}]]["Cartan"][[1]]], 4, TestID->"so5-spinor-dim4" ];
VerificationTest[ chevChecks[SO[5], {0, 1}], True, TestID->"so5-spinor-chevalley-relations" ];
VerificationTest[ chevChecks[Sp[4], {1, 0}], True, TestID->"sp4-defining-chevalley-relations" ];
VerificationTest[ chevChecks[SO[7], {0, 0, 1}], True, TestID->"so7-spinor-chevalley-relations" ];
(* Phase-1 tie-in: the fundamental irrep's H-spectra match the Chevalley defining rep's H diagonals *)
VerificationTest[
   Sort /@ (Diagonal /@ RepresentationMatrices[Irrep[SU[3], {1, 0}]]["Cartan"]) ===
   Sort /@ (Diagonal /@ Generators[SU[3], "Chevalley"]["Cartan"]),
   True, TestID->"su3-fund-matches-phase1-Hspectra" ];
