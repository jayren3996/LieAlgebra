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
