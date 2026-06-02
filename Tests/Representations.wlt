Needs["ClassicalLieAlgebra`"];
VerificationTest[ HighestWeight[Irrep[SU[3], {1, 1}]], {1, 1}, TestID->"hw" ];
VerificationTest[ Head[Irrep[SU[3], {1, 1}]], Irrep, TestID->"irrep-inert" ];
VerificationTest[ Irrep[SU[3], {1, 1, 0}], $Failed, {Irrep::badweight}, TestID->"wrong-length" ];
VerificationTest[ Irrep[SU[3], {1, -1}], $Failed, {Irrep::badweight}, TestID->"negative-label" ];
VerificationTest[ Irrep[SO[5], {0, 1}], Irrep[LieAlgebra["B", 2], {0, 1}], TestID->"sugar-canonicalizes" ];
