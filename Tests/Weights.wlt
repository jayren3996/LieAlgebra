Needs["ClassicalLieAlgebra`"];
VerificationTest[ ClassicalLieAlgebra`Weights`Private`toDynkin[SU[3], ClassicalLieAlgebra`Weights`Private`weylVector[SU[3]]], {1,1}, TestID->"rho-dynkin-A2" ];
VerificationTest[ ClassicalLieAlgebra`Weights`Private`toDynkin[SU[3], ClassicalLieAlgebra`Weights`Private`toEuclidean[SU[3], {2,1}]], {2,1}, TestID->"dynkin-euclid-roundtrip" ];
VerificationTest[ ClassicalLieAlgebra`Weights`Private`toDynkin[SO[5], ClassicalLieAlgebra`Weights`Private`toEuclidean[SO[5], {0,1}]], {0,1}, TestID->"roundtrip-B2-spinor" ];
