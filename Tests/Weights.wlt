Needs["ClassicalLieAlgebra`"];
VerificationTest[ ClassicalLieAlgebra`Weights`Private`toDynkin[SU[3], ClassicalLieAlgebra`Weights`Private`weylVector[SU[3]]], {1,1}, TestID->"rho-dynkin-A2" ];
VerificationTest[ ClassicalLieAlgebra`Weights`Private`toDynkin[SU[3], ClassicalLieAlgebra`Weights`Private`toEuclidean[SU[3], {2,1}]], {2,1}, TestID->"dynkin-euclid-roundtrip" ];
VerificationTest[ ClassicalLieAlgebra`Weights`Private`toDynkin[SO[5], ClassicalLieAlgebra`Weights`Private`toEuclidean[SO[5], {0,1}]], {0,1}, TestID->"roundtrip-B2-spinor" ];
(* A2: Weyl group S3; orbit of the fundamental weight omega_1 has 3 elements *)
VerificationTest[ Length[ClassicalLieAlgebra`Weights`Private`weylOrbit[SU[3], ClassicalLieAlgebra`Weights`Private`toEuclidean[SU[3], {1,0}]]], 3, TestID->"orbit-A2-omega1" ];
(* B2: orbit of the vector-rep highest weight (1,0) = {+-e1,+-e2} has 4 elements *)
VerificationTest[ Length[ClassicalLieAlgebra`Weights`Private`weylOrbit[SO[5], ClassicalLieAlgebra`Weights`Private`toEuclidean[SO[5], {1,0}]]], 4, TestID->"orbit-B2-vector" ];
(* A2 adjoint dominant weight (1,1): its orbit is the 6 long roots *)
VerificationTest[ Length[ClassicalLieAlgebra`Weights`Private`weylOrbit[SU[3], ClassicalLieAlgebra`Weights`Private`toEuclidean[SU[3], {1,1}]]], 6, TestID->"orbit-A2-omega1+omega2" ];
(* the zero weight is a fixed point: orbit size 1 *)
VerificationTest[ Length[ClassicalLieAlgebra`Weights`Private`weylOrbit[SU[3], {0,0,0}]], 1, TestID->"orbit-zero" ];
