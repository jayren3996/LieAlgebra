Needs["ClassicalLieAlgebra`"];

VerificationTest[ SU[3], LieAlgebra["A", 2], TestID -> "SU-sugar" ];
VerificationTest[ Sp[4], LieAlgebra["C", 2], TestID -> "Sp-sugar" ];
VerificationTest[ SO[5], LieAlgebra["B", 2], TestID -> "SO-odd-sugar" ];
VerificationTest[ SO[6], LieAlgebra["D", 3], TestID -> "SO-even-sugar" ];
VerificationTest[ Sp[3], $Failed, {Sp::evenrank}, TestID -> "Sp-odd-fails" ];
VerificationTest[ LieAlgebra["E", 2], $Failed, {LieAlgebra::badtype}, TestID -> "bad-type-fails" ];
VerificationTest[ LieAlgebra["A", 0], $Failed, {LieAlgebra::badrank}, TestID -> "bad-rank-fails" ];
VerificationTest[ MakeBoxes[LieAlgebra["A", 2], StandardForm], SubscriptBox["A", "2"], TestID -> "display-subscript" ];
