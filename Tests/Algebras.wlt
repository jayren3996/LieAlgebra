Needs["ClassicalLieAlgebra`"];

VerificationTest[ SU[3], LieAlgebra["A", 2], TestID -> "SU-sugar" ];
VerificationTest[ Sp[4], LieAlgebra["C", 2], TestID -> "Sp-sugar" ];
VerificationTest[ SO[5], LieAlgebra["B", 2], TestID -> "SO-odd-sugar" ];
VerificationTest[ SO[6], LieAlgebra["D", 3], TestID -> "SO-even-sugar" ];
VerificationTest[ Sp[3], $Failed, {Sp::evenrank}, TestID -> "Sp-odd-fails" ];
VerificationTest[ LieAlgebra["E", 2], $Failed, {LieAlgebra::badtype}, TestID -> "bad-type-fails" ];
VerificationTest[ LieAlgebra["A", 0], $Failed, {LieAlgebra::badrank}, TestID -> "bad-rank-fails" ];
VerificationTest[ MakeBoxes[LieAlgebra["A", 2], StandardForm], SubscriptBox["A", "2"], TestID -> "display-subscript" ];

VerificationTest[ Rank[SU[4]], 3, TestID -> "rank-A3" ];
VerificationTest[ LieAlgebraDimension[SU[3]], 8, TestID -> "dim-su3" ];
VerificationTest[ LieAlgebraDimension[SO[5]], 10, TestID -> "dim-so5" ];
VerificationTest[ LieAlgebraDimension[Sp[4]], 10, TestID -> "dim-sp4" ];
VerificationTest[ CartanMatrix[LieAlgebra["A", 2]], {{2, -1}, {-1, 2}}, TestID -> "cartan-A2" ];
VerificationTest[ CartanMatrix[LieAlgebra["B", 2]], {{2, -2}, {-1, 2}}, TestID -> "cartan-B2" ];
VerificationTest[ CartanMatrix[LieAlgebra["C", 2]], {{2, -1}, {-2, 2}}, TestID -> "cartan-C2" ];
VerificationTest[ CartanMatrix[LieAlgebra["D", 4]], {{2,-1,0,0},{-1,2,-1,-1},{0,-1,2,0},{0,-1,0,2}}, TestID -> "cartan-D4" ];
VerificationTest[ Length[PositiveRoots[LieAlgebra["A", 2]]], 3, TestID -> "posroots-A2-count" ];
VerificationTest[ FundamentalWeights[LieAlgebra["A", 2]], {{2/3, -1/3, -1/3}, {1/3, 1/3, -2/3}}, TestID -> "fundweights-A2" ];
