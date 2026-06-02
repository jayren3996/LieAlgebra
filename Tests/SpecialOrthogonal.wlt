Needs["ClassicalLieAlgebra`"];

dim[g_] := Length[First[Generators[g]]];
bracket[a_, b_] := a.b - b.a;
closesQ[gens_] := Module[{sp = Flatten /@ gens},
  AllTrue[Tuples[Range[Length@gens], 2],
    MatrixRank[Append[sp, Flatten[bracket[gens[[#[[1]]]], gens[[#[[2]]]]]]]] == MatrixRank[sp] &]];

VerificationTest[ dim[SO[5]], 5, TestID -> "so5-matdim" ];
VerificationTest[ dim[SO[6]], 6, TestID -> "so6-matdim" ];
VerificationTest[ Length[Generators[SO[5]]], 10, TestID -> "so5-gen-count" ];
VerificationTest[ closesQ[Generators[SO[5]]], True, TestID -> "so5-closes" ];
VerificationTest[ DiagonalMatrixQ[Generators[SO[6], "Chevalley"]["Cartan"][[1]]], True, TestID -> "so6-diag-cartan-default" ];
VerificationTest[ DiagonalMatrixQ[Generators[SO[6], "Chevalley", "Realization" -> "Antisymmetric"]["Cartan"][[1]]], False, TestID -> "so6-antisym-not-diag" ];
VerificationTest[ Det[BasisTransform[SO[6]]] != 0, True, TestID -> "so6-basistransform-invertible" ];
VerificationTest[ Generators[SO[6], "Chevalley", "Realization" -> "Nope"], $Failed, {Generators::badrealization}, TestID -> "so-bad-realization" ];
