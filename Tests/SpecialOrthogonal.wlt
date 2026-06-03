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

(* Batch B: BasisTransform conjugates the antisymmetric Chevalley realization into the
   diagonal one, on H, E and F -- for both the B and D families. *)
conjOK[b_, a_, d_, k_] := AllTrue[Flatten[((b . # . Inverse[b]) & /@ a[k]) - d[k]], Simplify[#] === 0 &];
soConjTest[g_] := With[{b = BasisTransform[g],
     a = Generators[g, "Chevalley", "Realization" -> "Antisymmetric"],
     d = Generators[g, "Chevalley", "Realization" -> "Diagonal"]},
   conjOK[b, a, d, "Cartan"] && conjOK[b, a, d, "Raising"] && conjOK[b, a, d, "Lowering"]];
VerificationTest[ soConjTest[SO[5]], True, TestID -> "so5-basistransform-conjugates-to-diagonal" ];
VerificationTest[ soConjTest[SO[6]], True, TestID -> "so6-basistransform-conjugates-to-diagonal" ];
VerificationTest[ soConjTest[SO[7]], True, TestID -> "so7-basistransform-conjugates-to-diagonal" ];
VerificationTest[ With[{b = BasisTransform[SO[6]]}, Simplify[b . ConjugateTranspose[b]] == IdentityMatrix[6]], True, TestID -> "so6-basistransform-unitary" ];
