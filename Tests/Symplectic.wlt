Needs["ClassicalLieAlgebra`"];

bracket[a_, b_] := a.b - b.a;
closesQ[gens_] := Module[{sp = Flatten /@ gens},
  AllTrue[Tuples[Range[Length@gens], 2],
    MatrixRank[Append[sp, Flatten[bracket[gens[[#[[1]]]], gens[[#[[2]]]]]]]] == MatrixRank[sp] &]];

VerificationTest[ Length[First[Generators[Sp[4]]]], 4, TestID -> "sp4-matdim" ];
VerificationTest[ Length[Generators[Sp[4]]], 10, TestID -> "sp4-gen-count" ];
VerificationTest[ closesQ[Generators[Sp[4]]], True, TestID -> "sp4-closes" ];

(* regression: the corrected standard basis is a genuine compact-symplectic usp basis *)
VerificationTest[ AllTrue[Generators[Sp[4]], # == -ConjugateTranspose[#] &], True, TestID -> "sp4-anti-hermitian" ];
VerificationTest[
  With[{g = Generators[Sp[4]], m = 4},
    Module[{jmat = Array[jj, {m, m}], vars, eqs, lin, forms},
      vars = Flatten[jmat];
      eqs = Flatten[Table[Transpose[x] . jmat + jmat . x, {x, g}]];
      lin = Normal@CoefficientArrays[eqs, vars][[2]];
      forms = NullSpace[lin];
      Length[forms] == 1 && With[{jj0 = ArrayReshape[forms[[1]], {m, m}]}, jj0 == -Transpose[jj0] && Det[jj0] =!= 0]]],
  True, TestID -> "sp4-symplectic-form" ];
VerificationTest[ Det[BasisTransform[Sp[4]]] != 0, True, TestID -> "sp4-basistransform-invertible" ];
VerificationTest[ DiagonalMatrixQ[Generators[Sp[6], "Chevalley"]["Cartan"][[1]]], True, TestID -> "sp6-diag-cartan" ];
VerificationTest[ Generators[Sp[4], "Chevalley", "Realization" -> "Nope"], $Failed, {Generators::badrealization}, TestID -> "sp-bad-realization" ];
