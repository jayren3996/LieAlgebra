Needs["ClassicalLieAlgebra`"];

VerificationTest[ Length[Generators[SU[2]]], 3, TestID -> "su2-gen-count" ];
VerificationTest[ AllTrue[Generators[SU[2]], # == ConjugateTranspose[#] &], True, TestID -> "su2-hermitian" ];
VerificationTest[ AllTrue[Generators[SU[2]], Tr[#] == 0 &], True, TestID -> "su2-traceless" ];
VerificationTest[ Length[Generators[SU[3]]], 8, TestID -> "su3-gen-count" ];
VerificationTest[ Keys[Generators[SU[3], "Chevalley"]], {"Cartan", "Raising", "Lowering"}, TestID -> "cheval-assoc-keys" ];
VerificationTest[ Diagonal[Generators[SU[3], "Chevalley"]["Cartan"][[1]]], {1, -1, 0}, TestID -> "su3-chevalley-H1" ];
VerificationTest[ Generators[SU[3], "Chevalley"]["Lowering"][[1]], Transpose[Generators[SU[3], "Chevalley"]["Raising"][[1]]], TestID -> "su3-F-is-Etranspose" ];
VerificationTest[ Length[Generators[SU[3], "CartanWeyl"]["Raising"]], 3, TestID -> "su3-cw-raising-count" ];
VerificationTest[ Generators[SU[3], "Nope"], $Failed, {Generators::badscheme}, TestID -> "bad-scheme" ];
