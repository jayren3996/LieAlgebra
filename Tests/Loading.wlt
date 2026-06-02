Needs["ClassicalLieAlgebra`"];

VerificationTest[ MemberQ[$Packages, "ClassicalLieAlgebra`"], True, TestID -> "package-loads" ];
VerificationTest[ Head[ClassicalLieAlgebra`Common`matrixUnit[2, 1, 2]], List, TestID -> "common-helper-visible-internally" ];

VerificationTest[ MemberQ[Attributes[Generators], Protected], True, TestID -> "generators-protected" ];
VerificationTest[ CartanWeyl[SU[3]], Generators[SU[3], "CartanWeyl"], TestID -> "cartanweyl-alias" ];
VerificationTest[ Chevalley[SU[3]], Generators[SU[3], "Chevalley"], TestID -> "chevalley-alias" ];
VerificationTest[ StringQ[Generators::usage], True, TestID -> "generators-has-usage" ];
