Needs["ClassicalLieAlgebra`"];

VerificationTest[ MemberQ[$Packages, "ClassicalLieAlgebra`"], True, TestID -> "package-loads" ];
VerificationTest[ Head[ClassicalLieAlgebra`Common`matrixUnit[2, 1, 2]], List, TestID -> "common-helper-visible-internally" ];
