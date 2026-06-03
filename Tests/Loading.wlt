Needs["ClassicalLieAlgebra`"];

VerificationTest[ MemberQ[$Packages, "ClassicalLieAlgebra`"], True, TestID -> "package-loads" ];
VerificationTest[ ClassicalLieAlgebra`Common`classicalTypeQ["A"], True, TestID -> "common-helper-visible-internally" ];

VerificationTest[ MemberQ[Attributes[Generators], Protected], True, TestID -> "generators-protected" ];
VerificationTest[ CartanWeyl[SU[3]], Generators[SU[3], "CartanWeyl"], TestID -> "cartanweyl-alias" ];
VerificationTest[ Chevalley[SU[3]], Generators[SU[3], "Chevalley"], TestID -> "chevalley-alias" ];
VerificationTest[ StringQ[Generators::usage], True, TestID -> "generators-has-usage" ];

(* Batch B: every public symbol is Protected after load (not just Generators) *)
VerificationTest[ AllTrue[Names["ClassicalLieAlgebra`*"], MemberQ[Attributes[#], Protected] &], True, TestID -> "all-public-symbols-protected" ];
