BeginPackage["ClassicalLieAlgebra`"];

(* Public symbols + usage are filled in across later tasks; declared here. *)

LieAlgebra::usage = "LieAlgebra[type, rank] is the simple Lie algebra of Cartan type \"A\"|\"B\"|\"C\"|\"D\" and given rank.";
SU::usage = "SU[n] represents su(n), the type A_{n-1} algebra (n>=2).";
SO::usage = "SO[n] represents so(n): type B for odd n>=3, type D for even n>=4.";
Sp::usage = "Sp[n] represents sp(n), the type C_{n/2} algebra (even n>=2).";

Rank::usage = "Rank[g] gives the rank of the Lie algebra g.";
LieAlgebraDimension::usage = "LieAlgebraDimension[g] gives the dimension of the Lie algebra g.";
CartanMatrix::usage = "CartanMatrix[g] gives the Cartan matrix of g.";
SimpleRoots::usage = "SimpleRoots[g] gives the simple roots of g in the Euclidean basis.";
PositiveRoots::usage = "PositiveRoots[g] gives the positive roots of g in the Euclidean basis.";
FundamentalWeights::usage = "FundamentalWeights[g] gives the fundamental weights of g in the Euclidean basis.";

LieAlgebra::badtype = "`1` is not a valid Cartan type; use \"A\", \"B\", \"C\" or \"D\".";
LieAlgebra::badrank = "`1` is not a valid rank for type `2`.";
SU::baddim = "SU[`1`] requires an integer n>=2.";
SO::baddim = "SO[`1`] requires an integer n>=3.";
Sp::evenrank = "Sp[`1`] requires an even integer n>=2.";

Generators::usage = "Generators[g] gives the defining-representation generators of g. Generators[g,\"CartanWeyl\"] and Generators[g,\"Chevalley\"] give associations <|\"Cartan\"->..,\"Raising\"->..,\"Lowering\"->..|>. Option \"Realization\"->\"Diagonal\"|\"Antisymmetric\" applies to so/sp.";
BasisTransform::usage = "BasisTransform[g] gives the matrix conjugating the antisymmetric realization of so/sp into the diagonal one (identity for su).";
Generators::badscheme = "`1` is not a valid scheme; use \"Standard\", \"CartanWeyl\" or \"Chevalley\".";
Generators::badrealization = "`1` is not a valid \"Realization\"; use \"Diagonal\" or \"Antisymmetric\".";
Options[Generators] = {"Realization" -> "Diagonal"};
SyntaxInformation[Generators] = {"ArgumentsPattern" -> {_, _., OptionsPattern[]}};

Begin["`Private`"];
Needs["ClassicalLieAlgebra`Common`"];
Get["ClassicalLieAlgebra`Algebras`"];
Get["ClassicalLieAlgebra`SpecialUnitary`"];
Get["ClassicalLieAlgebra`SpecialOrthogonal`"];
(* << subfiles added in later tasks *)
Generators[HoldPattern[_LieAlgebra], s_, OptionsPattern[]] := (Message[Generators::badscheme, s]; $Failed);
End[];

EndPackage[];
