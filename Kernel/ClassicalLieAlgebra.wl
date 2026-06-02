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

Irrep::usage = "Irrep[g, w] represents the irreducible representation of the classical Lie algebra g with highest weight given by the Dynkin labels w (a list of Rank[g] non-negative integers).";
HighestWeight::usage = "HighestWeight[Irrep[g, w]] returns the highest weight w (Dynkin labels).";
Irrep::badweight = "`1` is not a valid highest weight for `2`; expected a list of `3` non-negative integers.";

Generators::usage = "Generators[g] gives the defining-representation generators of g. Generators[g,\"CartanWeyl\"] and Generators[g,\"Chevalley\"] give associations <|\"Cartan\"->..,\"Raising\"->..,\"Lowering\"->..|>. Option \"Realization\"->\"Diagonal\"|\"Antisymmetric\" applies to so/sp.";

CartanWeyl::usage = "CartanWeyl[g] is shorthand for Generators[g, \"CartanWeyl\"].";
Chevalley::usage = "Chevalley[g] is shorthand for Generators[g, \"Chevalley\"].";

Tableau::usage = "Tableau[rows] represents a Young tableau given as a list of row-lists.";
TensorTableau::usage = "TensorTableau[rows] represents a tensor-product Young tableau with explicit index entries.";
Psi::usage = "Psi[i1,i2,...] represents a basis tensor |i1,i2,...> in the tensor product space.";
TableauForm::usage = "TableauForm[t] displays a TensorTableau (or linear combination thereof) as a grid.";
ToTensor::usage = "ToTensor[t] converts a TensorTableau (or linear combination) to a linear combination of Psi basis tensors.";
TableauPermute::usage = "TableauPermute[t,v] applies the Young symmetrizer of Tableau t to a Psi tensor v.";
TableauDot::usage = "TableauDot[t1,t2] computes the inner product of two TensorTableau expressions.";
TensorDot::usage = "TensorDot[p1,p2] computes the Hermitian inner product of two linear combinations of Psi tensors.";
TensorNorm::usage = "TensorNorm[p] gives the norm of a linear combination of Psi tensors.";
TableauNormalization::usage = "TableauNormalization[t] normalizes a TensorTableau expression to unit norm.";
TableauOrthogonalization::usage = "TableauOrthogonalization[t1,t2] returns {t1, t2 - <t2,t1>/<t1,t1> t1}, the Gram-Schmidt step.";
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
Get["ClassicalLieAlgebra`Symplectic`"];
Get["ClassicalLieAlgebra`YoungTableaux`"];
Get["ClassicalLieAlgebra`Weights`"];
Get["ClassicalLieAlgebra`Representations`"];
Generators[HoldPattern[_LieAlgebra], s_, OptionsPattern[]] := (Message[Generators::badscheme, s]; $Failed);

CartanWeyl[g_] := Generators[g, "CartanWeyl"];
Chevalley[g_] := Generators[g, "Chevalley"];

End[];

Protect[Evaluate[Names["ClassicalLieAlgebra`*"]]];

EndPackage[];
