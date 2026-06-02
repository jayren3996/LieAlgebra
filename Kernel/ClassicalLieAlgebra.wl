BeginPackage["ClassicalLieAlgebra`"];

(* Public symbols + usage are filled in across later tasks; declared here. *)

LieAlgebra::usage = "LieAlgebra[type, rank] is the simple Lie algebra of Cartan type \"A\"|\"B\"|\"C\"|\"D\" and given rank.";
SU::usage = "SU[n] represents su(n), the type A_{n-1} algebra (n>=2).";
SO::usage = "SO[n] represents so(n): type B for odd n>=3, type D for even n>=4.";
Sp::usage = "Sp[n] represents sp(n), the type C_{n/2} algebra (even n>=2).";

LieAlgebra::badtype = "`1` is not a valid Cartan type; use \"A\", \"B\", \"C\" or \"D\".";
LieAlgebra::badrank = "`1` is not a valid rank for type `2`.";
SU::baddim = "SU[`1`] requires an integer n>=2.";
SO::baddim = "SO[`1`] requires an integer n>=3.";
Sp::evenrank = "Sp[`1`] requires an even integer n>=2.";

Begin["`Private`"];
Needs["ClassicalLieAlgebra`Common`"];
Get["ClassicalLieAlgebra`Algebras`"];
(* << subfiles added in later tasks *)
End[];

EndPackage[];
