BeginPackage["ClassicalLieAlgebra`Representations`"];
Needs["ClassicalLieAlgebra`"];
Needs["ClassicalLieAlgebra`Common`"];
Needs["ClassicalLieAlgebra`Weights`"];
Begin["`Private`"];

validWeightQ[g_, w_] := ListQ[w] && Length[w] === Rank[g] && AllTrue[w, IntegerQ[#] && # >= 0 &];

(* fire only on INVALID input; valid Irrep[g,w] stays inert. Guard on g_LieAlgebra so the
   literal-LieAlgebra-pattern validation problem (Phase 1) does not bite. *)
Irrep[g_LieAlgebra, w_] /; ! validWeightQ[g, w] := (Message[Irrep::badweight, w, g, Rank[g]]; $Failed);

HighestWeight[Irrep[g_LieAlgebra, w_]] := w;

RepresentationDimension[Irrep[g_LieAlgebra, w_]] := ClassicalLieAlgebra`Weights`Private`weylDim[g, w];

End[];
EndPackage[];
