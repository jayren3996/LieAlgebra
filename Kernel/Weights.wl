BeginPackage["ClassicalLieAlgebra`Weights`"];
Needs["ClassicalLieAlgebra`"];
Needs["ClassicalLieAlgebra`Common`"];
Begin["`Private`"];

(* lambda (Dynkin labels) -> Euclidean coordinates: sum_i a_i * omega_i *)
toEuclidean[g_, a_List] := a . FundamentalWeights[g];

(* Euclidean mu -> Dynkin labels: a_i = 2 (mu . alpha_i)/(alpha_i . alpha_i) *)
toDynkin[g_, mu_] := With[{sr = SimpleRoots[g]},
   Table[2 (mu . sr[[i]])/(sr[[i]] . sr[[i]]), {i, Length[sr]}]];

weylVector[g_] := Total[FundamentalWeights[g]];   (* rho = sum of fundamental weights *)

End[];
EndPackage[];
