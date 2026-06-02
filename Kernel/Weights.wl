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

(* Weyl dimension formula: Prod_{alpha>0} (lambda+rho, alpha)/(rho, alpha), exact. *)
weylDim[g_, lambda_] := Module[{lam = toEuclidean[g, lambda], rho = weylVector[g], pos = PositiveRoots[g]},
   Times @@ Table[((lam + rho) . a)/(rho . a), {a, pos}]];

(* form normalized to long-root^2 = 2: raw Dot already correct for A,B,D; halve for C. *)
casimirForm[g_, u_, v_] := If[g[[1]] === "C", (u . v)/2, u . v];
casimir[g_, lambda_] := Module[{lam = toEuclidean[g, lambda], rho = weylVector[g]},
   casimirForm[g, lam, lam + 2 rho]];

End[];
EndPackage[];
