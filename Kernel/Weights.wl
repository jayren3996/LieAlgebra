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

dynkinLabelOf[g_, mu_, i_] := With[{a = SimpleRoots[g][[i]]}, 2 (mu . a)/(a . a)];
weylReflect[g_, mu_, i_] := mu - dynkinLabelOf[g, mu, i] SimpleRoots[g][[i]];

(* BFS with descent-only dedup: from weight w, reflect at node i only when its
   i-th Dynkin label is positive; collect every distinct image. *)
weylOrbit[g_, mu0_] := Module[{r = Rank[g], orbit = {mu0}, frontier = {mu0}, next, c, nu},
   While[frontier =!= {},
     next = {};
     Do[ Do[ c = dynkinLabelOf[g, w, i];
             If[c > 0, nu = weylReflect[g, w, i];
                If[! MemberQ[orbit, nu], AppendTo[orbit, nu]; AppendTo[next, nu]]],
          {i, r}], {w, frontier}];
     frontier = next];
   orbit];

End[];
EndPackage[];
