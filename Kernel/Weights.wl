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

(* form normalized to long-root^2 = 2: raw Dot is correct for A, B_{r>=2}, D;
   halve for C, and double B1 since its only root has raw length squared 1. *)
casimirForm[g_, u_, v_] := Which[
   g[[1]] === "C", (u . v)/2,
   g[[1]] === "B" && g[[2]] === 1, 2 (u . v),
   True, u . v];
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

(* ------------------------------------------------------------------ *)
(* Freudenthal multiplicity recursion                                  *)
(* ------------------------------------------------------------------ *)

(* freudenthalMults[g, lambda] returns Association[euclideanWeight -> multiplicity]
   for ALL weights in the irrep with highest weight lambda (Dynkin labels).
   Uses BFS layered by level = (number of simple roots subtracted from lam). *)

freudenthalMults[g_, lambda_] := Module[
  {lam, rho, pos, simple, r, mult, frontier, nextFrontier, nu, alpha,
   denom, accum, mu2, k, multVal},

  lam    = toEuclidean[g, lambda];
  rho    = weylVector[g];
  pos    = PositiveRoots[g];
  simple = SimpleRoots[g];
  r      = Length[simple];

  (* Start: highest weight has multiplicity 1 *)
  mult     = <| lam -> 1 |>;
  frontier = {lam};

  While[frontier =!= {},
    nextFrontier = {};

    (* For each weight mu in current frontier, try subtracting each simple root *)
    Do[
      Do[
        nu = mu - simple[[i]];

        (* Only process nu if not already computed *)
        If[! KeyExistsQ[mult, nu],

          (* Freudenthal denominator *)
          denom = (lam + rho) . (lam + rho) - (nu + rho) . (nu + rho);

          (* If denom == 0, nu is not in the weight system (same Casimir shell as lam
             but not the highest weight — can only happen when accum = 0 too, so
             multVal = 0; skip to avoid 0/0 messages) *)
          If[denom == 0, Continue[]];

          (* denom > 0 for genuine descendant weights *)
          accum = 0;
          Do[
            k   = 1;
            mu2 = nu + alpha;
            While[KeyExistsQ[mult, mu2],
              accum += mult[mu2] * (mu2 . alpha);
              k++;
              mu2 = nu + k * alpha
            ],
          {alpha, pos}];

          multVal = 2 * accum / denom;

          (* Integrality check: must be a positive integer *)
          If[multVal > 0,
            Assert[IntegerQ[multVal]];
            mult[nu] = multVal;
            AppendTo[nextFrontier, nu]
          ]
        ],
      {i, r}],
    {mu, frontier}];

    (* Deduplicate nextFrontier before next BFS level *)
    frontier = DeleteDuplicates[nextFrontier]
  ];

  mult
];

(* weightSystemDynkin: convert Euclidean-keyed association to Dynkin-label keys *)
weightSystemDynkin[g_, lambda_] :=
  KeyMap[toDynkin[g, #] &, freudenthalMults[g, lambda]];

End[];
EndPackage[];
