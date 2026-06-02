BeginPackage["ClassicalLieAlgebra`Algebras`"];
Needs["ClassicalLieAlgebra`"];
Needs["ClassicalLieAlgebra`Common`"];
Begin["`Private`"];

validRankQ["A", r_] := IntegerQ[r] && r >= 1;
validRankQ["B", r_] := IntegerQ[r] && r >= 1;
validRankQ["C", r_] := IntegerQ[r] && r >= 1;
validRankQ["D", r_] := IntegerQ[r] && r >= 2;
validRankQ[_, _] := False;

LieAlgebra[t_, r_] /; (! classicalTypeQ[t]) := (Message[LieAlgebra::badtype, t]; $Failed);
LieAlgebra[t_, r_] /; (classicalTypeQ[t] && ! validRankQ[t, r]) := (Message[LieAlgebra::badrank, r, t]; $Failed);

SU[n_Integer] /; n >= 2 := LieAlgebra["A", n - 1];
SU[n_] := (Message[SU::baddim, n]; $Failed);
SO[n_Integer] /; (n >= 3 && OddQ[n]) := LieAlgebra["B", (n - 1)/2];
SO[n_Integer] /; (n >= 4 && EvenQ[n]) := LieAlgebra["D", n/2];
SO[n_] := (Message[SO::baddim, n]; $Failed);
Sp[n_Integer] /; (n >= 2 && EvenQ[n]) := LieAlgebra["C", n/2];
Sp[n_] := (Message[Sp::evenrank, n]; $Failed);

LieAlgebra /: MakeBoxes[LieAlgebra[t_String, r_Integer], StandardForm] := SubscriptBox[ToString[t], ToString[r]];

Rank[HoldPattern[LieAlgebra[_, r_]]] := r;

LieAlgebraDimension[HoldPattern[LieAlgebra["A", r_]]] := r (r + 2);
LieAlgebraDimension[HoldPattern[LieAlgebra["B", r_]]] := r (2 r + 1);
LieAlgebraDimension[HoldPattern[LieAlgebra["C", r_]]] := r (2 r + 1);
LieAlgebraDimension[HoldPattern[LieAlgebra["D", r_]]] := r (2 r - 1);

ee[d_, i_] := UnitVector[d, i];
SimpleRoots[HoldPattern[LieAlgebra["A", r_]]] := Table[ee[r + 1, i] - ee[r + 1, i + 1], {i, r}];
SimpleRoots[HoldPattern[LieAlgebra["B", r_]]] := Append[Table[ee[r, i] - ee[r, i + 1], {i, r - 1}], ee[r, r]];
SimpleRoots[HoldPattern[LieAlgebra["C", r_]]] := Append[Table[ee[r, i] - ee[r, i + 1], {i, r - 1}], 2 ee[r, r]];
SimpleRoots[HoldPattern[LieAlgebra["D", r_]]] := Append[Table[ee[r, i] - ee[r, i + 1], {i, r - 1}], ee[r, r - 1] + ee[r, r]];

CartanMatrix[g_LieAlgebra] := Module[{a = SimpleRoots[g]},
  Table[2 (a[[i]] . a[[j]])/(a[[j]] . a[[j]]), {i, Length@a}, {j, Length@a}]];

FundamentalWeights[g_LieAlgebra] := Inverse[CartanMatrix[g]] . SimpleRoots[g];

PositiveRoots[HoldPattern[LieAlgebra["A", r_]]] := Module[{d = r + 1},
  Flatten[Table[ee[d, i] - ee[d, j], {i, d}, {j, i + 1, d}], 1]];
PositiveRoots[HoldPattern[LieAlgebra["B", r_]]] := Join[
  Flatten[Table[ee[r, i] - ee[r, j], {i, r}, {j, i + 1, r}], 1],
  Flatten[Table[ee[r, i] + ee[r, j], {i, r}, {j, i + 1, r}], 1],
  Table[ee[r, i], {i, r}]];
PositiveRoots[HoldPattern[LieAlgebra["C", r_]]] := Join[
  Flatten[Table[ee[r, i] - ee[r, j], {i, r}, {j, i + 1, r}], 1],
  Flatten[Table[ee[r, i] + ee[r, j], {i, r}, {j, i + 1, r}], 1],
  Table[2 ee[r, i], {i, r}]];
PositiveRoots[HoldPattern[LieAlgebra["D", r_]]] := Join[
  Flatten[Table[ee[r, i] - ee[r, j], {i, r}, {j, i + 1, r}], 1],
  Flatten[Table[ee[r, i] + ee[r, j], {i, r}, {j, i + 1, r}], 1]];

End[];
EndPackage[];
