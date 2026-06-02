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

End[];
EndPackage[];
