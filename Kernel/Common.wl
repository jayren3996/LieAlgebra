BeginPackage["ClassicalLieAlgebra`Common`"];

matrixUnit;       (* matrixUnit[n,i,j] = n x n matrix with 1 at (i,j) *)
classicalTypeQ;   (* classicalTypeQ["A"] etc. *)

Begin["`Private`"];

matrixUnit[n_Integer, i_Integer, j_Integer] := Normal[SparseArray[{{i, j} -> 1}, {n, n}]];
classicalTypeQ[t_] := MemberQ[{"A", "B", "C", "D"}, t];

End[];
EndPackage[];
