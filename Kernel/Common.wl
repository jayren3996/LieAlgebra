BeginPackage["ClassicalLieAlgebra`Common`"];

classicalTypeQ;   (* classicalTypeQ["A"] etc. *)

Begin["`Private`"];

classicalTypeQ[t_] := MemberQ[{"A", "B", "C", "D"}, t];

End[];
EndPackage[];
