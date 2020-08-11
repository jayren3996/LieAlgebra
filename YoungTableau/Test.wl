(* ::Package:: *)

Import[NotebookDirectory[]<>"Tableau.wl"]


Test[a_,b_]:=If[a==b,Print["pass"],Print["Failed: ",a," \[NotEqual] ",b]];


(*Test Cycles*)
c0=Cycles[{{}}];
c1=Cycles[{{1,2}}];
c2=Cycles[[{{1,2,3}}]];
Test[c1*c2,Cycles[{{2,3}}]]




