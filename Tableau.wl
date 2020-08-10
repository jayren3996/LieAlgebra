(* ::Package:: *)

BeginPackage["Tableau`"];
Quiet@Needs["Combinatorica`"];
(*Public functions*)
Times
Tableau
YoungSymmetrizer
Transpose
TableauForm
TensorTableau
Begin["`Private`"];


(*Remove conflict functions in Combinatorica`*)
namelist={
   "Combinatorica`Cycles",
   "Combinatorica`Permute",
   "Combinatorica`TableauQ"
};
functionlist={
   "Times",
   "Transpose"
};
Unprotect/@namelist;
Remove/@namelist;
(*Unprotect function names*)
Unprotect/@functionlist;


(*Permutation*)
CycleQ[a_]:=SameQ[Head[a],Cycles];
NCycleQ[a_]:=SameQ[Head[a],Times]&&SameQ[Head[a[[2]]],Cycles];
AbstractCycleQ[a_]:=CycleQ[a]||NCycleQ[a];
PermQ[a_]:=SameQ[Head[a],Plus]&&AbstractCycleQ[a[[1]]];
AbstractPermQ[a_]:=PermQ[a]||AbstractCycleQ[a];
Times[a_?CycleQ,b_?CycleQ]:=b\[PermutationProduct]a;
Times[a_?NCycleQ,b_?CycleQ]:=a[[1]]*b\[PermutationProduct]a[[2]];
Times[a_?CycleQ,b_?NCycleQ]:=b[[1]]*b[[2]]\[PermutationProduct]a;
Times[a_?NCycleQ,b_?NCycleQ]:=a[[1]]*b[[1]]*b[[2]]\[PermutationProduct]a[[2]];
Times[a_?PermQ,b_?AbstractPermQ]:=Module[{f},Distribute[f[a,b]]/.{f->Times}];
Times[a_?AbstractPermQ,b_?PermQ]:=Module[{f},Distribute[f[a,b]]/.{f->Times}];


(*Young tableax*)
TableauQ[t_]:=SameQ[Head[t],Tableau];
Transpose[t_?TableauQ]:=Tableau@TransposeTableau[t[[1]]];
TableauForm[x_?TableauQ]:=Grid[Map[Item[#,Frame->True]&,#]&/@x[[1]]];


(*Young symmetrizer*)
Parity[p_]:=(-1)^(Total[Length/@p[[1]]]+Length[p[[1]]]);
PermutationGroupElements[l_]:=GroupElements@PermutationGroup[Cycles[{#}]&/@{l[[1;;2]],l}];
PermuteElements[l_,1]:=Total@PermutationGroupElements[l];
PermuteElements[l_,2]:=(Dot[Parity/@#,#]&)@PermutationGroupElements[l];
YoungSymmetrizer[tab_?TableauQ]:=Module[{t=tab[[1]],res,l,tt,i},
   res=Cycles[{{}}];
   For[i=1,i<=Length[t],i++,
      If[Length[t[[i]]]<2,Break[]];
      res=res*PermuteElements[t[[i]],1];
   ];
   tt=TransposeTableau[t];
   For[i=1,i<=Length[tt],i++,
      If[Length[tt[[i]]]<2,Break[]];
      res=res*PermuteElements[tt[[i]],2];
   ];
   res
];
TensorTableau[tab_?TableauQ,ind_]:=Module[{sym,len,states,coeffs,s,c,p,i},
   sym=YoungSymmetrizer[tab];
   len=Length[sym];
   coeffs={};
   states={};
   For[i=1,i<=len,i++,
      s=sym[[i]];
      {c,p}=If[CycleQ[s],
         {1,Permute[ind,s]},
         {s[[1]],Permute[ind,s[[2]]]}
      ];
      If[MemberQ[states,p],
         coeffs[[FirstPosition[states,p]]]+=c,
         AppendTo[states,p];
         AppendTo[coeffs,c];
      ];
   ];
   {coeffs,states}
];


(*Protection*)
Protect/@functionlist;
Protect[Tableau];
End[];

EndPackage[];
