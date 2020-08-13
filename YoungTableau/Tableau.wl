(* ::Package:: *)

BeginPackage["Tableau`"];

(* Public functions *)
CycleDot
CycleProduct

Tableau
Symmetrizer
TableauForm
TableauTranspose

Psi
TensorDot
TensorNorm
TensorPermute
TableauPermute

ToTensor
TensorTableau

ColumnCanonicalize


Begin["`Private`"];

protectlist={
    "Tableau",
    "Psi",
    "TensorTableau",
    "Symmetrizer",
    "TableauForm",
    "ToTensor",
    "TensorDot",
    "TensorNorm"
};


(*Permutation*)
CycleDotAtom[c1_Cycles, c2_Cycles]:= PermutationProduct[c2, c1];
CycleDotAtom[a_*c1_Cycles, c2_Cycles]:= a * PermutationProduct[c2, c1];
CycleDotAtom[c1_Cycles, b_*c2_Cycles]:= b * PermutationProduct[c2, c1];
CycleDotAtom[a_*c1_Cycles, b_*c2_Cycles]:= a*b * PermutationProduct[c2, c1];

CycleDot[a_, b_]:= Distribute[CycleDotAtom[a,b]];
CycleProduct[a__]:= Fold[CycleDot, {a}];

(*Young tableax*)
TableauTranspose[t_List]:=Block[{l=Length[t], li=Length/@t, p=1, out={}, i, j},
    For[i=1, i<=l, i++,
        j = l-i+1;
        out = Join[out, Transpose@t[[1;;j, p;;li[[j]]]]];
        p = li[[j]]+1;
    ];
    out
];
ShowTableau[t_List]:= Grid[Map[Item[#,Frame->True]&,#]&/@t];
TableauTranspose[t_Tableau]:= Tableau @ TableauTranspose @ t[[1]];
TableauForm[t_Tableau]:= ShowTableau @ t[[1]];

(* Young symmetrizer *)
Parity[p_]:= (-1)^(Total[Length/@p[[1]]]+Length[p[[1]]]);
PermutationGroupElements[l_]:= GroupElements @ PermutationGroup[Cycles[{#}]&/@{l[[1;;2]],l}];
PermuteElements[l_,1]:= Total @ PermutationGroupElements[l];
PermuteElements[l_,2]:= (Dot[Parity/@#,#]&) @ PermutationGroupElements[l];
YoungSymmetrizer[t_List]:=Block[{p, q, tt=TableauTranspose[t], l, i},
    p = Cycles[{{}}];
    q = Cycles[{{}}];
    For[i=1, i<=Length[t], i++,
        If[Length[t[[i]]]<2, Break[]];
        p = CycleDot[PermuteElements[t[[i]],1], p];
    ];
    For[i=1, i<=Length[tt], i++,
        If[Length[tt[[i]]]<2, Break[]];
        q = CycleDot[PermuteElements[tt[[i]],2], q];
    ];
    CycleDot[p , q]
];
Symmetrizer[t_Tableau]:= YoungSymmetrizer @ t[[1]];

(* Tensor Permutation *)
TensorPermuteAtom[c_Cycles, v_Psi]:=Permute[v, c];
TensorPermuteAtom[a_*c_Cycles, v_Psi]:= a * Permute[v, c];
TensorPermuteAtom[c_Cycles, b_*v_Psi]:= b * Permute[v, c];
TensorPermuteAtom[a_*c_Cycles, b_*v_Psi]:= a*b * Permute[v, c];
TensorPermute[a_, b_Psi]:= Distribute @ TensorPermuteAtom[a, b];
TableauPermute[t_Tableau, v_Psi]:= TensorPermute[Symmetrizer[t], v];

(* Tensor Normalization *)
SquarePsi[p_Psi]:= 1;
SquarePsi[a_*p_Psi]:= Abs[a]^2;
TensorNorm[p_Psi]:= 1;
TensorNorm[a_*p_Psi]:= Abs[a];
TensorNorm[p_]:= Sqrt @ (SquarePsi /@ p)

DotPsi[p1_Psi, p2_Psi]:= If[SameQ[p1,p2], 1, 0];
DotPsi[a_*p1_Psi, p2_Psi]:= If[SameQ[p1,p2], Conjugate[a], 0];
DotPsi[p1_Psi, b_*p2_Psi]:= If[SameQ[p1,p2], b, 0];
DotPsi[a_*p1_Psi, b_*p2_Psi]:= If[SameQ[p1,p2], Conjugate[a]*b, 0];
TensorDot[p1_, p2_]:= Distribute @ DotPsi[p1, p2];

(* Tensor Tableau *)
TableauForm[t_TensorTableau]:= ShowTableau @ t[[1]];
TableauForm[a_*t_TensorTableau]:= a * ShowTableau @ t[[1]];
TableauForm[t1_+t2_]:= TableauForm[t1] + TableauForm[t2];

ListToTensor[t_List]:=Block[{l=Length[t], li=Length/@t, i, p=1, tab={}},
    For[i=1,i<=l,i++, 
        AppendTo[tab, Range[p,p+li[[i]]-1]];
        p+=li[[i]];
    ];
    TableauPermute[Tableau[tab], Psi @@ Flatten[t]]
];
TableauToTensor[t_TensorTableau]:= ListToTensor @ t[[1]];
TableauToTensor[a_*t_TensorTableau]:= a * ListToTensor @ t[[2]];

ToTensor[t_TensorTableau]:= ListToTensor @ t[[1]];
ToTensor[a_*t_TensorTableau]:= a * ListToTensor @ t[[2]];
ToTensor[t1_+t2_]:= ToTensor[t1] + ToTensor[t2];

(* Canonicalize Tensor Tableau *)
ColumnCanonicalize[tab_]:=Block[{t=tab, l, li, i, c, s, p=1},
    l = Length @ t;
    li = Length /@ t;
    l1 = li[[1]];
    cl = ConstantArray[0, l1];
    For[i=1,i<=l,i++, cl[[1;;li[[i]] ]] += 1];
    For[i=1,i<=l1, i++,
        c = t[[1;;cl[[i]], i]];
        s = Sort[c];
        t[[1;;cl[[i]], i]] = s;
        p *= Parity @ FindPermutation[c, s]
    ];
    {p, t}
];


(*Protection*)
Protect/@protectlist;
End[];

EndPackage[];
