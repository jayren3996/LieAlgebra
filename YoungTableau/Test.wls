#!/usr/bin/env wolframscript
(* ::Package:: *)

Import["Tableau.wl"];
PASS = 0;
FAIL = 0;
Test[a_, b_]:= If[SameQ[a,b], 
    PASS += 1,
    Print["Error: ", a, " \[NotEqual] ", b];
    FAIL += 1
];
Result[]:=Print["Test Ended: ", PASS," passed, ", FAIL, " failed."];

(* Test Cycles *)
Module[{a,b,c,d, x,y,z},
    a = Cycles[{{1,2}}];
    b = Cycles[{{2,3}}];
    c = Cycles[{{2,3,4}}];
    d = Cycles[{{1,2},{3,4}}];

    x = 2a + 3b;
    y = b + c;
    z = -c + d;

    Test[
        CycleDot[ a,  b],
        1*Cycles[{{1,2,3}}] 
    ];
    Test[
        CycleDot[2a,  d],  
        2*Cycles[{{3,4}}] 
    ];
    Test[ 
        CycleDot[-a, 3b], 
        -3*Cycles[{{1,2,3}}] 
    ];
    Test[ 
        CycleDot[ a, 2b],  
        2*Cycles[{{1,2,3}}] 
    ];

    Test[
        CycleDot[a,x], 
        2*Cycles[{{}}] + 3*Cycles[{{1,2,3}}]
    ];
    Test[
        CycleDot[x,a], 
        2*Cycles[{{}}] + 3*Cycles[{{3,2,1}}]
    ];
    Test[
        CycleDot[x,y], 
        2*Cycles[{{1,2,3}}] + 2*Cycles[{{1,2,3,4}}] + 3*Cycles[{{}}] + 3*Cycles[{{3,4}}]
    ];
    Test[
        CycleProduct[x,y,z], 
        out = 0;
        out -= 2*Cycles[{{1,2},{3,4}}];
        out -= 2*Cycles[{{1,2,4,3}}];
        out -= 3*Cycles[{{2,3,4}}];
        out -= 3*Cycles[{{4,2}}];
        out += 2*Cycles[{{1,3,4}}]; 
        out += 2*Cycles[{{1,3}}]; 
        out += 3*Cycles[{{1,2},{3,4}}]; 
        out += 3*Cycles[{{1,2}}]; 
        out
    ];
];

(* Test Young Tableau *)
Module[{sym},
    Test[
        TableauTranspose @ Tableau[{{1,2},{3}}],
        Tableau[{{1,3},{2}}]
    ];
    Test[
        TableauTranspose @ Tableau[{{1,1,1},{2,2,2},{3},{4}}],
        Tableau[{{1,2,3,4},{1,2},{1,2}}]
    ];
    Test[
        Symmetrizer @ Tableau[{{1,2},{3}}],
        Cycles[{{}}]+Cycles[{{1,2}}]-Cycles[{{1,3}}]-Cycles[{{2,1,3}}]
    ];
    Test[
        Symmetrizer @ Tableau[{{1,2,3},{4}}],
        sym = Cycles[{{}}];
        sym += Cycles[{{1,2}}];
        sym += Cycles[{{1,3}}];
        sym += Cycles[{{2,3}}];
        sym += Cycles[{{1,2,3}}];
        sym += Cycles[{{3,2,1}}];
        sym -= Cycles[{{1,4}}];
        sym -= Cycles[{{2,1,4}}];
        sym -= Cycles[{{3,1,4}}];
        sym -= Cycles[{{2,3},{1,4}}];
        sym -= Cycles[{{2,3,1,4}}];
        sym -= Cycles[{{3,2,1,4}}];
        sym
    ];
];

(* Test Tensor Permute *)
Module[{v,t},
    Test[
        t = Tableau[{{1,2,3}}];
        v = Psi[1,1,1];
        TableauPermute[t,v],
        6*v
    ];
    Test[
        t = Tableau[{{1,2},{3}}];
        v = Psi[1,1,2];
        TableauPermute[t,v],
        2*Psi[1,1,2] - Psi[2,1,1] - Psi[1,2,1]
    ];
];

(* Test Tensor Norm *)
Module[{v,v1,v2},
    Test[
        v = Psi[1,2,3];
        TensorNorm[v],
        1
    ];
    Test[
        v = -Psi[1,2,3];
        TensorNorm[v],
        1
    ];
    Test[
        v = 3*Psi[1,2,3];
        TensorNorm[v],
        3
    ];
    Test[
        v = Psi[1,2,3] - 2*Psi[3,2,1] + Psi[4,5,6] - 2*Psi[3,2,1];
        TensorNorm[v],
        Sqrt[18]
    ];
    Test[
        v = Psi[1,2,3] - 2*Psi[3,2,1] + Psi[4,5,6] - 2*Psi[3,2,1];
        TensorDot[v,v],
        18
    ];
    Test[
        v1 = Psi[1,2,3] - 2*Psi[3,2,1] + Psi[4,5,6];
        v2 = 3*Psi[1,2,3] - Psi[4,2,1] + 6*Psi[3,2,1] - 2*Psi[3,8,1];
        TensorDot[v1,v2],
        -9
    ];
];

(* Test Tensor Tableau *)
Module[{t,v},
    Test[
        t = TensorTableau[{{1,1},{2}}];
        ToTensor[t],
        2*Psi[1,1,2] - Psi[2,1,1] - Psi[1,2,1]
    ];
    Test[
        t = TensorTableau[{{1,2,3},{4}}];
        ToTensor[t],
        v = Psi[1,2,3,4];
        v += Psi[2,1,3,4];
        v += Psi[3,2,1,4];
        v += Psi[1,3,2,4];
        v += Psi[3,1,2,4];
        v += Psi[2,3,1,4];
        v -= Psi[4,2,3,1];
        v -= Psi[2,4,3,1];
        v -= Psi[4,3,2,1];
        v -= Psi[3,4,2,1];
        v -= Psi[2,3,4,1];
        v -= Psi[3,2,4,1];
        v
    ];
];

(* Test Canonical *)
Module[{t,p,tp},
    t = {{4,3,1},{2,2},{3,1},{1}};
    {p, tp} = ColumnCanonicalize[t];
    Test[tp, {{1,1,1},{2,2},{3,3},{4}}];
    Test[p, 1];
];

Result[];