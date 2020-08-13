(* ::Package:: *)
BeginPackage["SpecialUnitary`"]

SUT

SUCWH
SUCWE
SUCWF
SUCW

SUH
SUE
SUF
SUn

Begin["Private`"];

(* Standard Basis *)
SUT[{n_,1},{i_, j_}]:=Block[{m=ConstantArray[0,{n,n}], a, b}, 
    {a, b} = If[i<j, {i,j}, {j,i}];
    m[[a,b]] += 1/2; 
    m[[b,a]] += 1/2;
    m
];
SUT[{n_,2},{i_, j_}]:=Block[{m=ConstantArray[0,{n,n}], a, b}, 
    {a, b} = If[i<j, {i,j}, {j,i}];
    m[[a,b]] -= I/2; 
    m[[b,a]] += I/2;
    m
];
SUT[{n_,3},{i_}]:=Block[{d=ConstantArray[0,n]}, 
    d[[1;;i-1]] = 1/Sqrt[2i*(i-1)];
    d[[i]] = -Sqrt[(i-1)/(2i)];
    DiagonalMatrix[d]
];
SUT[n_]:=Block[{list={}, i, j},
    For[j=2, j<=n, j++,
        For[i=1, i<j, i++,
            AppendTo[list, SUT[{n,1},{i,j}]];
            AppendTo[list, SUT[{n,2},{i,j}]];
        ];
        AppendTo[list, SUT[{n,3},{j}]];
    ];
    list
];

(* Cartan-Weyl Basis *)
SUCWH[n_, i_]:=Block[{d=ConstantArray[0,n]}, 
    j = n+1-i;
    d[[1;;j-1]] = 1 / Sqrt[j*(j-1)];
    d[[j]]  = -Sqrt[(j-1)/j];
    DiagonalMatrix[d]
];
SUCWE[n_, {i_, j_}]:=Block[{m=ConstantArray[0,{n,n}], a, b}, 
    {a, b} = If[i<j, {i,j}, {j,i}];
    m[[a,b]] = 1;
    m
];
SUCWF[n_, {i_, j_}]:= Transpose @ SUCWE[n, {i,j}];
SUCW[n_]:=Block[{Hl, El={}, Fl, i, j},
    Hl = Table[SUCWH[n, i], {i,n-1}];
    For[j=2, j<=n, j++,
        For[i=1, i<j, i++,
            AppendTo[El, SUCWE[n, {i,j}]];
        ];
    ];
    Fl = Transpose /@ El;
    {Hl,El,Fl}
];

(* Chevalley Basis *)
SUH[n_, i_]:= Block[{d=ConstantArray[0,n]}, 
    d[[i]] = 1;
    d[[i+1]] = -1;
    DiagonalMatrix[d]
];
SUE[n_, i_]:= SUCWE[n, {i, i+1}];
SUF[n_, i_]:= SUCWF[n, {i, i+1}];

SUn[n_]:=Block[{Hl,El,Fl},
    Hl = Table[SUH[n, i], {i,n-1}];
    El = Table[SUE[n, i], {i,n-1}];
    Fl = Transpose /@ El;
    {Hl, El, Fl}
];


End[];
EndPackage[];