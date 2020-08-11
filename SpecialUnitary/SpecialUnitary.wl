(* ::Package:: *)
BeginPackage["SpecialUnitary`"]

SU::usage="SU[N, string, input] return properties of SU(N)";

Begin["Private`"];
(* Helper FUnctiond *)
Commutation[a_, b_]:=a.b-b.a;
InnerProduct[a_, b_]:=Tr[a.b];

(* Standard Basis *)
T1[i_, j_, n_]:=Block[{m=ConstantArray[0,{n,n}], a, b}, 
    {a, b} = If[i<j, {i,j}, {j,i}];
    m[[a,b]] += 1/2; 
    m[[b,a]] += 1/2;
    m
];
T2[i_, j_, n_]:=Block[{m=ConstantArray[0,{n,n}], a, b}, 
    {a, b} = If[i<j, {i,j}, {j,i}];
    m[[a,b]] -= I/2; 
    m[[b,a]] += I/2;
    m
];
T3[i_, n_]:=Block[{d=ConstantArray[0,n]}, 
    d[[1;;i-1]] = 1/Sqrt[2i*(i-1)];
    d[[i]] = -Sqrt[(i-1)/(2i)];
    DiagonalMatrix[d]
];
AllT[n_]:=Block[{list={}, i, j},
    For[j=2, j<=n, j++,
        For[i=1, i<j, i++,
            AppendTo[list, T1[i,j,n]];
            AppendTo[list, T2[i,j,n]];
        ];
        AppendTo[list, T3[j,n]];
    ];
    list
];

(* Cartan-Weyl Basis *)
CWE[i_, j_ , n_]:=Block[{m=ConstantArray[0,{n,n}], a, b}, 
    {a, b} = If[i<j, {i,j}, {j,i}];
    m[[a,b]] = 1;
    m
];
CWF[i_, j_ , n_]:=Block[{m=ConstantArray[0,{n,n}], a, b}, 
    {a, b} = If[i<j, {i,j}, {j,i}];
    m[[b,a]] = 1;
    m
];
CWH[i_, n_]:=Sqrt[2] * T3[n-i+1, n];
AllCW[n_]:=Block[{hl, epl={}, eml={}, i, j},
    hl = Table[CWH[i, n], {i,n-1}];
    For[j=2, j<=n, j++,
        For[i=1, i<j, i++,
            AppendTo[epl, CWE[i,j,n]];
            AppendTo[eml, CWF[i,j,n]];
        ];
    ];
    {hl,epl,eml}
];

(* Chevalley Basis *)
ChE[i_, n_]:=CWE[i, i+1, n];
ChF[i_, n_]:=CWF[i, i+1, n];
ChH[i_, n_]:=Block[{d=ConstantArray[0,n]}, 
    d[[i]] = 1;
    d[[i+1]] = -1;
    DiagonalMatrix[d]
];
AllCh[n_]:={
    Table[ChE[i,n], {i,n-1}],
    Table[ChF[i,n], {i,n-1}],
    Table[ChH[i,n], {i,n-1}]
};

(* Root System *)
RootSystem[n_]:=Block[{m=ConstantArray[0, {n-1,n-1}]},
    m[[1, n-1]] = Sqrt[2];
    For[i=2, i< n, i++,
        m[[i,   n-i]] =  Sqrt[(i+1)/i];
        m[[i, n-i+1]] = -Sqrt[(i-1)/i];
    ];
    m
];

(* Cartan Matrix *)
CartanMatrix[n_]:=Block[{d,offd},
    d = DiagonalMatrix@ConstantArray[2, n];
    offd = DiagonalMatrix[ConstantArray[-1, n-1],1];
    d + offd + Transpose[offd]
];

(* Outer Function *)
SU[n_Integer]:=AllT[n];
SU[n_Integer,"CartanWeyl"]:=AllCW[n];
SU[n_Integer,"Simple"]:=AllCh[n];

SU[n_Integer, "Roots"]:=RootSystem[n];
SU[n_Integer, "SimpleRoots"]:=RootSystem[n];
SU[n_Integer, "CartanMatrix"]:=CartanMatrix[n];

SU[n_Integer, "E", i_]:=ChE[i,n];
SU[n_Integer, "F", i_]:=ChF[i,n];
SU[n_Integer, "H", i_]:=ChH[i,n];
SU[n_Integer, "T1", {i_, j_}]:=T1[i, j, n];
SU[n_Integer, "T2", {i_, j_}]:=T2[i, j, n];
SU[n_Integer, "T3", i_]:=T3[i, n];
SU[n_Integer, "CWE", {i_, j_}]:=CWE[i, j, n];
SU[n_Integer, "CWF", {i_, j_}]:=CWF[i, j, n];
SU[n_Integer, "CWH", i_]:=CWH[i, n];


End[];
EndPackage[];