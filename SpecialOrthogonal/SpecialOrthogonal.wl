(* ::Package:: *)
BeginPackage["SpecialOrthogonal`"];

SO::usage="SO[n, string, input] returns properties of SO(N)"

Begin["Private`"];
(* Helper FUnctiond *)
Commutation[a_, b_]:=a.b-b.a;
InnerProduct[a_, b_]:=Tr[a.b];

(* Standard Basis *)
T[i_, j_, n_]:=Block[{m=ConstantArray[0,{n,n}], a, b}, 
    {a, b} = If[i<j, {i,j}, {j,i}];
    m[[a,b]] -= I; 
    m[[b,a]] += I;
    m
];
AllT[n_]:=Block[{list={}, i, j},
    For[j=2, j<=n, j++,
        For[i=1, i<j, i++,
            AppendTo[list, T[i,j,n]];
        ];
    ];
    list
];

(* Cartan-Weyl Basis *)
CWE1[i_, j_, n_]:=Block[{m=ConstantArray[0,{n,n}], a, b}, 
    {a, b} = If[i<j, {i,j}, {j,i}];
    m[[2a  ,2b-1]] =  1/2;
    m[[2a-1,2b-1]] = -I/2;
    m[[2a  ,2b  ]] = -I/2;
    m[[2a-1,2b  ]] = -1/2;
    -I*m + I*Transpose[m]
];
CWE2[i_, j_, n_]:=Block[{m=ConstantArray[0,{n,n}], a, b}, 
    {a, b} = If[i<j, {i,j}, {j,i}];
    m[[2a  ,2b-1]] =  1/2;
    m[[2a-1,2b-1]] =  I/2;
    m[[2a  ,2b  ]] =  I/2;
    m[[2a-1,2b  ]] = -1/2;
    -I*m + I*Transpose[m]
];
CWE3[i_, j_, n_]:=Block[{m=ConstantArray[0,{n,n}], a, b}, 
    {a, b} = If[i<j, {i,j}, {j,i}];
    m[[2a  ,2b-1]] =  1/2;
    m[[2a-1,2b-1]] = -I/2;
    m[[2a  ,2b  ]] =  I/2;
    m[[2a-1,2b  ]] =  1/2;
    -I*m + I*Transpose[m]
];
CWE4[i_, j_, n_]:=Block[{m=ConstantArray[0,{n,n}], a, b}, 
    {a, b} = If[i<j, {i,j}, {j,i}];
    m[[2a  ,2b-1]] =  1/2;
    m[[2a-1,2b-1]] =  I/2;
    m[[2a  ,2b  ]] = -I/2;
    m[[2a-1,2b  ]] =  1/2;
    -I*m + I*Transpose[m]
];
CWE5[i_, n_]:=Block[{m=ConstantArray[0,{n,n}]}, 
    m[[2i  ,n]] =  1/Sqrt[2];
    m[[2i-1,n]] = -I/Sqrt[2];
    -I*m + I*Transpose[m]
];
CWE6[i_, n_]:=Block[{m=ConstantArray[0,{n,n}]}, 
    m[[2i  ,n]] = 1/Sqrt[2];
    m[[2i-1,n]] = I/Sqrt[2];
    -I*m + I*Transpose[m]
];
CWH[i_, n_]:=T[2i-1,2i,n];

(* Chevalley Basis *)
ChHo[i_, n_]:=Block[{m=ConstantArray[0,{n,n}], l=(n-1)/2},
    If[i==l,
        m[[2l-1, 2l]] -= 2I; 
        m[[2l, 2l-1]] += 2I;
        ,
        m[[2i-1, 2i]] -= I; 
        m[[2i, 2i-1]] += I;
        m[[2i+1, 2i+2]] += I; 
        m[[2i+2, 2i+1]] -= I;
    ];
    m
];
ChHe[i_, n_]:=Block[{m=ConstantArray[0,{n,n}], l=n/2},
    If[i==l,
        m[[2l-1, 2l]] -= I; 
        m[[2l, 2l-1]] += I;
        m[[2l-3, 2l-2]] -= I; 
        m[[2l-2, 2l-3]] += I;
        ,
        m[[2i-1, 2i]] -= I; 
        m[[2i, 2i-1]] += I;
        m[[2i+1, 2i+2]] += I; 
        m[[2i+2, 2i+1]] -= I;
    ];
    m
];
ChEo[i_, n_]:=If[i==(n-1)/2, Sqrt[2]*CWE5[i,n], CWE1[i,i+1,n]];
ChEe[i_, n_]:=If[i==n/2, CWE3[n/2-1,n/2,n], CWE1[i,i+1,n]];
ChFo[i_, n_]:=If[i==(n-1)/2, Sqrt[2]*CWE6[i,n], CWE2[i,i+1,n]];
ChFe[i_, n_]:=If[i==n/2, CWE4[n/2-1,n/2,n], CWE2[i,i+1,n]];
AllCh[n_]:=Block[{h,e,f},
    If[EvenQ[n], 
        l = n/2;
        h = Table[ChHe[i, n], {i,l}];
        e = Table[ChEe[i, n], {i,l}];
        f = Table[ChFe[i, n], {i,l}];
        , 
        l = (n-1)/2;
        h = Table[ChHo[i, n], {i,l}];
        e = Table[ChEo[i, n], {i,l}];
        f = Table[ChFo[i, n], {i,l}];
    ];
    {h,e,f}
];

(* Root System *)
RootSystem[i_, n_]:=Block[{l,rt=ConstantArray[0, {n,n}]},
    If[EvenQ[n],
        If[i==n/2, 
            rt[[i,i  ]] = 1, 
            rt[[i,i  ]] = 1; 
            rt[[i,i+1]] = -1
        ],
        If[i==(n-1)/2,
            rt[[i,i-1]] = 1; 
            rt[[i,i  ]] = -1,
            rt[[i,i  ]] = 1; 
            rt[[i,i+1]] = -1
        ];
    ];
    rt
];

(* Spherical Harmonic Basis *)
Basiso[n_]:=Block[{l=(n-1)/2, b},
    b = ConstantArray[0, {n,n}];
    For[i=1, i<=l, i++,
        b[[i,2i-1]] = (-1)^(l-i+1)/Sqrt[2];
        b[[i,2i  ]] = (-1)^(l-i+1)/Sqrt[2] * I;
    ];
    b[[l+1,n]] = 1;
    For[i=l+2, i<=n, i++,
        b[[i, 4l-2i+3]] = 1/Sqrt[2];
        b[[i, 4l-2i+4]] = -I/Sqrt[2];
    ];
    b
];
Basise[n_]:=Block[{l=n/2, b},
    b = ConstantArray[0, {n,n}];
    For[i=1, i<=l, i++,
        b[[i,2i-1]] = (-1)^(l-i)/Sqrt[2];
        b[[i,2i  ]] = (-1)^(l-i)/Sqrt[2] * I;
    ];
    For[i=l+1, i<=n, i++,
        b[[i,4l-2i+1]] = 1/Sqrt[2];
        b[[i,4l-2i+2]] = -I/Sqrt[2];
    ];
    b
];

(* Outer Function *)
SO[n_Integer]:=AllT[n];
SO[n_Integer, "Simple"]:=AllCh[n];
SO[n_Integer, "T", l_]:=T[l[[1]], l[[2]], n];
SO[n_Integer, "Basis"]:=If[EvenQ[n], Basise[n], Basiso[n]];

End[];
EndPackage[];
