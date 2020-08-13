(* ::Package:: *)
BeginPackage["SpecialOrthogonal`"];

SOT

SOCWH
SOCWE
SOCWF
SOCW

SOChH
SOChE
SOChF
SOCh

SOBasis

SOH
SOE
SOF
SOn

Begin["Private`"];

(* Standard Basis *)
SOT[n_, {i_, j_}]:=Block[{m=ConstantArray[0,{n,n}], a, b}, 
    {a, b} = If[i<j, {i,j}, {j,i}];
    m[[a,b]] -= I; 
    m[[b,a]] += I;
    m
];
SOT[n_]:=Block[{list={}, i, j},
    For[j=2, j<=n, j++,
        For[i=1, i<j, i++,
            AppendTo[list, SOT[n,{i,j}]];
        ];
    ];
    list
];

(* Cartan-Weyl Basis *)
SOCWH[n_,i_]:=Block[{m=ConstantArray[0,{n,n}]},
    m[[2i-1,2i  ]] = -I;
    m[[2i  ,2i-1]] = +I;
    m
];
SOCWE[{n_,1},{i_, j_}]:=Block[{m=ConstantArray[0,{n,n}], a, b}, 
    {a, b} = If[i<j, {i,j}, {j,i}];
    m[[2a  ,2b-1]] = +1/2;
    m[[2a-1,2b-1]] = -I/2;
    m[[2a  ,2b  ]] = -I/2;
    m[[2a-1,2b  ]] = -1/2;
    -I*m + I*Transpose[m]
];
SOCWE[{n_,2},{i_, j_}]:=Block[{m=ConstantArray[0,{n,n}], a, b}, 
    {a, b} = If[i<j, {i,j}, {j,i}];
    m[[2a  ,2b-1]] = +1/2;
    m[[2a-1,2b-1]] = -I/2;
    m[[2a  ,2b  ]] = +I/2;
    m[[2a-1,2b  ]] = +1/2;
    -I*m + I*Transpose[m]
];
SOCWE[{n_,3},{i_}]:=Block[{m=ConstantArray[0,{n,n}]}, 
    m[[2i  ,n]] = +1/Sqrt[2];
    m[[2i-1,n]] = -I/Sqrt[2];
    -I*m + I*Transpose[m]
];
SOCWF[{n_,1},{i_, j_}]:= Conjugate @ Transpose @ SOCWE[{n,1},{i,j}];
SOCWF[{n_,2},{i_, j_}]:= Conjugate @ Transpose @ SOCWE[{n,2},{i,j}];
SOCWF[{n_,3},{i_}    ]:= Conjugate @ Transpose @ SOCWE[{n,3},{i}  ];
SOCW[n_]:=Block[{l,p,i,j,h,e={},f},
    {l,p} = QuotientRemainder[n,2];
    h = Table[SOCWH[n,i], {i,l}];
    For[j=2, j<=l, j++,
        For[i=1, i<j, i++,
            AppendTo[e, SOCWE[{n,1},{i,j}]];
            AppendTo[e, SOCWE[{n,2},{i,j}]];
        ];
    ];
    If[p==1, 
        For[i=1,i<=l,i++,
            AppendTo[e, SOCWE[{n,3}, {i}]];
        ];
    ];
    f = Conjugate @* Transpose /@ e;
    {h,e,f}
];

(* Chevalley Basis *)
SOChH[2,1]:={{0,-I},{I,0}};
SOChH[n_?OddQ,i_]:=Block[{m=ConstantArray[0,{n,n}], l=(n-1)/2},
    If[i==l,
        m[[2l-1, 2l]] = -2I; 
        m[[2l, 2l-1]] = +2I;
        ,
        m[[2i-1, 2i  ]] = -I; 
        m[[2i  , 2i-1]] = +I;
        m[[2i+1, 2i+2]] = +I; 
        m[[2i+2, 2i+1]] = -I;
    ];
    m
];
SOChH[n_?EvenQ,i_]:=Block[{m=ConstantArray[0,{n,n}], l=n/2},
    If[i==l,
        m[[2l-1, 2l  ]] = -I; 
        m[[2l  , 2l-1]] = +I;
        m[[2l-3, 2l-2]] = -I; 
        m[[2l-2, 2l-3]] = +I;
        ,
        m[[2i-1, 2i  ]] = -I; 
        m[[2i  , 2i-1]] = +I;
        m[[2i+1, 2i+2]] = +I; 
        m[[2i+2, 2i+1]] = -I;
    ];
    m
];
SOChE[2,1]:= Nothing;
SOChE[n_?OddQ,i_]:=Block[{m=ConstantArray[0,{n,n}], l=(n-1)/2},
    If[i==l,
        m[[2i  ,n]] = +1;
        m[[2i-1,n]] = -I;
        ,
        m[[2i  ,2i+1]] = +1/2;
        m[[2i-1,2i+1]] = -I/2;
        m[[2i  ,2i+2]] = -I/2;
        m[[2i-1,2i+2]] = -1/2;
    ];
    -I*m + I*Transpose[m]
];
SOChE[n_?EvenQ,i_]:=Block[{m=ConstantArray[0,{n,n}], l=n/2},
    If[i==l,
        m[[2i-2,2i-1]] = +1/2;
        m[[2i-3,2i-1]] = -I/2;
        m[[2i-2,2i  ]] = +I/2;
        m[[2i-3,2i  ]] = +1/2;
        ,
        m[[2i  ,2i+1]] =  1/2;
        m[[2i-1,2i+1]] = -I/2;
        m[[2i  ,2i+2]] = -I/2;
        m[[2i-1,2i+2]] = -1/2;
    ];
    -I*m + I*Transpose[m]
];
SOChF[n_,i_]:= Conjugate @ Transpose @ SOChE[n,i];
SOCh[n_]:=Block[{l,h,e,f},
    l = Quotient[n,2];
    h = Table[SOChH[n,i], {i,l}];
    e = Table[SOChE[n,i], {i,l}];
    f = Conjugate @* Transpose /@ e;
    {h,e,f}
];

(* Spherical Harmonic Basis *)
SOBasis[n_?OddQ]:=Block[{l=(n-1)/2, b},
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
SOBasis[n_?OddQ]:=Block[{l=n/2, b},
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

(* Chevalley in new Basis *)
SOH[n_?OddQ, i_]:=Block[{m=ConstantArray[0,{n,n}], l=(n-1)/2},
    If[i==l,
        m[[i  , i  ]] = +2; 
        m[[i+2, i+2]] = -2;
        ,
        m[[i  , i  ]] = +1; 
        m[[i+1, i+1]] = -1;
        m[[n-i, n-i]] = +1; 
        m[[n-i+1, n-i+1]] = -1;
    ];
    m
];
SOE[n_?OddQ, i_]:=Block[{m=ConstantArray[0,{n,n}], l=(n-1)/2},
    If[i==l,
        m[[i  , i+1]] = +Sqrt[2]; 
        m[[i+1, i+2]] = +Sqrt[2];
        ,
        m[[i  ,   i+1]] = +1; 
        m[[n-i, n-i+1]] = +1;
    ];
    m
];
SOH[2,1]:={{1,0},{0,-1}};
SOH[n_?EvenQ, i_]:=Block[{m=ConstantArray[0,{n,n}], l=n/2},
    If[i==l,
        m[[i-1, i-1]] = +1; 
        m[[i  , i  ]] = +1;
        m[[i+1, i+1]] = -1; 
        m[[i+2, i+2]] = -1;
        ,
        m[[i  , i  ]] = +1; 
        m[[i+1, i+1]] = -1;
        m[[n-i  , n-i  ]] = +1; 
        m[[n-i+1, n-i+1]] = -1;
    ];
    m
];
SOE[2, 1]:= Nothing;
SOE[n_?EvenQ, i_]:=Block[{m=ConstantArray[0,{n,n}], l=n/2},
    If[i==l,
        m[[i-1, i+1]] = +1; 
        m[[i  , i+2]] = +1;
        ,
        m[[i  ,   i+1]] = +1; 
        m[[n-i, n-i+1]] = +1;
    ];
    m
];
SOF[n_, i_]:= Transpose @ SOE[n, i];
SOn[n_]:=Block[{l,h,e,f},
    l = Quotient[n,2];
    h = Table[SOH[n,i], {i,l}];
    e = Table[SOE[n,i], {i,l}];
    f = Transpose /@ e;
    {h,e,f}
];


End[];
EndPackage[];
