(* ::Package:: *)

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
ChH[i_, n_]:=Block[{m=ConstantArray[0,{n,n}]},
    If[EvenQ[n],
        l = n/2;
        If[i==l,
            m[[2l-1, 2l]] -= 2I; 
            m[[2l, 2l-1]] += 2I;
            m[[2l-3, 2l-2]] -= I; 
            m[[2l-2, 2l-3]] += I;
            ,
            m[[2i-1, 2i]] -= I; 
            m[[2i, 2i-1]] += I;
            m[[2i+1, 2i+2]] += I; 
            m[[2i+2, 2i+1]] -= I;
        ];
        ,
        l = (n-1)/2;
        If[i==l,
            m[[2l-1, 2l]] -= 2I; 
            m[[2l, 2l-1]] += 2I;
            ,
            m[[2i-1, 2i]] -= I; 
            m[[2i, 2i-1]] += I;
            m[[2i+1, 2i+2]] += I; 
            m[[2i+2, 2i+1]] -= I;
        ];
    ];
    m
];
ChE[i_, n_]:=Block[{m},
    l = If[EvenQ[n], n/2, (n-1)/2];
    m = If[i==n,
        If[EvenQ[n], Sqrt[2]*CWE5[l,n], ]
        ,
        CWE1[i,i+1,n]
    ];
];
ChF[i_, n_]:=Block[{m},
    l = If[EvenQ[n], n/2, (n-1)/2];
    m = If[i==n,
        l = n/2;
        ,
        CWE2[i,i+1,n]
    ];
];
SO[n_]:=Module[{l, Hi, Ei, Fi, i},
    l = If[EvenQ[n], n/2, (n-1)/2];
    Hi = Table[ConstantArray[0, {n,n}], l];
    Ei = Table[ConstantArray[0, {n,n}], l];
    Fi = Table[ConstantArray[0, {n,n}], l];
    For[i=1, i<l, i++,
        Hi[[i]][[2i-1,2i  ]] = 1;
        Ei[[i]][[2i  ,2i+1]] = 1/2;
        Ei[[i]][[2i-1,2i+1]] = -I/2;
        Ei[[i]][[2i  ,2i+2]] = -I/2;
        Ei[[i]][[2i-1,2i+2]] = -1/2;
        Fi[[i]][[2i  ,2i+1]] = 1/2;
        Fi[[i]][[2i-1,2i+1]] = I/2;
        Fi[[i]][[2i  ,2i+2]] = I/2;
        Fi[[i]][[2i-1,2i+2]] = -1/2;
    ];
    If[EvenQ[n],
        Hi[[l]][[2l-1,2l  ]] = 1;
        Ei[[l]][[2l-2,2l-1]] = 1/2;
        Ei[[l]][[2l-3,2l-1]] = -I/2;
        Ei[[l]][[2l-2,2l  ]] = I/2;
        Ei[[l]][[2l-3,2l  ]] = 1/2;
        Fi[[l]][[2l-2,2l-1]] = 1/2;
        Fi[[l]][[2l-3,2l-1]] = I/2;
        Fi[[l]][[2l-2,2l  ]] = -I/2;
        Fi[[l]][[2l-3,2l  ]] = 1/2;
        ,
        Hi[[l]][[2l-1,2l  ]] = 1;
        Ei[[l]][[2l  ,2l+1]] = 1/Sqrt[2];
        Ei[[l]][[2l-1,2l+1]] = -I/Sqrt[2];
        Fi[[l]][[2l  ,2l+1]] = 1/Sqrt[2];
        Fi[[l]][[2l-1,2l+1]] = I/Sqrt[2];
    ];
    Hi = (-I*# + I*Transpose[#] &)/@Hi;
    Ei = (-I*# + I*Transpose[#] &)/@Ei;
    Fi = (-I*# + I*Transpose[#] &)/@Fi;
    {Hi, Ei, Fi}
];

Commute[a_, b_]:=a.b - b.a;

RootSystem[Hi_, Ei_]:=Module[{l, roots, root, i, j},
    l = Length[Hi];
    roots = ConstantArray[0, {l, l}];
    For[i = 1, i <= l, i++,
        root = Commute[Ei[[i]], Conjugate@Transpose[Ei[[i]]]];
        For[j = 1, j <= l, j++,
            roots[[i, j]] = Simplify@Tr[root.Hi[[j]]/2];
        ];
    ];
    roots
];

CartanMatrix[rs_]:=Module[{l},
    l = Length[rs];
    Table[Dot[rs[[i]], rs[[j]]], {i, l}, {j, l}]
];

StandardBasis[n_]:=Module[{b, l, i},
    b = ConstantArray[0, {n,n}];
    If[EvenQ[n],
        l = n/2;
        For[i=1, i<=l, i++,
            b[[i,2i-1]] = (-1)^(l-i)/Sqrt[2];
            b[[i,2i  ]] = (-1)^(l-i)/Sqrt[2] * I;
        ];
        For[i=l+1; i<=2l, i++,
            b[[i,4l-2i+1]] = 1/Sqrt[2];
            b[[i,4l-2i+2]] = -I/Sqrt[2];
        ];
        ,
        l = (n-1)/2;
        For[i=1, i<=l, i++,
            b[[i,2i-1]] = (-1)^(l-i+1)/Sqrt[2];
            b[[i,2i  ]] = (-1)^(l-i+1)/Sqrt[2] * I;
        ];
        b[[l+1,2l+1]] = 1;
        For[i=l+2, i<=2l+1, i++,
            b[[i, 4l-2i+3]] = 1/Sqrt[2];
            b[[i, 4l-2i+4]] = -I/Sqrt[2];
        ];
    ];
    b
];
n=3;
{Hn,En,Fn}=SO[3];
TraditionalForm/@Hn
TraditionalForm/@En
TraditionalForm/@Fn
RootSystem[Hn,En]
nb=StandardBasis[n]
chb=Conjugate@Transpose[nb].#.nb &;
chb/@Hn
chb/@En
chb/@Fn



