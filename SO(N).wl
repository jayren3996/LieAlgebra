(* ::Package:: *)

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
        For[i=l+1, i<=2l, i++,
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
n=6;
{Hn,En,Fn}=SO[n];
TraditionalForm/@Hn
TraditionalForm/@En
TraditionalForm/@Fn
RootSystem[Hn,En]
nb=StandardBasis[n]
chb=TraditionalForm[Conjugate[nb].#.Transpose[nb]] &;
chb/@Hn
chb/@En
chb/@Fn



