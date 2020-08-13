(* ::Package:: *)

BeginPackage["SymplecticUnitary`"]

SpT

SpCWH
SpCWE
SpCWF
SpCW

SpChH
SpChE
SpChF
SpCh

SpBasis

SpH
SpE
SpF
Spn

Begin["Private`"];
(* Standard Generators *)
SpT[{n_,1},{i_, j_}]:=Block[{m=ConstantArray[0,{n,n}], a, b}, 
    {a, b} = If[i<j, {i,j}, {j,i}];
    m[[2a-1,2b-1]] = -I/2;
    m[[2b-1,2a-1]] = +I/2;
    m[[2a  ,2b  ]] = -I/2;
    m[[2b  ,2a  ]] = -I/2;
    m
];
SpT[{n_,2},{i_, j_}]:=Block[{m=ConstantArray[0,{n,n}], a, b}, 
    {a, b} = If[i<j, {i,j}, {j,i}];
    m[[2a-1,2b  ]] = 1/2;
    m[[2b-1,2a  ]] = 1/2;
    m[[2a  ,2b-1]] = 1/2;
    m[[2b  ,2a-1]] = 1/2;
    m
];
SpT[{n_,3},{i_, j_}]:=Block[{m=ConstantArray[0,{n,n}], a, b}, 
    {a, b} = If[i<j, {i,j}, {j,i}];
    m[[2a-1,2b  ]] = -I/2;
    m[[2b-1,2a  ]] = -I/2;
    m[[2a  ,2b-1]] = +I/2;
    m[[2b  ,2a-1]] = +I/2;
    m
];
SpT[{n_,4},{i_, j_}]:=Block[{m=ConstantArray[0,{n,n}], a, b}, 
    {a, b} = If[i<j, {i,j}, {j,i}];
    m[[2a-1,2b-1]] = +1/2;
    m[[2b-1,2a-1]] = +1/2;
    m[[2a  ,2b  ]] = -1/2;
    m[[2b  ,2a  ]] = -1/2;
    m
];
SpT[{n_,5},{i_}]:=Block[{m=ConstantArray[0,{n,n}]}, 
    m[[2i-1,2i  ]] = +1/Sqrt[2];
    m[[2i  ,2i-1]] = +1/Sqrt[2];
    m
];
SpT[{n_,6},{i_}]:=Block[{m=ConstantArray[0,{n,n}]}, 
    m[[2i-1,2i  ]] = -I/Sqrt[2];
    m[[2i  ,2i-1]] = +I/Sqrt[2];
    m
];
SpT[{n_,7},{i_}]:=Block[{m=ConstantArray[0,{n,n}]}, 
    m[[2i-1,2i-1]] = +1/Sqrt[2];
    m[[2i  ,2i  ]] = -1/Sqrt[2];
    m
];
SpT[n_]:=Block[{l=n/2, Tl={}, i, j},
    For[j=2, j<=l, j++,
        For[i=1, i<j, i++,
            AppendTo[Tl, SpT[{n,1},{i,j}]];
            AppendTo[Tl, SpT[{n,2},{i,j}]];
            AppendTo[Tl, SpT[{n,3},{i,j}]];
            AppendTo[Tl, SpT[{n,4},{i,j}]];
        ];
    ];
    For[i=1,i<=l,i++,
        AppendTo[Tl, SpT[{n,5},{i}]];
        AppendTo[Tl, SpT[{n,6},{i}]];
        AppendTo[Tl, SpT[{n,7},{i}]];
    ];
    Tl
];

(* Cartan-Weyl Basis *)
SpCWH[n_,i_]:=Block[{m=ConstantArray[0,{n,n}]}, 
    m[[2i-1,2i-1]] = 1/Sqrt[2];
    m[[2i  ,2i  ]] = -1/Sqrt[2];
    m
];
SpCWE[{n_,1},{i_, j_}]:=Block[{m=ConstantArray[0,{n,n}], a, b}, 
    {a, b} = If[i<j, {i,j}, {j,i}];
    m[[2a-1,2b-1]] = +1/Sqrt[2];
    m[[2b  ,2a  ]] = -1/Sqrt[2];
    m
];
SpCWE[{n_,2},{i_, j_}]:=Block[{m=ConstantArray[0,{n,n}], a, b}, 
    {a, b} = If[i<j, {i,j}, {j,i}];
    m[[2a-1,2b]] = 1/Sqrt[2];
    m[[2b-1,2a]] = 1/Sqrt[2];
    m
];
SpCWE[{n_,3},{i_}]:=Block[{m=ConstantArray[0,{n,n}]}, 
    m[[2i-1,2i]] = 1;
    m
];
SpCWF[{n_,i_}, index]:= Conjugate @ Transpose @ SpCWE[{n,i}, index];
SpCW[n_]:=Block[{l=n/2, Hl, El={}, Fl, i, j},
    Hl = Table[SpCWH[n, i], {i,l}];
    For[j=2, j<=l, j++,
        For[i=1, i<j, i++,
            AppendTo[El, SpCWE[{n,1},{i,j}]];
            AppendTo[El, SpCWE[{n,2},{i,j}]];
        ];
    ];
    For[i=1,i<=l,i++,
        AppendTo[El, SpCWE[{n,3},{i}]];
    ];
    Fl = Conjugate @* Transpose /@ El;
    {Hl, El, Fl}
];

(* Chevalley Basis *)
SpChH[n_,i_]:=Block[{m=ConstantArray[0,{n,n}]}, 
    If[i==n/2,
        m[[2i-1,2i-1]] = +1;
        m[[2i  ,2i  ]] = -1;
        ,
        m[[2i-1,2i-1]] = +1;
        m[[2i  ,2i  ]] = -1;
        m[[2i+1,2i+1]] = -1;
        m[[2i+2,2i+2]] = +1;
    ];
    m
];
SpChE[n_,i_]:=Block[{m=ConstantArray[0,{n,n}]}, 
    If[i==n/2,
        m[[2i-1,2i]] = 1;
        ,
        m[[2i-1,2i+1]] = +1;
        m[[2i+2,2i  ]] = -1;
    ];
    m
];
SpChF[n_,i_]:= Conjugate @ Transpose @ SpChE[n, i];
SpCh[n_]:=Block[{l=n/2, Hl, El={}, Fl, i},
    Hl = Table[SpChH[n,i], {i,l}];
    El = Table[SpChE[n,i], {i,l}];
    Fl = Conjugate @* Transpose /@ El;
    {Hl, El, Fl}
];

(* Permuted Basis *)
SpBasis[n_]:= Block[{m=ConstantArray[0, {n,n}], l=n/2, i},
    For[i=1,i<=l,i++,
        m[[2i-1, i]] = 1;
        m[[n+2-2i, l+i]] = 1;
    ];
    m
];

(* Chevalley in new Basis *)
SpH[n_,i_]:=Block[{m=ConstantArray[0,{n,n}]}, 
    If[i==n/2,
        m[[i,i]] = 1;
        m[[i+1,i+1]] = -1;
        ,
        m[[i,i]] = 1;
        m[[i+1,i+1]] = -1;
        m[[n-i,n-i]] = 1;
        m[[n-i+1,n-i+1]] = -1;
    ];
    m
];
SpE[n_,i_]:=Block[{m=ConstantArray[0,{n,n}]}, 
    If[i==n/2,
        m[[i,i+1]] = 1;
        ,
        m[[i,i+1]] = 1;
        m[[n-i,n-i+1]] = -1;
    ];
    m
];
SpF[n_,i_]:= Conjugate @ Transpose @ SpE[n, i];
Spn[n_]:=Block[{l=n/2, Hl, El={}, Fl},
    Hl = Table[SpH[n,i], {i,l}];
    El = Table[SpE[n,i], {i,l}];
    Fl = Conjugate @* Transpose /@ El;
    {Hl, El, Fl}
];

End[];
EndPackage[];
