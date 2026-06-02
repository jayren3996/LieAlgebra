#!/usr/bin/env wolframscript
(* ::Package:: *)

(* MakeFigures.wls
   Regenerates every figure used in README.md from the ClassicalLieAlgebra
   package, so the figures stay accurate to the code.

   Run with a kernel that can reach a front end (needed to rasterize
   typeset matrices / Young tableaux), e.g. on macOS:

     /Applications/Wolfram.app/Contents/MacOS/wolframscript -file pics/MakeFigures.wls

   Output: PNG files written next to this script (the pics/ folder).

   The SU(3) (1,1)-representation expressions below are not hand-typed
   constants: each one is checked against the package in MakeFigures
   verification (see the Print[...] accuracy report at the end). *)

scriptDir = If[$InputFileName === "", Directory[], DirectoryName[$InputFileName]];
root = ParentDirectory[scriptDir];
Get[FileNameJoin[{root, "ClassicalLieAlgebra.wl"}]];
out[name_] := FileNameJoin[{scriptDir, name}];

(* ----- generator action on tensor (Psi) states: F_total = sum_sites G ----- *)
actPsi[g_, Psi[a___]] := Module[{l = {a}, n = Length[g], s, j},
   Sum[g[[j, l[[s]]]] (Psi @@ ReplacePart[l, s -> j]), {s, Length[l]}, {j, n}]];
act[g_, expr_] := Expand[expr /. p_Psi :> actPsi[g, p]];

(* ----- small helpers for the schematic diagrams ----- *)
tabSize[rows_, c_] := {Max[Length /@ rows] c, Length[rows] c};

(* draw a Young tableau as white-filled boxes, top-left corner at {x0,y0} *)
drawTab[rows_, {x0_, y0_}, c_, fs_] := Flatten @ Table[
   With[{r = rr, k = kk, val = rows[[rr, kk]]},
    {EdgeForm[{Black, Thickness[0.004]}], FaceForm[White],
     Rectangle[{x0 + (k - 1) c, y0 - r c}, {x0 + k c, y0 - (r - 1) c}],
     Text[Style[val, fs, FontFamily -> "Times"],
      {x0 + (k - 0.5) c, y0 - (r - 0.5) c}]}],
   {rr, Length[rows]}, {kk, Length[rows[[rr]]]}];

(* same, but centered at {cx,cy} *)
drawTabC[rows_, {cx_, cy_}, c_, fs_] := Module[{w, h},
   {w, h} = tabSize[rows, c];
   drawTab[rows, {cx - w/2, cy + h/2}, c, fs]];

(* two parallel lines between a and b (the "F2" double-line convention) *)
doubleLine[a_, b_, d_] := Module[{dir, perp},
   dir = Normalize[b - a]; perp = d {-dir[[2]], dir[[1]]};
   {Line[{a + perp/2, b + perp/2}], Line[{a - perp/2, b - perp/2}]}];

green = Darker[Green, 0.2];

(* =====================================================================
   Schematic diagrams (pure Graphics)
   ===================================================================== *)

(* ---- su3-levels.png : the SU(3) three-level system (was P1) ---- *)
yOf = <|1 -> 2, 2 -> 1, 3 -> 0|>;
su3levels = Graphics[{
    Thickness[0.006], Black,
    Line[{{0, 2}, {1.8, 2}}], Line[{{0, 1}, {1.8, 1}}], Line[{{0, 0}, {1.8, 0}}],
    Text[Style[Row[{"\[LeftBracketingBar]", #, "\[RightAngleBracket]"}], 22], {-0.25, yOf[#]}, {1, 0}] & /@ {1, 2, 3},
    (* E_i : raising, red, up *)
    Red, Arrowheads[0.06], Thickness[0.007],
    Arrow[{{2.15, 1.12}, {2.15, 1.88}}], Arrow[{{2.15, 0.12}, {2.15, 0.88}}],
    Text[Style[Subscript["E", 1], 20], {1.92, 1.5}],
    Text[Style[Subscript["E", 2], 20], {1.92, 0.5}],
    (* F_i : lowering, green, down *)
    green, Arrowheads[0.06],
    Arrow[{{2.7, 1.88}, {2.7, 1.12}}], Arrow[{{2.7, 0.88}, {2.7, 0.12}}],
    Text[Style[Subscript["F", 1], 20], {2.95, 1.5}],
    Text[Style[Subscript["F", 2], 20], {2.95, 0.5}],
    (* weights *)
    Black,
    MapThread[
     Text[Style[Row[{Subscript["H", 1], " = ", #2, ",  ", Subscript["H", 2], " = ", #3}], 18], {3.9, yOf[#1]}, {-1, 0}] &,
     {{1, 2, 3}, {1, -1, 0}, {0, 1, -1}}]
    },
   ImageSize -> 560, PlotRange -> {{-0.7, 7.2}, {-0.5, 2.5}}, Background -> White];

(* ---- rep11-first-level.png : weight-lowering tree (was T1) ---- *)
cTop = {0, 2}; cL = {-1.7, 0}; cR = {1.7, 0}; cc = 0.5; fs = 18;
rep11tree = Graphics[{
    Thickness[0.005], Black,
    Line[{cTop + {-0.45, -0.35}, cL + {0.3, 0.45}}],          (* F1 single *)
    doubleLine[cTop + {0.45, -0.35}, cR + {-0.3, 0.45}, 0.07], (* F2 double *)
    drawTabC[{{1, 1}, {2}}, cTop, cc, fs],
    drawTabC[{{1, 2}, {2}}, cL, cc, fs],
    drawTabC[{{1, 1}, {3}}, cR, cc, fs],
    Text[Style["(1,1)", fs], cTop + {0.95, 0.2}, {-1, 0}],
    Text[Style["(-1,2)", fs], cL + {-0.7, 0}, {1, 0}],
    Text[Style["(2,-1)", fs], cR + {0.7, 0}, {-1, 0}],
    Text[Style[Subscript["F", 1], fs], {-1.15, 1.2}],
    Text[Style[Subscript["F", 2], fs], {1.15, 1.2}]
    },
   ImageSize -> 460, PlotRange -> {{-3, 3}, {-0.7, 2.7}}, Background -> White];

(* ---- rep11-transition.png : the |a> -> {|c>,|d>} transition (was T4) ---- *)
aTop = {0, 2}; cNode = {-1.5, 0}; dNode = {1.9, 0};
rep11trans = Graphics[{
    Thickness[0.005], Black,
    doubleLine[aTop + {-0.35, -0.35}, cNode + {0.25, 0.45}, 0.07],
    doubleLine[aTop + {0.4, -0.35}, {1.0, 0.5}, 0.07],
    drawTabC[{{1, 2}, {2}}, aTop, cc, fs],
    drawTabC[{{1, 3}, {2}}, cNode, cc, fs],
    Text[Style[Subscript["t", 1], fs], {-1.15, 1.2}],
    Text[Style[Subscript["t", 2], fs], {1.25, 1.25}],
    (* |d> = 2[[1,2],[3]] - [[1,3],[2]] in a frame *)
    EdgeForm[{Black, Thickness[0.004]}], FaceForm[None],
    Rectangle[{1.05, -0.7}, {3.8, 0.7}],
    Text[Style["2", fs], {1.35, 0}],
    drawTabC[{{1, 2}, {3}}, {1.95, 0.02}, cc, fs],
    Text[Style["-", fs], {2.55, 0.02}],
    drawTabC[{{1, 3}, {2}}, {3.15, 0.02}, cc, fs]
    },
   ImageSize -> 640, PlotRange -> {{-2.4, 4.0}, {-0.9, 2.7}}, Background -> White];

(* =====================================================================
   Typeset figures (matrices / tableaux) -- need a front end to rasterize
   ===================================================================== *)
tb[rows_] := TableauForm[TensorTableau[rows]];
matRow[lbl_, mats_] := Row[Join[{Style[lbl <> " = ", Bold]}, Riffle[MatrixForm /@ mats, "  "]]];

cw = CartanWeyl[SU[3]];
ch = Chevalley[SU[3]];

(* (1,1) rep tableau-basis objects *)
T = TensorTableau;
o4 = TableauForm /@ TableauOrthogonalization[
    2 T[{{1, 2}, {3}}] - T[{{1, 3}, {2}}], T[{{1, 3}, {2}}] + T[{{1, 2}, {3}}]];
naState = TableauNormalization[T[{{1, 2}, {2}}]];
(* distributed form so TableauForm / ToTensor (which do not auto-expand a
   scalar over a sum) match on each TensorTableau term *)
bState = T[{{1, 2}, {3}}]/Sqrt[6] + T[{{1, 3}, {2}}]/Sqrt[6];
eq = Style["  =  ", 18];

typeset = <|
   "su3-generators.png"      -> Grid[Partition[MatrixForm /@ Generators[SU[3]], 4], Spacings -> {2, 1.5}],
   "su3-cartan-weyl.png"     -> Column[{matRow["H", cw[[1]]], matRow["E", cw[[2]]], matRow["F", cw[[3]]]}, Spacings -> 1.4],
   "su3-chevalley.png"       -> Column[{matRow["H", ch[[1]]], matRow["E", ch[[2]]], matRow["F", ch[[3]]]}, Spacings -> 1.4],
   "rep11-highest-weight.png"-> tb[{{1, 1}, {2}}],
   "rep11-state-a.png"       -> tb[{{1, 2}, {2}}],
   "rep11-f2-superposition.png" -> TableauForm[T[{{1, 2}, {3}}] + T[{{1, 3}, {2}}]],
   "rep11-f1-result.png"     -> TableauForm[T[{{1, 2}, {3}}] + T[{{2, 1}, {3}}]],
   "rep11-symmetrize.png"    -> Row[{TableauForm[T[{{2, 1}, {3}}]], eq,
                                     TableauForm[T[{{1, 2}, {3}}] + T[{{2, 3}, {1}}]], eq,
                                     TableauForm[T[{{1, 2}, {3}}] - T[{{1, 3}, {2}}]]}],
   "rep11-regular.png"       -> Row[{TableauForm[T[{{1, 2}, {3}}] + T[{{2, 1}, {3}}]], eq,
                                     TableauForm[2 T[{{1, 2}, {3}}] - T[{{1, 3}, {2}}]]}],
   "rep11-orthogonalized.png"-> Row[{"{  ", o4[[1]], "  ,  ", o4[[2]], "  }"}],
   "rep11-state-a-normalized.png" -> TableauForm[naState],
   "rep11-f2-action.png"     -> TableauForm[bState]
   |>;

graphicsFigs = <|
   "su3-levels.png"        -> su3levels,
   "rep11-first-level.png" -> rep11tree,
   "rep11-transition.png"  -> rep11trans
   |>;

UsingFrontEnd[
  KeyValueMap[Export[out[#1], #2, ImageResolution -> 150] &, typeset];
  KeyValueMap[Export[out[#1], #2, ImageResolution -> 144] &, graphicsFigs];
];

(* =====================================================================
   Accuracy report: confirm the (1,1) expressions match the package
   ===================================================================== *)
hw = ToTensor[T[{{1, 1}, {2}}]];
chk[label_, b_] := Print[If[b, "  ok  ", " FAIL "], label];
Print["--- accuracy checks ---"];
chk["F1.hw == {{1,2},{2}}", Expand[act[ch[[3, 1]], hw] - ToTensor[T[{{1, 2}, {2}}]]] === 0];
chk["F2.hw == {{1,1},{3}}", Expand[act[ch[[3, 2]], hw] - ToTensor[T[{{1, 1}, {3}}]]] === 0];
chk["F2.{{1,2},{2}} == {{1,2},{3}}+{{1,3},{2}}", Expand[act[ch[[3, 2]], ToTensor[T[{{1, 2}, {2}}]]] - ToTensor[T[{{1, 2}, {3}}] + T[{{1, 3}, {2}}]]] === 0];
chk["{{2,1},{3}} == {{1,2},{3}}-{{1,3},{2}}", Expand[ToTensor[T[{{2, 1}, {3}}]] - ToTensor[T[{{1, 2}, {3}}] - T[{{1, 3}, {2}}]]] === 0];
chk["F2.|a> == |b>", Expand[act[ch[[3, 2]], ToTensor[naState]] - ToTensor[bState]] === 0];
chk["<a|a> == 1", TableauDot[naState, naState] === 1];
Print["files written to ", scriptDir];
Print["exported: ", FileNameTake /@ Join[Keys[typeset], Keys[graphicsFigs]]];
