#!/usr/bin/env wolframscript
(* generate-reference.wls -- generate the theme-grouped Markdown API reference
   for the ClassicalLieAlgebra documentation site.

   Each entry combines the symbol's own ::usage (the code-sourced contract) with
   a curated detail layer from scripts/reference-details.wl: a signature, a
   fuller description, an argument table, options, worked examples, notes, and
   "see also" links. Every example is EVALUATED here and its output embedded, so
   the reference cannot drift from the code without the build failing.

   Output:
     MkDocs (default): <repo>/site-src/reference/{index,algebras,
       root-system,representations,young-tableaux}.md
     Sphinx (--sphinx): <repo>/docs-site/reference/{index,algebras,
       root-system,representations,young-tableaux}.md

   Self-checks print ok / FAIL and exit 1 on any mismatch. Run from anywhere:

     wolframscript -file scripts/generate-reference.wls
*)

root = DirectoryName[$InputFileName, 2];
scriptDir = DirectoryName[$InputFileName];
target = If[MemberQ[$ScriptCommandLine, "--sphinx"], "Sphinx", "MkDocs"];
PacletDirectoryLoad[root];
Needs["ClassicalLieAlgebra`"];

outRoot = If[target === "Sphinx", "docs-site", "site-src"];
outDir = FileNameJoin[{root, outRoot, "reference"}];
If[! DirectoryQ[outDir], CreateDirectory[outDir, CreateIntermediateDirectories -> True]];

details = Get[FileNameJoin[{scriptDir, "reference-details.wl"}]];

shortName[name_] := Last[StringSplit[name, "`"]];
(* MessageName is HoldFirst, so look the usage up by parsing the qualified
   name -- always against the paclet context, to avoid any path shadowing. *)
usageOf[name_] := Quiet@ToExpression["ClassicalLieAlgebra`" <> shortName[name] <> "::usage"];

(* Neutralize Markdown/arithmatex specials in code-sourced usage prose (the
   fallback when a symbol has no curated description). Brackets become HTML
   entities, not backslash escapes: arithmatex reads \[ ... \] as display math,
   so escaping "[" as "\[" would turn signatures into math. Curated prose is
   authored Markdown-safe and is not passed through this. *)
mdEscape[s_String] := StringReplace[s, {
   "\\" -> "\\\\", "`" -> "\\`", "*" -> "\\*", "_" -> "\\_",
   "[" -> "&#91;", "]" -> "&#93;", "<" -> "&lt;", ">" -> "&gt;", "|" -> "&#124;"}];

docLink[url_String] := If[target === "Sphinx", StringReplace[url, ".md" -> ""], url];
pageLink[slug_String] := If[target === "Sphinx", slug, slug <> ".md"];

(* page slug -> title, one-line intro, and the symbols it documents *)
groups = {
   <|"slug" -> "algebras", "title" -> "Algebras & generators",
     "intro" -> "Construct the classical Lie algebras and read off their generators in the standard, Cartan-Weyl, and Chevalley bases.",
     "symbols" -> {"LieAlgebra", "SU", "SO", "Sp", "Generators", "CartanWeyl", "Chevalley", "BasisTransform"}|>,
   <|"slug" -> "root-system", "title" -> "Root system",
     "intro" -> "Root-system data of an algebra: its rank and dimension, Cartan matrix, simple and positive roots, and fundamental weights.",
     "symbols" -> {"Rank", "LieAlgebraDimension", "CartanMatrix", "SimpleRoots", "PositiveRoots", "FundamentalWeights"}|>,
   <|"slug" -> "representations", "title" -> "Representations",
     "intro" -> "Build an irreducible representation from its highest weight, then read off its dimension, weight system, Casimir eigenvalue, and explicit generator matrices.",
     "symbols" -> {"Irrep", "HighestWeight", "RepresentationDimension", "RepresentationMatrices", "CasimirEigenvalue", "WeightSystem"}|>,
   <|"slug" -> "young-tableaux", "title" -> "Young tableaux",
     "intro" -> "The Young-tableau toolkit for many-body wavefunctions: tableaux and their tensor form, the Young symmetrizer, inner products, normalization, and orthogonalization of degenerate weight spaces.",
     "symbols" -> {"Tableau", "TensorTableau", "TableauForm", "TableauPermute", "TableauNormalization", "TableauOrthogonalization", "TableauDot", "ToTensor", "TensorDot", "TensorNorm", "Psi"}|>
};

(* symbol -> page slug, so "see also" can link across pages *)
symbolPage = Association @@ Flatten[Function[g, (# -> g["slug"]) & /@ g["symbols"]] /@ groups];
anchor[name_] := ToLowerCase[name];
(* a "see also" entry is either a symbol name (linked to its reference anchor) or
   a {label, url} pair linking out to a narrative page (tour, tutorials, ...) *)
seeAlsoLink[{label_String, url_String}] := "[" <> label <> "](" <> docLink[url] <> ")";
seeAlsoLink[name_] := If[KeyExistsQ[symbolPage, name],
   "[`" <> name <> "`](" <> pageLink[symbolPage[name]] <> "#" <> anchor[name] <> ")",
   "`" <> name <> "`"];

(* " middot " built from ASCII source: wolframscript reads the script as
   Latin-1, so a literal U+00B7 here would double-encode in the output. *)
dot = " " <> FromCharacterCode[183] <> " ";

(* run an example expression, recording any that error or emit messages *)
exErrors = {};
runExample[in_String] := Module[{res, bad = False},
   res = Check[ToExpression[in], bad = True; $Failed];
   If[bad, AppendTo[exErrors, in]];
   res
];
fmtOut[res_] := ToString[res, InputForm, PageWidth -> 76];

rawHtml[s_String] := "```{raw} html\n" <> s <> "\n```\n\n";
sphinxInputCell[n_Integer, expr_String] := StringJoin[
   rawHtml["<div class=\"cla-example\"><div class=\"cla-example-row cla-example-in\"><div class=\"cla-example-label\">In[" <> ToString[n] <> "]:=</div><div class=\"cla-example-cell cla-example-input\">"],
   "```mathematica\n", expr, "\n```\n\n",
   rawHtml["</div></div>"]];
sphinxOutputOpen[n_Integer] := rawHtml["<div class=\"cla-example-row cla-example-out\"><div class=\"cla-example-label\">Out[" <> ToString[n] <> "]=</div><div class=\"cla-example-cell cla-example-output\">"];
sphinxOutputClose[] := rawHtml["</div></div></div>"];
sphinxNoteOpen[] := rawHtml["<div class=\"cla-example-row cla-example-out\"><div class=\"cla-example-label\"></div><div class=\"cla-example-cell\"><div class=\"cla-example-note\">"];
sphinxNoteClose[] := rawHtml["</div></div></div></div>"];

renderSphinxExample[ex_, n_Integer] := Which[
   ListQ[ex] && Length[ex] >= 3 && ex[[2]] === "image",
     (runExample[ex[[1]]];
      StringJoin[
         sphinxInputCell[n, ex[[1]]],
         sphinxOutputOpen[n],
         "![", If[Length[ex] >= 4, ex[[4]], "output"], "](", docLink[ex[[3]]], ")\n\n",
         sphinxOutputClose[]]),
   ListQ[ex],
     (runExample[ex[[1]]];
      StringJoin[
         sphinxInputCell[n, ex[[1]]],
         sphinxNoteOpen[], ex[[2]], sphinxNoteClose[]]),
   True,
     With[{res = runExample[ex]},
        StringJoin[
           sphinxInputCell[n, ex],
           sphinxOutputOpen[n],
           "```text\n", fmtOut[res], "\n```\n\n",
           sphinxOutputClose[]]]];

(* an example is one of:
     "expr"                       run it, embed the output;
     {"expr", "note"}             run to validate, show the italic note instead
                                  of output (for front-end display forms);
     {"expr", "image", "path"}    run to validate, embed an image at path (a
                                  committed rasterization, e.g. from MakeFigures.wls). *)
renderExample[ex_, n_: 1] := If[target === "Sphinx", renderSphinxExample[ex, n],
  Which[
   ListQ[ex] && Length[ex] >= 3 && ex[[2]] === "image",
     (runExample[ex[[1]]];
      StringJoin["```mathematica\n", ex[[1]], "\n```\n\n",
         "![", If[Length[ex] >= 4, ex[[4]], "output"], "](", ex[[3]], "){ .center }\n\n"]),
   ListQ[ex],
     (runExample[ex[[1]]];
      StringJoin["```mathematica\n", ex[[1]], "\n```\n\n", "*", ex[[2]], "*\n\n"]),
   True,
     With[{res = runExample[ex]},
        StringJoin["```mathematica\n", ex, "\n```\n\n", "```\n", fmtOut[res], "\n```\n\n"]]]];

argRow[{a_, d_}] := "| `" <> a <> "` | " <> d <> " |\n";
optRow[{o_, def_, d_}] := "| `" <> o <> "` | `" <> def <> "` | " <> d <> " |\n";
msgRow[{m_, w_}] := "| `" <> m <> "` | " <> w <> " |\n";

renderSymbol[name_] := Module[
   {d = Lookup[details, name, <||>], u = usageOf[name], sigs, desc, args, opts, msgs, exs, notes, see},
   sigs = Lookup[d, "signatures", {name}];
   desc = Lookup[d, "description", If[StringQ[u], mdEscape[u], ""]];
   args = Lookup[d, "args", {}];
   opts = Lookup[d, "options", {}];
   msgs = Lookup[d, "messages", {}];
   exs = Lookup[d, "examples", {}];
   notes = Lookup[d, "notes", {}];
   see = Lookup[d, "seealso", {}];
   StringJoin[
      "### ", name, "\n\n",
      StringRiffle[("`" <> # <> "`") & /@ sigs, "<br>\n"], "\n\n",
      desc, "\n\n",
      If[args =!= {}, "**Arguments**\n\n| Argument | Description |\n| :-- | :-- |\n" <> StringJoin[argRow /@ args] <> "\n", ""],
      If[opts =!= {}, "**Options**\n\n| Option | Default | Description |\n| :-- | :-- | :-- |\n" <> StringJoin[optRow /@ opts] <> "\n", ""],
      If[msgs =!= {}, "**Messages**\n\nEmits a message and returns `$Failed` on invalid input.\n\n| Message | When |\n| :-- | :-- |\n" <> StringJoin[msgRow /@ msgs] <> "\n", ""],
      If[exs =!= {}, "**Examples**\n\n" <> StringJoin[MapIndexed[renderExample[#1, First[#2]] &, exs]], ""],
      If[notes =!= {}, "**Notes**\n\n" <> StringJoin[("- " <> # <> "\n") & /@ notes] <> "\n", ""],
      If[see =!= {}, "**See also:** " <> StringRiffle[seeAlsoLink /@ see, dot] <> "\n", ""]
   ]
];

header = "<!-- Generated by scripts/generate-reference.wls. Do not edit by hand. -->\n\n";

renderPage[g_] := StringJoin[
   header, "# ", g["title"], "\n\n", g["intro"], "\n\n",
   StringRiffle[renderSymbol /@ g["symbols"], "\n---\n\n"]
];

(* per-group page files (this is where examples get evaluated) *)
Scan[Export[FileNameJoin[{outDir, #["slug"] <> ".md"}], renderPage[#], "Text", CharacterEncoding -> "UTF-8"] &, groups];

(* the reference index: the enriched overview, then the grouped symbol index *)
overviewFile = FileNameJoin[{scriptDir, "reference-overview.md"}];
intro = If[FileExistsQ[overviewFile], Import[overviewFile, "Text", CharacterEncoding -> "UTF-8"], "# API reference"];
groupIndexBlock[g_] := StringJoin[
   "### [", g["title"], "](", pageLink[g["slug"]], ")\n\n",
   g["intro"], "\n\n",
   StringRiffle[("[`" <> # <> "`](" <> pageLink[g["slug"]] <> "#" <> anchor[#] <> ")") & /@ g["symbols"], dot], "\n"];
indexContent = StringJoin[
   header, intro, "\n\n## Symbols by theme\n\n",
   If[target === "Sphinx",
      "```{toctree}\n:hidden:\n:maxdepth: 1\n\nalgebras\nroot-system\nrepresentations\nyoung-tableaux\n```\n\n",
      ""],
   StringRiffle[groupIndexBlock /@ groups, "\n"]];
Export[FileNameJoin[{outDir, "index.md"}], indexContent, "Text", CharacterEncoding -> "UTF-8"];

(* --- self-checks --------------------------------------------------------- *)
declared = Join @@ (#["symbols"] & /@ groups);
present = shortName /@ Select[Names["ClassicalLieAlgebra`*"], StringQ[usageOf[#]] &];

check[label_, ok_] := (Print[If[ok, "[ok]   ", "[FAIL] "], label]; ok);

dupes = Cases[Tally[declared], {s_, n_} /; n > 1 :> s];
missingUsage = Select[declared, ! StringQ[Quiet@usageOf[#]] &];
ungrouped = Complement[present, declared];
phantom = Complement[declared, present];
missingDetail = Complement[declared, Keys[details]];

results = {
   check["all " <> ToString[Length[declared]] <> " grouped symbols have a usage message", missingUsage === {}],
   check["no symbol appears in more than one group", dupes === {}],
   check["every public symbol (" <> ToString[Length[present]] <> ") is grouped", ungrouped === {}],
   check["no grouped symbol is missing from the paclet", phantom === {}],
   check["every grouped symbol has a curated detail entry", missingDetail === {}],
   check["every example evaluates cleanly", exErrors === {}]
};
If[missingUsage =!= {}, Print["       missing usage: ", missingUsage]];
If[dupes =!= {}, Print["       duplicated: ", dupes]];
If[ungrouped =!= {}, Print["       ungrouped: ", ungrouped]];
If[phantom =!= {}, Print["       not in paclet: ", phantom]];
If[missingDetail =!= {}, Print["       missing detail: ", missingDetail]];
If[exErrors =!= {}, Print["       example errors:"]; Scan[Print["         ", #] &, exErrors]];

Print["wrote ", Length[groups] + 1, " files to ", outDir];
Exit[If[AllTrue[results, TrueQ], 0, 1]]
