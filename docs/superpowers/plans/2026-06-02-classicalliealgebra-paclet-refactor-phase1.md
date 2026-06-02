# ClassicalLieAlgebra Phase 1 — paclet refactor + API redesign — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Turn the duplicated, monolithic `ClassicalLieAlgebra.wl` into a conventional Wolfram **paclet** with a redesigned API (canonical `LieAlgebra[type,rank]` + `SU/SO/Sp` sugar, a root-system data layer, one unified `Generators[g,scheme]` basis interface), robust error handling, automated `.wlt` tests, and CI — preserving the current matrix/tableau math.

**Architecture:** A single public context `ClassicalLieAlgebra`` whose main `Kernel/` file declares all usage/messages and `Get`s feature subfiles; shared private helpers live in a `ClassicalLieAlgebra`Common`` sub-context. Algebra objects canonicalize to `LieAlgebra[type,rank]`; matrix realizations and Young-tableau machinery are migrated from the current code with `Block`→`Module` and the `` `Private` `` backtick fix. Tests are `VerificationTest` `.wlt` files run headless with an `Exit[1]`-on-failure runner; CI uses the free Wolfram Engine container.

**Tech Stack:** Wolfram Language (paclet, `BeginPackage`/`Needs`, `VerificationTest`/`TestReport`), GitHub Actions + `wolframresearch/wolframengine` image.

---

## Conventions for the implementer

- **Running Wolfram on this machine:** the system `wolframscript` is broken; use the app-bundled one. Define for your shell session:
  `WS=/Applications/Wolfram.app/Contents/MacOS/wolframscript`
  Run a single test file: `"$WS" -code 'r=TestReport["Tests/Algebras.wlt"]; Print[r["AllTestsSucceeded"]," ",r["TestsFailedCount"]]'`
  The committed `scripts/runTests.wls` and CI use plain `wolframscript` (correct inside the engine container).
- **Dev load loop:** `"$WS" -code 'PacletDirectoryLoad["'"$PWD"'"]; Needs["ClassicalLieAlgebra`"]; <expr>'`. Editing source then re-running a fresh `wolframscript` always loads current source (new process), so no reload dance is needed for test runs.
- **TDD:** every task writes the `.wlt` test first, runs it red, implements, runs it green, commits. A test "fails red" if `AllTestsSucceeded` is `False` (or the file errors on load).
- **Commits:** conventional-commit style, on `master` only with the user's say-so — otherwise commit to a `paclet-refactor` branch. End commit messages with the Co-Authored-By trailer.
- **Each `.wlt` begins with** loading the package:
  ```wolfram
  PacletDirectoryLoad[DirectoryName[$InputFileName, 2]]; (* repo root from Tests/X.wlt *)
  Needs["ClassicalLieAlgebra`"];
  ```
  (When run via `TestReport[file]`, `$InputFileName` is the file; `DirectoryName[#,2]` is the repo root.)

---

## File structure

| File | Responsibility |
| :--- | :--- |
| `PacletInfo.wl` | paclet metadata + `Kernel` extension |
| `Kernel/ClassicalLieAlgebra.wl` | public context: all `::usage` + `::tag` messages, `SyntaxInformation`, `Get` subfiles, `Protect`, `EndPackage` |
| `Kernel/Common.wl` | `ClassicalLieAlgebra`Common`` sub-context: matrix-unit builder, shared validation predicates |
| `Kernel/Algebras.wl` | canonical `LieAlgebra[]`, `SU/SO/Sp` sugar, display, validation, root-system data (`Rank`, `LieAlgebraDimension`, `CartanMatrix`, `SimpleRoots`, `PositiveRoots`, `FundamentalWeights`) |
| `Kernel/SpecialUnitary.wl` | su(n) defining/Cartan–Weyl/Chevalley matrices → wired into `Generators[g,scheme]` |
| `Kernel/SpecialOrthogonal.wl` | so(n) matrices + `"Realization"` + `BasisTransform` |
| `Kernel/Symplectic.wl` | sp(2n) matrices + `"Realization"` + `BasisTransform` |
| `Kernel/YoungTableaux.wl` | `Tableau`/`TensorTableau`/`Psi`/… (robust linear combinations) |
| `Tests/*.wlt` | one per feature module |
| `scripts/runTests.wls` | run all `.wlt`, `Exit[1]` on failure |
| `.github/workflows/test.yml` | CI |

Deleted at the end: `ClassicalLieAlgebra.wl` (old monolith, root), `SpecialUnitary/`, `SpecialOrthogonal/`, `SymplecticUnitary/`, `YoungTableau/`.

---

## Task 1: Paclet scaffold loads

**Files:**
- Create: `PacletInfo.wl`, `Kernel/ClassicalLieAlgebra.wl`, `Kernel/Common.wl`
- Test: `Tests/Loading.wlt`

- [ ] **Step 1: Write the failing test** — `Tests/Loading.wlt`

```wolfram
PacletDirectoryLoad[DirectoryName[$InputFileName, 2]];
Needs["ClassicalLieAlgebra`"];

VerificationTest[ MemberQ[$Packages, "ClassicalLieAlgebra`"], True, TestID -> "package-loads" ];
VerificationTest[ Head[ClassicalLieAlgebra`Common`matrixUnit[2, 1, 2]], List, TestID -> "common-helper-visible-internally" ];
```

- [ ] **Step 2: Run, verify red**

Run: `"$WS" -code 'r=TestReport["Tests/Loading.wlt"]; Print[r["AllTestsSucceeded"]]'`
Expected: `False` (package not found / errors).

- [ ] **Step 3: Create `PacletInfo.wl`**

```wolfram
PacletObject[<|
  "Name" -> "ClassicalLieAlgebra",
  "Version" -> "1.0.0",
  "WolframVersion" -> "13.0+",
  "Description" -> "Generators, bases, and Young-tableau machinery for the classical Lie algebras.",
  "Creator" -> "Jie Ren",
  "License" -> "MIT",
  "SourceControlURL" -> "https://github.com/jayren3996/LieAlgebra",
  "Extensions" -> {
    {"Kernel", "Root" -> "Kernel", "Context" -> "ClassicalLieAlgebra`"}
  }
|>]
```

- [ ] **Step 4: Create `Kernel/Common.wl`**

```wolfram
BeginPackage["ClassicalLieAlgebra`Common`"];

matrixUnit;       (* matrixUnit[n,i,j] = n x n matrix with 1 at (i,j) *)
classicalTypeQ;   (* classicalTypeQ["A"] etc. *)

Begin["`Private`"];

matrixUnit[n_Integer, i_Integer, j_Integer] := SparseArray[{{i, j} -> 1}, {n, n}] // Normal;
classicalTypeQ[t_] := MemberQ[{"A", "B", "C", "D"}, t];

End[];
EndPackage[];
```

- [ ] **Step 5: Create `Kernel/ClassicalLieAlgebra.wl`** (shell that loads Common; subfiles added in later tasks)

```wolfram
BeginPackage["ClassicalLieAlgebra`"];

(* Public symbols + usage are filled in across later tasks; declared here. *)

Begin["`Private`"];
Needs["ClassicalLieAlgebra`Common`"];
(* << subfiles added in later tasks *)
End[];

EndPackage[];
```

- [ ] **Step 6: Run, verify green; commit**

Run: `"$WS" -code 'r=TestReport["Tests/Loading.wlt"]; Print[r["AllTestsSucceeded"]]'` → `True`
```bash
git add PacletInfo.wl Kernel/ Tests/Loading.wlt
git commit -m "feat(paclet): scaffold ClassicalLieAlgebra paclet that loads"
```

---

## Task 2: Canonical algebra objects + SU/SO/Sp sugar + display

**Files:**
- Create: `Kernel/Algebras.wl`; Modify: `Kernel/ClassicalLieAlgebra.wl` (add usage/messages, `<< Algebras`)
- Test: `Tests/Algebras.wlt`

- [ ] **Step 1: Write the failing test** — `Tests/Algebras.wlt`

```wolfram
PacletDirectoryLoad[DirectoryName[$InputFileName, 2]];
Needs["ClassicalLieAlgebra`"];

VerificationTest[ SU[3], LieAlgebra["A", 2], TestID -> "SU-sugar" ];
VerificationTest[ Sp[4], LieAlgebra["C", 2], TestID -> "Sp-sugar" ];
VerificationTest[ SO[5], LieAlgebra["B", 2], TestID -> "SO-odd-sugar" ];
VerificationTest[ SO[6], LieAlgebra["D", 3], TestID -> "SO-even-sugar" ];
VerificationTest[ Sp[3], $Failed, {Sp::evenrank}, TestID -> "Sp-odd-fails" ];
VerificationTest[ LieAlgebra["E", 2], $Failed, {LieAlgebra::badtype}, TestID -> "bad-type-fails" ];
VerificationTest[ LieAlgebra["A", 0], $Failed, {LieAlgebra::badrank}, TestID -> "bad-rank-fails" ];
VerificationTest[ ToString[LieAlgebra["A", 2] // StandardForm], "A\!\(\*SubscriptBox[\(\), \(2\)]\)" =!= "", TestID -> "has-display" ];
```

- [ ] **Step 2: Run, verify red**

Run: `"$WS" -code 'r=TestReport["Tests/Algebras.wlt"]; Print[r["AllTestsSucceeded"]]'` → `False`.

- [ ] **Step 3: Add public declarations to `Kernel/ClassicalLieAlgebra.wl`** (between `BeginPackage` and `Begin["`Private`"]`)

```wolfram
LieAlgebra::usage = "LieAlgebra[type, rank] is the simple Lie algebra of Cartan type \"A\"|\"B\"|\"C\"|\"D\" and given rank.";
SU::usage = "SU[n] represents su(n), the type A_{n-1} algebra (n>=2).";
SO::usage = "SO[n] represents so(n): type B for odd n>=3, type D for even n>=4.";
Sp::usage = "Sp[n] represents sp(n), the type C_{n/2} algebra (even n>=2).";

LieAlgebra::badtype = "`1` is not a valid Cartan type; use \"A\", \"B\", \"C\" or \"D\".";
LieAlgebra::badrank = "`1` is not a valid rank for type `2`.";
SU::baddim = "SU[`1`] requires an integer n>=2.";
SO::baddim = "SO[`1`] requires an integer n>=3.";
Sp::evenrank = "Sp[`1`] requires an even integer n>=2.";
```

- [ ] **Step 4: Create `Kernel/Algebras.wl`**

```wolfram
BeginPackage["ClassicalLieAlgebra`Algebras`"];
Needs["ClassicalLieAlgebra`"];
Needs["ClassicalLieAlgebra`Common`"];
Begin["`Private`"];

validRankQ["A", r_] := IntegerQ[r] && r >= 1;
validRankQ["B", r_] := IntegerQ[r] && r >= 1;
validRankQ["C", r_] := IntegerQ[r] && r >= 1;
validRankQ["D", r_] := IntegerQ[r] && r >= 2;
validRankQ[_, _] := False;

(* validation downvalues fire ONLY for invalid input; valid LieAlgebra[t,r] stays inert *)
LieAlgebra[t_, r_] /; (! classicalTypeQ[t]) := (Message[LieAlgebra::badtype, t]; $Failed);
LieAlgebra[t_, r_] /; (classicalTypeQ[t] && ! validRankQ[t, r]) := (Message[LieAlgebra::badrank, r, t]; $Failed);

SU[n_Integer] /; n >= 2 := LieAlgebra["A", n - 1];
SU[n_] := (Message[SU::baddim, n]; $Failed);
SO[n_Integer] /; (n >= 3 && OddQ[n]) := LieAlgebra["B", (n - 1)/2];
SO[n_Integer] /; (n >= 4 && EvenQ[n]) := LieAlgebra["D", n/2];
SO[n_] := (Message[SO::baddim, n]; $Failed);
Sp[n_Integer] /; (n >= 2 && EvenQ[n]) := LieAlgebra["C", n/2];
Sp[n_] := (Message[Sp::evenrank, n]; $Failed);

LieAlgebra /: MakeBoxes[LieAlgebra[t_String, r_Integer], StandardForm] :=
  SubscriptBox[ToString[t], ToString[r]];

End[];
EndPackage[];
```

- [ ] **Step 5: Wire `<< Algebras` into the main file** — in `Kernel/ClassicalLieAlgebra.wl`, inside the `Begin["`Private`"]` block, after the `Needs["ClassicalLieAlgebra`Common`"]` line add:
```wolfram
Get["ClassicalLieAlgebra`Algebras`"];
```

- [ ] **Step 6: Run green; commit**

Run: `"$WS" -code 'r=TestReport["Tests/Algebras.wlt"]; Print[r["AllTestsSucceeded"]," ",r["TestsFailedCount"]]'` → `True 0`
```bash
git add Kernel/ Tests/Algebras.wlt
git commit -m "feat(algebra): canonical LieAlgebra objects with SU/SO/Sp sugar, validation, display"
```

---

## Task 3: Root-system data layer

**Files:** Modify `Kernel/ClassicalLieAlgebra.wl` (usage), `Kernel/Algebras.wl` (definitions); Test: append to `Tests/Algebras.wlt`.

Math facts used (Bourbaki/Humphreys, simple roots in the Euclidean `e_i` basis):
- `SimpleRoots`: A_r → `e_i - e_{i+1}` (i=1..r) in `R^{r+1}`; B_r → `e_i-e_{i+1}` (i<r) and `e_r` in `R^r`; C_r → `e_i-e_{i+1}` (i<r) and `2 e_r`; D_r → `e_i-e_{i+1}` (i<r) and `e_{r-1}+e_r`.
- `CartanMatrix[g][[i,j]] = 2 (a_i·a_j)/(a_j·a_j)`.
- `FundamentalWeights = Inverse[CartanMatrix].SimpleRoots` (rows = ω_i in the `e_i` basis).
- `LieAlgebraDimension`: A_r → `r(r+2)`, B_r → `r(2r+1)`, C_r → `r(2r+1)`, D_r → `r(2r-1)`.
- `PositiveRoots`: A_r → `e_i-e_j (i<j)`; B_r → `{e_i±e_j (i<j)} ∪ {e_i}`; C_r → `{e_i±e_j (i<j)} ∪ {2 e_i}`; D_r → `{e_i±e_j (i<j)}`.

- [ ] **Step 1: Write failing tests** — append to `Tests/Algebras.wlt`

```wolfram
VerificationTest[ Rank[SU[4]], 3, TestID -> "rank-A3" ];
VerificationTest[ LieAlgebraDimension[SU[3]], 8, TestID -> "dim-su3" ];
VerificationTest[ LieAlgebraDimension[SO[5]], 10, TestID -> "dim-so5" ];
VerificationTest[ LieAlgebraDimension[Sp[4]], 10, TestID -> "dim-sp4" ];
VerificationTest[ CartanMatrix[LieAlgebra["A", 2]], {{2, -1}, {-1, 2}}, TestID -> "cartan-A2" ];
VerificationTest[ CartanMatrix[LieAlgebra["B", 2]], {{2, -2}, {-1, 2}}, TestID -> "cartan-B2" ];
VerificationTest[ CartanMatrix[LieAlgebra["C", 2]], {{2, -1}, {-2, 2}}, TestID -> "cartan-C2" ];
VerificationTest[ CartanMatrix[LieAlgebra["D", 4]], {{2,-1,0,0},{-1,2,-1,-1},{0,-1,2,0},{0,-1,0,2}}, TestID -> "cartan-D4" ];
VerificationTest[ Length[PositiveRoots[LieAlgebra["A", 2]]], 3, TestID -> "posroots-A2-count" ];
VerificationTest[ FundamentalWeights[LieAlgebra["A", 2]], {{2/3, -1/3, -1/3}, {1/3, 1/3, -2/3}}, TestID -> "fundweights-A2" ];
```

- [ ] **Step 2: Run, verify red** (`AllTestsSucceeded` → `False`).

- [ ] **Step 3: Add usage messages** to `Kernel/ClassicalLieAlgebra.wl`:

```wolfram
Rank::usage = "Rank[g] gives the rank of the Lie algebra g.";
LieAlgebraDimension::usage = "LieAlgebraDimension[g] gives the dimension of the Lie algebra g.";
CartanMatrix::usage = "CartanMatrix[g] gives the Cartan matrix of g.";
SimpleRoots::usage = "SimpleRoots[g] gives the simple roots of g in the Euclidean basis.";
PositiveRoots::usage = "PositiveRoots[g] gives the positive roots of g in the Euclidean basis.";
FundamentalWeights::usage = "FundamentalWeights[g] gives the fundamental weights of g in the Euclidean basis.";
```

- [ ] **Step 4: Implement** in `Kernel/Algebras.wl` (inside `` `Private` ``)

```wolfram
Rank[LieAlgebra[_, r_]] := r;

LieAlgebraDimension[LieAlgebra["A", r_]] := r (r + 2);
LieAlgebraDimension[LieAlgebra["B", r_]] := r (2 r + 1);
LieAlgebraDimension[LieAlgebra["C", r_]] := r (2 r + 1);
LieAlgebraDimension[LieAlgebra["D", r_]] := r (2 r - 1);

ee[d_, i_] := UnitVector[d, i];
SimpleRoots[LieAlgebra["A", r_]] := Table[ee[r + 1, i] - ee[r + 1, i + 1], {i, r}];
SimpleRoots[LieAlgebra["B", r_]] := Append[Table[ee[r, i] - ee[r, i + 1], {i, r - 1}], ee[r, r]];
SimpleRoots[LieAlgebra["C", r_]] := Append[Table[ee[r, i] - ee[r, i + 1], {i, r - 1}], 2 ee[r, r]];
SimpleRoots[LieAlgebra["D", r_]] := Append[Table[ee[r, i] - ee[r, i + 1], {i, r - 1}], ee[r, r - 1] + ee[r, r]];

CartanMatrix[g : LieAlgebra[_, _]] := Module[{a = SimpleRoots[g]},
  Table[2 (a[[i]] . a[[j]])/(a[[j]] . a[[j]]), {i, Length@a}, {j, Length@a}]];

FundamentalWeights[g : LieAlgebra[_, _]] := Inverse[CartanMatrix[g]] . SimpleRoots[g];

PositiveRoots[LieAlgebra["A", r_]] := Module[{d = r + 1},
  Flatten[Table[ee[d, i] - ee[d, j], {i, d}, {j, i + 1, d}], 1]];
PositiveRoots[LieAlgebra["B", r_]] := Join[
  Flatten[Table[ee[r, i] - ee[r, j], {i, r}, {j, i + 1, r}], 1],
  Flatten[Table[ee[r, i] + ee[r, j], {i, r}, {j, i + 1, r}], 1],
  Table[ee[r, i], {i, r}]];
PositiveRoots[LieAlgebra["C", r_]] := Join[
  Flatten[Table[ee[r, i] - ee[r, j], {i, r}, {j, i + 1, r}], 1],
  Flatten[Table[ee[r, i] + ee[r, j], {i, r}, {j, i + 1, r}], 1],
  Table[2 ee[r, i], {i, r}]];
PositiveRoots[LieAlgebra["D", r_]] := Join[
  Flatten[Table[ee[r, i] - ee[r, j], {i, r}, {j, i + 1, r}], 1],
  Flatten[Table[ee[r, i] + ee[r, j], {i, r}, {j, i + 1, r}], 1]];
```

- [ ] **Step 5: Run green; commit**

Run: `"$WS" -code 'r=TestReport["Tests/Algebras.wlt"]; Print[r["AllTestsSucceeded"]," ",r["TestsFailedCount"]]'` → `True 0`
```bash
git add Kernel/ Tests/Algebras.wlt
git commit -m "feat(algebra): root-system data (rank, dim, Cartan matrix, roots, weights)"
```

---

## Task 4: su(n) matrices wired to `Generators[g, scheme]`

Migrate the SU constructors from the current `ClassicalLieAlgebra.wl:76-143` (`SUT`, `SUCWH`, `SUCWE`, `SUCWF`, `SUCW`, `SUH`, `SUE`, `SUF`, `SUn`) into `Kernel/SpecialUnitary.wl`. **Transformation rules applied to every migrated function:** rename to lowerCamelCase private (`suStandard`, `suCartanWeyl`, `suChevalley`, …); `Block`→`Module`; **declare every local including loop counters** (fixes the `j` leak at current `:106-107`); use `matrixUnit` from `Common` where a `{i,j}->1` matrix is built.

**Files:** Create `Kernel/SpecialUnitary.wl`; Modify main file (usage for `Generators`; `<< SpecialUnitary`); Test `Tests/SpecialUnitary.wlt`.

- [ ] **Step 1: Write failing test** — `Tests/SpecialUnitary.wlt`

```wolfram
PacletDirectoryLoad[DirectoryName[$InputFileName, 2]];
Needs["ClassicalLieAlgebra`"];

(* su(2): 3 generators, each 2x2 traceless Hermitian *)
VerificationTest[ Length[Generators[SU[2]]], 3, TestID -> "su2-gen-count" ];
VerificationTest[ AllTrue[Generators[SU[2]], # == ConjugateTranspose[#] &], True, TestID -> "su2-hermitian" ];
VerificationTest[ AllTrue[Generators[SU[2]], Tr[#] == 0 &], True, TestID -> "su2-traceless" ];

(* su(3) defining rep: 8 generators *)
VerificationTest[ Length[Generators[SU[3]]], 8, TestID -> "su3-gen-count" ];

(* Chevalley association shape + diagonal Cartan with the expected H1 *)
VerificationTest[ Keys[Generators[SU[3], "Chevalley"]], {"Cartan", "Raising", "Lowering"}, TestID -> "cheval-assoc-keys" ];
VerificationTest[ Diagonal[Generators[SU[3], "Chevalley"]["Cartan"][[1]]], {1, -1, 0}, TestID -> "su3-chevalley-H1" ];
VerificationTest[ Generators[SU[3], "Chevalley"]["Lowering"][[1]], Transpose[Generators[SU[3], "Chevalley"]["Raising"][[1]]], TestID -> "su3-F-is-Etranspose" ];

(* Cartan–Weyl association shape *)
VerificationTest[ Length[Generators[SU[3], "CartanWeyl"]["Raising"]], 3, TestID -> "su3-cw-raising-count" ];

(* bad scheme *)
VerificationTest[ Generators[SU[3], "Nope"], $Failed, {Generators::badscheme}, TestID -> "bad-scheme" ];
```

- [ ] **Step 2: Run red.**

- [ ] **Step 3: Add to main file** — usage + message + the dispatcher signature lives here so all algebra files extend one symbol:

```wolfram
Generators::usage = "Generators[g] gives the defining-representation generators of g. Generators[g,\"CartanWeyl\"] and Generators[g,\"Chevalley\"] give associations <|\"Cartan\"->..,\"Raising\"->..,\"Lowering\"->..|>. Option \"Realization\"->\"Diagonal\"|\"Antisymmetric\" applies to so/sp.";
Generators::badscheme = "`1` is not a valid scheme; use \"Standard\", \"CartanWeyl\" or \"Chevalley\".";
Generators::badrealization = "`1` is not a valid \"Realization\"; use \"Diagonal\" or \"Antisymmetric\".";
Options[Generators] = {"Realization" -> "Diagonal"};
SyntaxInformation[Generators] = {"ArgumentsPattern" -> {_, _., OptionsPattern[]}};
```

- [ ] **Step 4: Create `Kernel/SpecialUnitary.wl`** — migrate constructors and add the SU dispatch downvalues. Skeleton (migrate bodies from `ClassicalLieAlgebra.wl:76-143` applying the transformation rules above):

```wolfram
BeginPackage["ClassicalLieAlgebra`SpecialUnitary`"];
Needs["ClassicalLieAlgebra`"];
Needs["ClassicalLieAlgebra`Common`"];
Begin["`Private`"];

(* --- migrated, Module-ized, locals declared (example shows the pattern) --- *)
suStandard[n_] := Module[{list = {}, i, j},
  For[j = 2, j <= n, j++,
    For[i = 1, i < j, i++,
      AppendTo[list, suT[{n, 1}, {i, j}]];
      AppendTo[list, suT[{n, 2}, {i, j}]];
    ];
    AppendTo[list, suT[{n, 3}, {j}]];
  ]; list];
(* suT[{n,1|2|3},..], suCartanWeyl[n], suChevalley[n] migrated verbatim from
   ClassicalLieAlgebra.wl:77-143 with Block->Module and all locals (incl. the
   leaked `j` in SUCWH:106) declared. suChevalley returns {Hl,El,Fl}. *)

(* --- dispatch: g already canonical LieAlgebra["A", r] (SU sugar resolved) --- *)
toAssoc[{h_, e_, f_}] := <|"Cartan" -> h, "Raising" -> e, "Lowering" -> f|>;

Generators[LieAlgebra["A", r_]] := suStandard[r + 1];
Generators[LieAlgebra["A", r_], "Standard"] := suStandard[r + 1];
Generators[LieAlgebra["A", r_], "CartanWeyl", OptionsPattern[]] := toAssoc[suCartanWeyl[r + 1]];
Generators[LieAlgebra["A", r_], "Chevalley", OptionsPattern[]] := toAssoc[suChevalley[r + 1]];

End[];
EndPackage[];
```

Add the generic bad-scheme fallthrough **in the main file's Private block, after all subfiles load** (so it is least-specific): `Generators[_LieAlgebra, s_, OptionsPattern[]] := (Message[Generators::badscheme, s]; $Failed);`

- [ ] **Step 5: Wire `<< SpecialUnitary` into main file**; run green; commit.

Run: `"$WS" -code 'r=TestReport["Tests/SpecialUnitary.wlt"]; Print[r["AllTestsSucceeded"]," ",r["TestsFailedCount"]]'` → `True 0`
```bash
git add Kernel/ Tests/SpecialUnitary.wlt
git commit -m "feat(su): su(n) generators via unified Generators[g,scheme]"
```

---

## Task 5: so(n) matrices + "Realization" + BasisTransform

Migrate SO constructors from `ClassicalLieAlgebra.wl:152-363` (`SOT`, `SOCW*`, `SOCh*`, `SOBasis`, `SO[H/E/F]`, `SOn`) into `Kernel/SpecialOrthogonal.wl`, same transformation rules. Map: `"Chevalley","Realization"->"Diagonal"` ← `SOn` (diagonal); `"Realization"->"Antisymmetric"` ← `SOCh`; `"CartanWeyl"` ← `SOCW`; `BasisTransform` ← `SOBasis`. Canonical input: B_r ↔ matrix dim `2r+1`, D_r ↔ `2r`.

**Files:** Create `Kernel/SpecialOrthogonal.wl`; Modify main (usage `BasisTransform`, `<< SpecialOrthogonal`); Test `Tests/SpecialOrthogonal.wlt`.

- [ ] **Step 1: Write failing test** — `Tests/SpecialOrthogonal.wlt`

```wolfram
PacletDirectoryLoad[DirectoryName[$InputFileName, 2]];
Needs["ClassicalLieAlgebra`"];

dim[g_] := Length[First[Generators[g]]];
bracket[a_, b_] := a.b - b.a;
closesQ[gens_] := Module[{sp = Flatten[#] & /@ gens},
  AllTrue[Tuples[Range[Length@gens], 2],
    MatrixRank[Append[sp, Flatten[bracket[gens[[#[[1]]]], gens[[#[[2]]]]]]]] == MatrixRank[sp] &]];

VerificationTest[ dim[SO[5]], 5, TestID -> "so5-matdim" ];      (* B2 -> 5x5 *)
VerificationTest[ dim[SO[6]], 6, TestID -> "so6-matdim" ];      (* D3 -> 6x6 *)
VerificationTest[ Length[Generators[SO[5]]], 10, TestID -> "so5-gen-count" ];  (* dim so(5)=10 *)
VerificationTest[ closesQ[Generators[SO[5]]], True, TestID -> "so5-closes-under-bracket" ];

(* Realization: diagonal Cartan really diagonal; antisymmetric differs *)
VerificationTest[ DiagonalMatrixQ[Generators[SO[6], "Chevalley"]["Cartan"][[1]]], True, TestID -> "so6-diag-cartan-default" ];
VerificationTest[ DiagonalMatrixQ[Generators[SO[6], "Chevalley", "Realization" -> "Antisymmetric"]["Cartan"][[1]]], False, TestID -> "so6-antisym-not-diag" ];

(* BasisTransform conjugates antisymmetric Chevalley into the diagonal one *)
VerificationTest[
  With[{u = BasisTransform[SO[6]],
        d = Generators[SO[6], "Chevalley"]["Raising"][[1]],
        a = Generators[SO[6], "Chevalley", "Realization" -> "Antisymmetric"]["Raising"][[1]]},
    Chop[u.a.ConjugateTranspose[u] - d] == 0 ConstantArray[0, Dimensions[d]] // (Norm[Flatten[#]] == 0 &)],
  True, TestID -> "so6-basistransform-conjugation" ];

VerificationTest[ Generators[SO[6], "Chevalley", "Realization" -> "Nope"], $Failed, {Generators::badrealization}, TestID -> "so-bad-realization" ];
```

- [ ] **Step 2: Run red.**

- [ ] **Step 3: Add usage to main file**
```wolfram
BasisTransform::usage = "BasisTransform[g] gives the matrix conjugating the antisymmetric realization of so/sp into the diagonal one (identity for su).";
```

- [ ] **Step 4: Create `Kernel/SpecialOrthogonal.wl`** — migrate the bodies; add downvalues. The `"Realization"` option is read with `OptionValue`, validated, and selects the migrated constructor:

```wolfram
BeginPackage["ClassicalLieAlgebra`SpecialOrthogonal`"];
Needs["ClassicalLieAlgebra`"];
Needs["ClassicalLieAlgebra`Common`"];
Begin["`Private`"];

(* soStandard[m], soCartanWeyl[m], soChevalleyDiagonal[m] (<-SOn),
   soChevalleyAntisym[m] (<-SOCh), soBasis[m] (<-SOBasis) migrated from
   ClassicalLieAlgebra.wl:152-363, Block->Module, all locals declared
   (incl. leaked `i` near :281-305). m = matrix dimension. *)

soMatDim[LieAlgebra["B", r_]] := 2 r + 1;
soMatDim[LieAlgebra["D", r_]] := 2 r;
toAssoc[{h_, e_, f_}] := <|"Cartan" -> h, "Raising" -> e, "Lowering" -> f|>;

realization[OptionsPattern[Generators]] := Module[{rz = OptionValue[Generators, {}, "Realization"]},
  If[MemberQ[{"Diagonal", "Antisymmetric"}, rz], rz, $bad]];

Generators[g : LieAlgebra["B" | "D", _]] := soStandard[soMatDim[g]];
Generators[g : LieAlgebra["B" | "D", _], "Standard"] := soStandard[soMatDim[g]];
Generators[g : LieAlgebra["B" | "D", _], "CartanWeyl", OptionsPattern[]] := toAssoc[soCartanWeyl[soMatDim[g]]];
Generators[g : LieAlgebra["B" | "D", _], "Chevalley", opts : OptionsPattern[]] :=
  Switch[realization[opts],
    "Diagonal", toAssoc[soChevalleyDiagonal[soMatDim[g]]],
    "Antisymmetric", toAssoc[soChevalleyAntisym[soMatDim[g]]],
    _, Message[Generators::badrealization, OptionValue[Generators, {opts}, "Realization"]]; $Failed];
BasisTransform[g : LieAlgebra["B" | "D", _]] := soBasis[soMatDim[g]];

End[];
EndPackage[];
```

- [ ] **Step 5: Wire `<< SpecialOrthogonal`; run green; commit.**
```bash
git add Kernel/ Tests/SpecialOrthogonal.wlt
git commit -m "feat(so): so(n) generators, realizations, BasisTransform"
```

---

## Task 6: sp(2n) matrices + "Realization" + BasisTransform

Migrate Sp constructors from `ClassicalLieAlgebra.wl:373-543` (`SpT`, `SpCW*`, `SpCh*`, `SpBasis`, `Sp[H/E/F]`, `Spn`) into `Kernel/Symplectic.wl`, same rules. Map identical to Task 5 (`spChevalleyDiagonal` ← `Spn`, `spChevalleyAntisym` ← `SpCh`, `spCartanWeyl` ← `SpCW`, `spBasis` ← `SpBasis`). Canonical input C_r ↔ matrix dim `2r`.

**Files:** Create `Kernel/Symplectic.wl`; Modify main (`<< Symplectic`); Test `Tests/Symplectic.wlt`.

- [ ] **Step 1: Write failing test** — `Tests/Symplectic.wlt` (includes the `StandardBasis[USp[2n]]` regression that motivated the last fix)

```wolfram
PacletDirectoryLoad[DirectoryName[$InputFileName, 2]];
Needs["ClassicalLieAlgebra`"];

bracket[a_, b_] := a.b - b.a;
closesQ[gens_] := Module[{sp = Flatten[#] & /@ gens},
  AllTrue[Tuples[Range[Length@gens], 2],
    MatrixRank[Append[sp, Flatten[bracket[gens[[#[[1]]]], gens[[#[[2]]]]]]]] == MatrixRank[sp] &]];

VerificationTest[ Length[First[Generators[Sp[4]]]], 4, TestID -> "sp4-matdim" ];   (* C2 -> 4x4 *)
VerificationTest[ Length[Generators[Sp[4]]], 10, TestID -> "sp4-gen-count" ];        (* dim sp(4)=10 *)
VerificationTest[ closesQ[Generators[Sp[4]]], True, TestID -> "sp4-closes" ];

(* regression: BasisTransform (old StandardBasis) for sp is invertible & conjugates correctly *)
VerificationTest[ Det[BasisTransform[Sp[4]]] != 0, True, TestID -> "sp4-basistransform-invertible" ];
VerificationTest[
  With[{u = BasisTransform[Sp[4]],
        d = Generators[Sp[4], "Chevalley"]["Raising"][[1]],
        a = Generators[Sp[4], "Chevalley", "Realization" -> "Antisymmetric"]["Raising"][[1]]},
    Norm[Flatten[Chop[u.a.Inverse[u] - d]]] == 0],
  True, TestID -> "sp4-basistransform-conjugation" ];
VerificationTest[ DiagonalMatrixQ[Generators[Sp[6], "Chevalley"]["Cartan"][[1]]], True, TestID -> "sp6-diag-cartan" ];
```

- [ ] **Step 2: Run red.**
- [ ] **Step 3: Create `Kernel/Symplectic.wl`** — same shape as Task 5's file with `sp*` constructors and `spMatDim[LieAlgebra["C", r_]] := 2 r`, downvalues on `LieAlgebra["C", _]`.
- [ ] **Step 4: Wire `<< Symplectic`; run green; commit.**
```bash
git add Kernel/ Tests/Symplectic.wlt
git commit -m "feat(sp): sp(2n) generators, realizations, BasisTransform (USp regression)"
```

---

## Task 7: Young-tableau module (robust linear combinations)

Migrate the **`YoungTableau/Tableau.wl`** version (the more advanced one) into `Kernel/YoungTableaux.wl`. Apply: `Block`→`Module` with declared locals; keep `` Begin["`Private`"] ``. **Drop** `ColumnCanonicalize` (declared, never defined). **Robustness fix:** wrap the public entry points so a scalar-times-sum distributes — change `ToTensor`, `TableauForm`, `TensorNorm`, `TableauDot`, `TableauNormalization`, `TableauOrthogonalization` to `Expand` their argument first (so `(T[a]+T[b])/Sqrt[6]` is turned into `T[a]/Sqrt[6]+T[b]/Sqrt[6]`, which the existing `t1_+t2_` and `c_*t_` rules already handle).

**Files:** Create `Kernel/YoungTableaux.wl`; Modify main (usage for the 11 symbols, `<< YoungTableaux`); Test `Tests/YoungTableaux.wlt`.

- [ ] **Step 1: Write failing test** — `Tests/YoungTableaux.wlt` (ports `YoungTableau/Test.wls` plus the README (1,1)-rep facts plus the scalar-over-sum regression)

```wolfram
PacletDirectoryLoad[DirectoryName[$InputFileName, 2]];
Needs["ClassicalLieAlgebra`"];

(* ported from YoungTableau/Test.wls *)
VerificationTest[ TableauPermute[Tableau[{{1, 2}, {3}}], Psi[1, 1, 2]],
  2 Psi[1, 1, 2] - Psi[2, 1, 1] - Psi[1, 2, 1], TestID -> "permute-112" ];
VerificationTest[ TensorNorm[Psi[1, 2, 3] - 2 Psi[3, 2, 1] + Psi[4, 5, 6] - 2 Psi[3, 2, 1]],
  Sqrt[18], TestID -> "tensornorm" ];
VerificationTest[ ToTensor[TensorTableau[{{1, 1}, {2}}]],
  2 Psi[1, 1, 2] - Psi[1, 2, 1] - Psi[2, 1, 1], TestID -> "totensor-11hw" ];

(* README (1,1)-rep orthogonalization *)
VerificationTest[
  TableauOrthogonalization[
    2 TensorTableau[{{1, 2}, {3}}] - TensorTableau[{{1, 3}, {2}}],
    TensorTableau[{{1, 3}, {2}}] + TensorTableau[{{1, 2}, {3}}]],
  {2 TensorTableau[{{1, 2}, {3}}] - TensorTableau[{{1, 3}, {2}}], (3/2) TensorTableau[{{1, 3}, {2}}]},
  TestID -> "orthogonalization" ];

(* robustness: scalar times a sum must distribute, not stay unevaluated *)
VerificationTest[
  ToTensor[(TensorTableau[{{1, 2}, {3}}] + TensorTableau[{{1, 3}, {2}}])/Sqrt[6]],
  ToTensor[TensorTableau[{{1, 2}, {3}}]/Sqrt[6] + TensorTableau[{{1, 3}, {2}}]/Sqrt[6]],
  TestID -> "scalar-over-sum-distributes" ];
VerificationTest[ Head[TableauForm[(TensorTableau[{{1,2},{3}}]+TensorTableau[{{1,3},{2}}])/Sqrt[6]]] =!= TableauForm,
  True, TestID -> "tableauform-evaluates" ];
```

- [ ] **Step 2: Run red.**
- [ ] **Step 3: Add usage messages** to main file for `Tableau, TensorTableau, Psi, TableauForm, ToTensor, TableauPermute, TableauDot, TensorDot, TensorNorm, TableauNormalization, TableauOrthogonalization` (one short `::usage` each).
- [ ] **Step 4: Create `Kernel/YoungTableaux.wl`** — migrate `YoungTableau/Tableau.wl:31-120` bodies; wrap the six public entry points with a leading `Expand`, e.g.:
```wolfram
ToTensor[expr_] := iToTensor[Expand[expr]];      (* iToTensor has the t_TensorTableau, c_*t, t1_+t2_ rules *)
TableauNormalization[t_] := Module[{e = Expand[t]}, e/TensorNorm[Expand@ToTensor[e]]];
```
(Keep `iToTensor`/`iTableauForm` private with the existing three rules.)
- [ ] **Step 5: Wire `<< YoungTableaux`; run green; commit.**
```bash
git add Kernel/ Tests/YoungTableaux.wlt
git commit -m "feat(tableaux): migrate Young-tableau module with robust linear combinations"
```

---

## Task 8: Protect public API + alias convenience names

**Files:** Modify `Kernel/ClassicalLieAlgebra.wl` (aliases + `Protect` before `EndPackage[]`); Test `Tests/Loading.wlt` (append).

- [ ] **Step 1: Append failing tests** to `Tests/Loading.wlt`
```wolfram
VerificationTest[ MemberQ[Attributes[Generators], Protected], True, TestID -> "generators-protected" ];
VerificationTest[ CartanWeyl[SU[3]], Generators[SU[3], "CartanWeyl"], TestID -> "cartanweyl-alias" ];
VerificationTest[ Chevalley[SU[3]], Generators[SU[3], "Chevalley"], TestID -> "chevalley-alias" ];
VerificationTest[ StringQ[Generators::usage], True, TestID -> "generators-has-usage" ];
```
- [ ] **Step 2: Run red.**
- [ ] **Step 3: Implement** — add `CartanWeyl::usage`, `Chevalley::usage` declarations; in the Private block after all subfiles load, define aliases and protect:
```wolfram
CartanWeyl[g_] := Generators[g, "CartanWeyl"];
Chevalley[g_] := Generators[g, "Chevalley"];
```
Just before `EndPackage[]` (outside `` `Private` ``):
```wolfram
Protect[Evaluate[Names["ClassicalLieAlgebra`*"]]];
```
- [ ] **Step 4: Run green; commit.**
```bash
git add Kernel/ Tests/Loading.wlt
git commit -m "feat(api): CartanWeyl/Chevalley aliases and Protect public symbols"
```

---

## Task 9: Headless test runner

**Files:** Create `scripts/runTests.wls`.

- [ ] **Step 1: Create `scripts/runTests.wls`**
```wolfram
#!/usr/bin/env wolframscript
root = DirectoryName[$InputFileName, 2];
PacletDirectoryLoad[root];
Needs["ClassicalLieAlgebra`"];
files = FileNames["*.wlt", FileNameJoin[{root, "Tests"}]];
reports = AssociationMap[TestReport[#] &, files];
KeyValueMap[Print[FileNameTake[#1], ": ",
   If[#2["AllTestsSucceeded"], "PASS", "FAIL (" <> ToString[#2["TestsFailedCount"]] <> ")"]] &, reports];
allOK = AllTrue[Values[reports], #["AllTestsSucceeded"] &];
Print[If[allOK, "ALL PASS", "FAILURES"]];
Exit[If[allOK, 0, 1]]
```
- [ ] **Step 2: Run, verify it passes and exits 0**

Run: `"$WS" -file scripts/runTests.wls; echo "exit=$?"`
Expected: each test file `PASS`, final `ALL PASS`, `exit=0`.
- [ ] **Step 3: Commit.**
```bash
git add scripts/runTests.wls
git commit -m "test: headless runner over Tests/*.wlt with Exit[1] on failure"
```

---

## Task 10: CI workflow

**Files:** Create `.github/workflows/test.yml`.

- [ ] **Step 1: Create `.github/workflows/test.yml`**
```yaml
name: Test
on: [push, pull_request]
jobs:
  test:
    runs-on: ubuntu-latest
    container:
      image: wolframresearch/wolframengine:latest
      options: --user root
    env:
      WOLFRAMSCRIPT_ENTITLEMENTID: ${{ secrets.WOLFRAMSCRIPT_ENTITLEMENTID }}
    steps:
      - uses: actions/checkout@v4
      - name: Run tests
        run: wolframscript -file scripts/runTests.wls
```
- [ ] **Step 2: Verify YAML parses locally** (optional): `"$WS" -code 'Import[".github/workflows/test.yml","String"]; "ok"'` → `ok`.
- [ ] **Step 3: Commit.** (Note in PR body: add `WOLFRAMSCRIPT_ENTITLEMENTID` repo secret so CI can activate the engine.)
```bash
git add .github/workflows/test.yml
git commit -m "ci: run .wlt suite on the Wolfram Engine container"
```

---

## Task 11: README + remove old duplicated files

**Files:** Modify `README.md` (Usage section); Delete `ClassicalLieAlgebra.wl`, `SpecialUnitary/`, `SpecialOrthogonal/`, `SymplecticUnitary/`, `YoungTableau/`.

- [ ] **Step 1: Update README "Usage"** — replace the `Import[NotebookDirectory[]<>"ClassicalLieAlgebra.wl"]` block with:
````markdown
## Usage

Load from a local checkout:

```mathematica
PacletDirectoryLoad["/path/to/LieAlgebra"];
Needs["ClassicalLieAlgebra`"];
```

Or install a released build:

```mathematica
PacletInstall["https://github.com/jayren3996/LieAlgebra/releases/download/v1.0.0/ClassicalLieAlgebra-1.0.0.paclet"];
Needs["ClassicalLieAlgebra`"];
```
````
Add a short "Breaking changes in 1.0" note: `CartanWeyl`/`Chevalley` return associations; per-element `…H/E/F` accessors, `StandardChevalley`, and `StandardBasis` are removed (`StandardBasis`→`BasisTransform`).

- [ ] **Step 2: Delete the superseded files**
```bash
git rm ClassicalLieAlgebra.wl
git rm -r SpecialUnitary SpecialOrthogonal SymplecticUnitary YoungTableau
```
- [ ] **Step 3: Full suite green after deletion**

Run: `"$WS" -file scripts/runTests.wls; echo "exit=$?"` → `ALL PASS`, `exit=0`.
- [ ] **Step 4: Commit.**
```bash
git add README.md
git commit -m "docs: paclet install/usage; remove superseded duplicated sources"
```

---

## Self-Review

**Spec coverage:** §2 structure → Tasks 1,9,10,11. §3.1 algebra objects → Task 2. §3.2 root-system layer → Task 3. §3.3 unified `Generators`/`BasisTransform` + realization → Tasks 4,5,6. §3.4 tableaux + robustness, drop `ColumnCanonicalize` → Task 7. §3.5 error handling/Protect/SyntaxInformation → Tasks 2,4,8. §4 bug fixes (`` `Private` `` backtick, `Block`→`Module`, scalar-over-sum) → Tasks 4–7. §5 tests → every task + Task 9. §6 CI → Task 10. §7 docs/migration → Task 11. No spec section is unaddressed.

**Placeholder scan:** Migration tasks (4–7) reference exact source line ranges in the current `ClassicalLieAlgebra.wl` plus explicit transformation rules and full test code; the only "fill-in" is mechanically copying verified math bodies, which is intentional (don't re-type 700 verified lines). All new code (paclet, algebra objects, root-system data, dispatch, runner, CI) and all test code is shown in full.

**Type consistency:** `toAssoc` returns `<|"Cartan","Raising","Lowering"|>` in Tasks 4/5/6 (same keys as the `cheval-assoc-keys` test). `Generators[g,scheme,opts]` signature and `Options[Generators]={"Realization"->"Diagonal"}` are declared once (Task 4) and reused. `BasisTransform` used consistently (Tasks 5,6,11). `soMatDim`/`spMatDim` map canonical type→matrix dimension consistently with the defining-dim table in §3.1.

**Open risks flagged for execution:** (1) symbol-name collisions — first action when implementing Task 3 is `"$WS" -code 'Names["System`"<>#]&/@{"Rank","CartanMatrix","SimpleRoots","PositiveRoots","FundamentalWeights"}'`; if any is a `System`` symbol, rename (e.g. `Rank`→`LieAlgebraRank`) in the main file and tests before proceeding. (2) `Generators` is one symbol with downvalues added across four files — the least-specific bad-scheme/bad-input fallthroughs must be defined in the main file *after* all `<<` loads (noted in Task 4).
