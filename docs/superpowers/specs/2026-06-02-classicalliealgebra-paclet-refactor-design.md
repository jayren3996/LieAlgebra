# ClassicalLieAlgebra — Phase 1 design: paclet refactor + API redesign

**Date:** 2026-06-02
**Status:** Draft for review
**Scope:** Phase 1 of a two-phase program. This phase delivers a properly structured Wolfram paclet with a redesigned, robust API for the package's **current** capabilities (generators, bases, Young-tableau machinery) plus a new root-system data layer. Phase 2 (a high-level representation engine, `Irrep` & friends) is out of scope here and gets its own design cycle; this phase is built so Phase 2 can stand on it.

## 1. Goals

1. **Single source of truth.** Eliminate the ~50% duplication: the monolithic `ClassicalLieAlgebra.wl` and the four subfolder packages (`SpecialUnitary/`, `SpecialOrthogonal/`, `SymplecticUnitary/`, `YoungTableau/`) collapse into one paclet, reconciling the drift in favor of the more advanced `Tableau` versions.
2. **Conventional paclet structure** following the official WolframResearch reference layout.
3. **Redesigned API** (the "aggressive" option): canonical `LieAlgebra[type, rank]` objects with `SU/SO/Sp` sugar, a root-system data layer, and one association-returning basis interface.
4. **Robustness:** usage messages, argument validation, error messages, `Protect`; fix the known bugs (context-leak, `Block`/scoping, scalar-over-sum non-distribution, dangling `ColumnCanonicalize`).
5. **Automated tests + CI** on the free Wolfram Engine.

Non-goals (Phase 2 or later): the `Irrep`/weight-system engine; Wolfram documentation notebooks; Paclet Repository submission.

## 2. Target structure

```
LieAlgebra/                       repo root = paclet root (PacletInfo.wl here; Name = "ClassicalLieAlgebra")
├── PacletInfo.wl
├── Kernel/
│   ├── ClassicalLieAlgebra.wl    main: BeginPackage, all ::usage + ::messages, Get subfiles, Protect, EndPackage
│   ├── Common.wl                 sub-context ClassicalLieAlgebra`Common`: shared private helpers
│   ├── Algebras.wl               canonical LieAlgebra[], SU/SO/Sp sugar, validation, display, root-system data
│   ├── SpecialUnitary.wl         su(n) matrix constructors (standard / Cartan–Weyl / Chevalley)
│   ├── SpecialOrthogonal.wl      so(n) matrix constructors + realizations + BasisTransform
│   ├── Symplectic.wl             sp(2n) matrix constructors + realizations + BasisTransform
│   └── YoungTableaux.wl          Tableau / TensorTableau / Psi / …
├── Tests/
│   ├── Algebras.wlt
│   ├── SpecialUnitary.wlt
│   ├── SpecialOrthogonal.wlt
│   ├── Symplectic.wlt
│   └── YoungTableaux.wlt
├── scripts/runTests.wls          runs all .wlt, Exit[1] on failure
├── .github/workflows/test.yml    Wolfram Engine container CI
├── pics/ …                       existing figures + MakeFigures.wls (unchanged)
├── README.md                     install/usage updated to paclet flow
├── LICENSE
└── Demo.nb                        kept as an example (not loaded by the paclet)
```

The duplicated subfolder packages and their `Test.nb` notebooks are **deleted** (content migrated into `Kernel/` and `Tests/`).

### PacletInfo.wl

```wolfram
PacletObject[<|
  "Name"            -> "ClassicalLieAlgebra",
  "Version"         -> "1.0.0",
  "WolframVersion"  -> "13.0+",
  "Description"     -> "Generators, bases, and Young-tableau machinery for the classical Lie algebras.",
  "Creator"         -> "Jie Ren",
  "License"         -> "MIT",
  "SourceControlURL"-> "https://github.com/jayren3996/LieAlgebra",
  "Extensions"      -> {
    {"Kernel", "Root" -> "Kernel", "Context" -> "ClassicalLieAlgebra`"}
  }
|>]
```

No `Kernel/init.m` (discouraged). `Needs["ClassicalLieAlgebra`"]` auto-resolves to `Kernel/ClassicalLieAlgebra.wl`, which declares all public `::usage`/messages and then `Get`s the subfiles.

### Context discipline (three tiers)

- **Public API:** `::usage` declared in `Kernel/ClassicalLieAlgebra.wl` (context `ClassicalLieAlgebra`).
- **Cross-file private helpers:** declared bare in `ClassicalLieAlgebra`Common`` (above its `` `Private` ``), defined inside it; feature files `Needs` it.
- **File-local helpers:** undeclared symbols inside each feature file's `` `Private` `` block.

Every file uses `` Begin["`Private`"] `` (with the leading backtick — fixes the current bug) and `Module` (not `Block`) for locals.

## 3. Public API

### 3.1 Algebra objects

Canonical head: **`LieAlgebra[type, rank]`**, `type ∈ {"A","B","C","D"}`, `rank` a positive integer. Sugar constructors normalize to it:

| Sugar | Canonical | Constraint |
| :--- | :--- | :--- |
| `SU[n]` | `LieAlgebra["A", n-1]` | integer `n ≥ 2` |
| `SO[n]`, odd `n=2r+1` | `LieAlgebra["B", r]` | integer `n ≥ 3` |
| `SO[n]`, even `n=2r` | `LieAlgebra["D", r]` | integer `n ≥ 4` |
| `Sp[n]`, even `n=2r` | `LieAlgebra["C", r]` | even integer `n ≥ 2` |

- Defining-representation matrix dimension: `A_{n}→n+1`, `B_r→2r+1`, `C_r→2r`, `D_r→2r`.
- Sugar heads **evaluate to** the canonical form (so `SU[3]` returns `LieAlgebra["A", 2]`); all downstream functions dispatch on the canonical form, and accept either form as input.
- Low-rank coincidences (e.g. `D_2 ≅ A_1⊕A_1`, `C_1 ≅ A_1`, `B_1 ≅ A_1`) are **allowed** for matrix construction and documented; strict simplicity enforcement is a Phase-2 concern.
- **Display:** a `MakeBoxes` rule pretty-prints `LieAlgebra["A", 2]` as `A₂` (type letter with subscript rank). Sugar forms are not retained for display.
- **Validation:** invalid type, non-integer/out-of-range rank, or odd `Sp` argument → a `Message` + `$Failed` (see §3.5).

### 3.2 Root-system data layer (new — foundation for Phase 2)

All take a canonical algebra `g` (or sugar):

| Function | Returns |
| :--- | :--- |
| `Rank[g]` | the rank (integer) |
| `LieAlgebraDimension[g]` | dimension of the algebra (integer) |
| `CartanMatrix[g]` | the `rank × rank` Cartan matrix |
| `SimpleRoots[g]` | simple roots (list of vectors) |
| `PositiveRoots[g]` | positive roots (list of vectors) |
| `FundamentalWeights[g]` | fundamental weights (list of vectors) |

These are computed from `type`/`rank` (standard Cartan data), independent of the matrix realizations. Roots and weights are expressed in one documented, consistent convention (Bourbaki / Humphreys ordering of simple roots; weights in the same basis), so that the Cartan matrix, roots, and weights all agree. (Symbol names are provisional pending a `System`` collision check — see §8.)

### 3.3 Unified basis interface

One association-returning function replaces `CartanWeyl` / `Chevalley` / `StandardChevalley` and their nine per-element accessors:

```wolfram
Generators[g]                         (* "Standard": flat list of defining-rep generators *)
Generators[g, "CartanWeyl"]           (* -> <|"Cartan"->{Hᵢ}, "Raising"->{E_α}, "Lowering"->{F_α}|> *)
Generators[g, "Chevalley"]            (* -> same association shape *)
```

- **`"Realization"` option** (default `"Diagonal"`) selects the matrix realization for the **`"Chevalley"`** scheme on `so`/`sp`:
  - `"Diagonal"` — Cartan generators diagonal (current `Chevalley[SO]/[Sp]`; the rep-useful basis).
  - `"Antisymmetric"` — aligned with the standard antisymmetric generators (current `StandardChevalley`).
  - For `su(n)` the Cartan is already diagonal; the option is a no-op.
  - The `"CartanWeyl"` scheme is returned in the antisymmetric realization (matching the current code); the diagonal version, if ever needed, is just a conjugation by `BasisTransform[g]`. Phase 1 ships the realizations the current code already provides.
- **`BasisTransform[g]`** returns the change-of-basis matrix between the antisymmetric and diagonal realizations (current `StandardBasis`; defined for `so`/`sp`, identity for `su`).
- The matrices returned are **identical** to those the current code produces; only the packaging (association, option) changes. This keeps the math verified-equivalent.
- **Optional convenience aliases** (recommended, for `?`-discoverability): `CartanWeyl[g] := Generators[g, "CartanWeyl"]`, `Chevalley[g] := Generators[g, "Chevalley"]`. (Decision in §8.)

The niche multi-argument individual-generator forms (`Generators[SU[n], i, {i,j}]`, `Generators[SO[n], {i,j}]`, `Generators[Sp[n], i, {i,j}]`) are **dropped**; callers index the returned list/association.

### 3.4 Young-tableau module (kept, cleaned)

Public symbols retained: `Tableau`, `TensorTableau`, `Psi`, `TableauForm`, `ToTensor`, `TableauPermute`, `TableauDot`, `TensorDot`, `TensorNorm`, `TableauNormalization`, `TableauOrthogonalization`.

Changes:
- **Robust linear combinations.** `ToTensor`, `TableauForm`, `TensorNorm`, `TableauDot`, `TableauNormalization`, `TableauOrthogonalization` are made to handle a scalar times a sum of tableaux (e.g. `(T[a]+T[b])/Sqrt[6]`) by normalizing/expanding internally, fixing the current non-distribution that silently leaves expressions unevaluated.
- **`ColumnCanonicalize`** (declared public in the `Tableau` submodule but never defined) is **removed** unless we choose to implement it (Decision §8).
- Names standardized; behavior otherwise preserved and pinned by tests.

### 3.5 Error handling conventions

- Validate at **public entry points only**; private helpers assume valid input.
- Cheap structural checks via argument patterns (`_Integer?Positive`, `{i_Integer, j_Integer}`, …) so wrong-shaped calls stay unevaluated.
- Semantic errors (odd `Sp` rank, unknown type, bad option value) → `Message[f::tag, …]` then return `$Failed`. Named templates declared with the usage messages, e.g.:
  ```wolfram
  LieAlgebra::badtype = "`1` is not a valid Cartan type; use \"A\", \"B\", \"C\", or \"D\".";
  Sp::evenrank       = "Sp[`1`] requires an even positive integer.";
  Generators::badscheme = "`1` is not a valid basis scheme; use \"Standard\", \"CartanWeyl\", or \"Chevalley\".";
  ```
- Public symbols `Protect`ed at load end via `Protect[Evaluate[Names["ClassicalLieAlgebra`*"]]]` (replaces the hand-maintained `protectlist`). `SyntaxInformation` set for the main entry points.

## 4. Bug fixes folded in

- `` Begin["Private`"] `` → `` Begin["`Private`"] `` (4 files).
- `Block` → `Module` throughout; **every** local (including loop counters) declared — fixes the global `i`/`j` leaks in `SUCWH`, `SOBasis`, and similar.
- Tableau scalar-over-sum non-distribution (§3.4).
- Dangling `ColumnCanonicalize` (§3.4).

## 5. Tests

- `Tests/*.wlt` using `VerificationTest`, one file per `Kernel/` module.
- Coverage: port the existing `YoungTableau/Test.wls` cases; encode the README/Demo SU(3) (1,1)-rep computations; verify each algebra's generators close under the Lie bracket and that `BasisTransform` conjugates the antisymmetric realization into the diagonal one; add a **regression for `StandardBasis[USp[2n]]`** (the recent fix); validation tests asserting `$Failed` + the right message on bad input.
- `scripts/runTests.wls`: load the paclet via `PacletDirectoryLoad`, run `TestReport` over `Tests/*.wlt`, print a summary, and `Exit[1]` if any test fails (so CI goes red).

## 6. CI

`.github/workflows/test.yml`: run on the `wolframresearch/wolframengine:latest` container, activate via a `WOLFRAMSCRIPT_ENTITLEMENTID` repo secret (fallback: `WOLFRAM_ID`/`WOLFRAM_PASSWORD`), execute `scripts/runTests.wls`. **User action required:** add the secret to the GitHub repo for CI to run; the workflow is committed regardless.

## 7. Documentation & migration

- `::usage` on every public symbol.
- README "Usage" section updated from `Import[NotebookDirectory[]<>…]` to the paclet flow: dev = `PacletDirectoryLoad["<repo>"]; Needs["ClassicalLieAlgebra`"]`; release = `CreatePacletArchive` → GitHub release → `PacletInstall[url]`.
- **Breaking changes** (documented in README + a CHANGELOG note): `CartanWeyl`/`Chevalley` now return associations (or are aliases of `Generators[…]`); the per-element `…H/E/F` accessors, `StandardChevalley`, and `StandardBasis` are removed/renamed (`StandardBasis`→`BasisTransform`); the niche individual-generator forms are dropped.

### Implementation phasing (within Phase 1)

1. Scaffold the paclet (`PacletInfo.wl`, `Kernel/ClassicalLieAlgebra.wl`, `Common.wl`); get `Needs` loading an empty shell.
2. `Algebras.wl`: canonical objects, sugar, validation, display, root-system data.
3. Migrate `su`/`so`/`sp` matrix constructors into their files (`Block`→`Module`, backtick fix), wired to `Generators[g, scheme]` + `BasisTransform`.
4. `YoungTableaux.wl`: migrate, de-dup, robust linear combinations.
5. Usage messages, validation/messages, `Protect`, `SyntaxInformation`.
6. `Tests/*.wlt` + `scripts/runTests.wls`; get them green locally.
7. CI workflow; README/install updates; delete the old duplicated files.

## 8. Decisions for review

1. **Default `"Realization"`** for `Generators[g, "Chevalley"]` on `so`/`sp`: proposed **`"Diagonal"`** (the rep-useful basis the README uses).
2. **Convenience aliases** `CartanWeyl[g]` / `Chevalley[g]` alongside `Generators[g, "…"]`: proposed **keep** (discoverability) — or drop for a single interface.
3. **Symbol names** `Rank` / `LieAlgebraDimension` / `CartanMatrix` / `SimpleRoots` / `PositiveRoots` / `FundamentalWeights`: verify against `System`` in a clean kernel; rename any that collide (e.g. fall back to `LieAlgebraRank`). Proposed names as listed.
4. **`ColumnCanonicalize`**: proposed **remove** (never implemented). Alternative: implement it as column-wise canonicalization.
5. **Individual-generator accessor forms**: proposed **drop**.
6. **Paclet name**: keep `ClassicalLieAlgebra` (context `ClassicalLieAlgebra`) even though the repo folder is `LieAlgebra`. Proposed keep.

## 9. Phase 2 outline (out of scope here)

A high-level representation engine built on the root-system layer: `Irrep[g, highestWeightDynkin]` returning a representation object exposing dimension, weight system, basis states, and the generator matrices in the irrep. SU(n) first (via the tableau machinery, generalizing the README procedure), then B/C/D (traceless / symplectic-traceless constructions). Designed and specced separately after Phase 1 lands.
```
