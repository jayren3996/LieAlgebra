# Representation Engine (Phase 2) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add `Irrep[g, λ]` and accessors (`RepresentationDimension`, `WeightSystem`, `CasimirEigenvalue`, `RepresentationMatrices`, `HighestWeight`) giving exact dimension, weight multiplicities, Casimir, and explicit Chevalley-generator matrices in any irrep of the classical algebras A/B/C/D — including so/sp spinors.

**Architecture:** Two layers on the Phase 1 paclet. **Layer 1 (`Kernel/Weights.wl`)**: closed-form/recursive combinatorics (Weyl dimension, Freudenthal multiplicities, Weyl orbits, Casimir) in Euclidean weight coordinates. **Layer 2 (`Kernel/Representations.wl`)**: the abstract highest-weight/Shapovalov construction — build the module as lowering-words `f_{i_k}…f_{i_1} v_λ`, cap each weight space by the Layer-1 multiplicity, resolve dependence via the rank of the Shapovalov Gram matrix (computed from the Chevalley relations), then read off `H_i` (diagonal), `F_i`, and `E_i = F_i^†`.

**Tech Stack:** Wolfram Language paclet (extends Phase 1: `LieAlgebra`, `Rank`, `CartanMatrix`, `SimpleRoots`, `PositiveRoots`, `FundamentalWeights`, `Generators[g,"Chevalley"]`), `VerificationTest`/`.wlt`.

---

## Conventions for the implementer

- **Wolfram kernel (system one is broken):** `WS=/Applications/Wolfram.app/Contents/MacOS/wolframscript`. Run one suite: `"$WS" -code 'PacletDirectoryLoad["'"$PWD"'"]; r=TestReport["Tests/Weights.wlt"]; Print[r["AllTestsSucceeded"]," ",r["TestsFailedCount"]]'`. Run all: `"$WS" -file scripts/runTests.wls`.
- Each `.wlt` starts with only `Needs["ClassicalLieAlgebra`"];` (the caller does `PacletDirectoryLoad`).
- **Phase 1 idiom (reuse it):** definitions whose LHS contains a literal `LieAlgebra[type, r_]` pattern get corrupted by the validation downvalues — dispatch instead with a guard `f[g_LieAlgebra /; MatchQ[g[[1]], "A"|"B"|"C"|"D"]] := …` and read `g[[1]]` (type), `g[[2]]` (rank), as `Kernel/SpecialOrthogonal.wl` does. Use `Module` (not `Block`) with all locals declared. New files: `BeginPackage["ClassicalLieAlgebra`<Sub>`"]; Needs["ClassicalLieAlgebra`"]; Needs["ClassicalLieAlgebra`Common`"]; Begin["`Private`"]; … End[]; EndPackage[];` and `Get` them from the main file before the final `Protect`.
- TDD: write the `.wlt` test first, run red, implement, run green, commit. Conventional commits ending with `Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>`.
- **Weights internally in Euclidean coords** (use Phase 1's `SimpleRoots`/`PositiveRoots`/`FundamentalWeights`, inner product = `Dot`); convert to/from Dynkin labels only at the public boundary.

---

## File structure

| File | Responsibility |
| :--- | :--- |
| `Kernel/Weights.wl` | Dynkin↔Euclidean, `ρ`, Weyl reflections & orbits, Weyl dimension, Freudenthal multiplicities, Casimir — all on `(g, λ)` internally |
| `Kernel/Representations.wl` | `Irrep[g,λ]` object + validation/display; the Shapovalov engine; the five public accessors |
| `Kernel/ClassicalLieAlgebra.wl` | add `::usage`/messages for new public symbols; `Get` the two new files; extend `Protect` |
| `Tests/Weights.wlt`, `Tests/Representations.wlt` | suites |

Public symbols added: `Irrep`, `HighestWeight`, `RepresentationDimension`, `WeightSystem`, `CasimirEigenvalue`, `RepresentationMatrices`.

---

## Task 1: Weights.wl scaffold + Dynkin↔Euclidean + Weyl vector

**Files:** Create `Kernel/Weights.wl`; Modify main file (`Get`); Test `Tests/Weights.wlt`.

- [ ] **Step 1: collision check (first action).** Run `"$WS" -code 'Print[{#,Names["System`"<>#]}&/@{"Irrep","HighestWeight","RepresentationDimension","WeightSystem","CasimirEigenvalue","RepresentationMatrices"}]'`. Any non-empty `System`` list → rename that symbol (e.g. `WeightSystem`→`WeightMultiplicities`) consistently in plan + tests + code, and report it.

- [ ] **Step 2: write failing test** `Tests/Weights.wlt`:
```wolfram
Needs["ClassicalLieAlgebra`"];
(* internal helpers are tested via their public effects later; here test the boundary maps
   through the exposed RepresentationDimension once it exists. For now test Weyl vector ρ
   indirectly: ρ in Dynkin labels is all 1s, i.e. toDynkin[g, weylVector[g]] == {1,..,1}. *)
VerificationTest[ ClassicalLieAlgebra`Weights`Private`toDynkin[SU[3], ClassicalLieAlgebra`Weights`Private`weylVector[SU[3]]], {1,1}, TestID->"rho-dynkin-A2" ];
VerificationTest[ ClassicalLieAlgebra`Weights`Private`toEuclidean[SU[3], {1,0}], ClassicalLieAlgebra`Weights`Private`toEuclidean[SU[3], {1,0}], TestID->"toEuclidean-defined" ];
VerificationTest[ ClassicalLieAlgebra`Weights`Private`toDynkin[SU[3], ClassicalLieAlgebra`Weights`Private`toEuclidean[SU[3], {2,1}]], {2,1}, TestID->"dynkin-euclid-roundtrip" ];
```

- [ ] **Step 3: run red.**

- [ ] **Step 4: create `Kernel/Weights.wl`:**
```wolfram
BeginPackage["ClassicalLieAlgebra`Weights`"];
Needs["ClassicalLieAlgebra`"];
Needs["ClassicalLieAlgebra`Common`"];
Begin["`Private`"];

(* λ (Dynkin labels) -> Euclidean: sum a_i * omega_i *)
toEuclidean[g_, a_List] := a . FundamentalWeights[g];

(* Euclidean μ -> Dynkin labels: a_i = 2 (μ,α_i)/(α_i,α_i) *)
toDynkin[g_, mu_] := With[{sr = SimpleRoots[g]},
   Table[2 (mu . sr[[i]])/(sr[[i]] . sr[[i]]), {i, Length[sr]}]];

weylVector[g_] := Total[FundamentalWeights[g]];   (* ρ *)

End[];
EndPackage[];
```

- [ ] **Step 5: wire `Get["ClassicalLieAlgebra`Weights`"];` into the main file's Private block (before the Generators fallthrough / Protect). Run green; commit.**
```bash
git add Kernel/ Tests/Weights.wlt
git commit -m "feat(weights): Weights.wl scaffold with Dynkin<->Euclidean and Weyl vector"
```

---

## Task 2: Irrep object + HighestWeight + validation

**Files:** Create `Kernel/Representations.wl`; Modify main (usages/messages, `Get`); Test `Tests/Representations.wlt`.

- [ ] **Step 1: failing test** `Tests/Representations.wlt`:
```wolfram
Needs["ClassicalLieAlgebra`"];
VerificationTest[ HighestWeight[Irrep[SU[3], {1, 1}]], {1, 1}, TestID->"hw" ];
VerificationTest[ Head[Irrep[SU[3], {1, 1}]], Irrep, TestID->"irrep-inert" ];
VerificationTest[ Irrep[SU[3], {1, 1, 0}], $Failed, {Irrep::badweight}, TestID->"wrong-length" ];
VerificationTest[ Irrep[SU[3], {1, -1}], $Failed, {Irrep::badweight}, TestID->"negative-label" ];
```

- [ ] **Step 2: run red.**

- [ ] **Step 3: add to main file declarations:**
```wolfram
Irrep::usage = "Irrep[g, w] represents the irreducible representation of the classical Lie algebra g with highest weight given by the Dynkin labels w (a list of Rank[g] non-negative integers).";
HighestWeight::usage = "HighestWeight[Irrep[g, w]] returns the highest weight w (Dynkin labels).";
Irrep::badweight = "`1` is not a valid highest weight for `2`; expected a list of `3` non-negative integers.";
```

- [ ] **Step 4: create `Kernel/Representations.wl`:**
```wolfram
BeginPackage["ClassicalLieAlgebra`Representations`"];
Needs["ClassicalLieAlgebra`"];
Needs["ClassicalLieAlgebra`Common`"];
Needs["ClassicalLieAlgebra`Weights`"];
Begin["`Private`"];

(* validation: w must be length Rank[g], all non-negative integers. Guard on g_LieAlgebra. *)
validWeightQ[g_, w_] := ListQ[w] && Length[w] === Rank[g] && AllTrue[w, IntegerQ[#] && # >= 0 &];
Irrep[g_LieAlgebra, w_] /; !validWeightQ[g, w] := (Message[Irrep::badweight, w, g, Rank[g]]; $Failed);
(* valid Irrep[g,w] stays inert *)

HighestWeight[Irrep[g_LieAlgebra, w_]] := w;

End[];
EndPackage[];
```
(`SU[3]` sugar evaluates to `LieAlgebra["A",2]`, so `Irrep[SU[3],w]` becomes `Irrep[LieAlgebra["A",2],w]` and matches the `g_LieAlgebra` guard.)

- [ ] **Step 5: wire `Get["ClassicalLieAlgebra`Representations`"];` into main (after Weights). Run green; commit.**
```bash
git add Kernel/ Tests/Representations.wlt
git commit -m "feat(rep): Irrep object with validation and HighestWeight"
```

---

## Task 3: RepresentationDimension (Weyl dimension formula)

**Files:** Modify main (usage), `Kernel/Weights.wl` (internal `weylDim`), `Kernel/Representations.wl` (accessor); Test append to `Tests/Representations.wlt`.

- [ ] **Step 1: failing tests** (append):
```wolfram
VerificationTest[ RepresentationDimension[Irrep[SU[3], {1, 0}]], 3, TestID->"dim-su3-fund" ];
VerificationTest[ RepresentationDimension[Irrep[SU[3], {1, 1}]], 8, TestID->"dim-su3-adjoint" ];
VerificationTest[ RepresentationDimension[Irrep[SO[5], {1, 0}]], 5, TestID->"dim-so5-vector" ];
VerificationTest[ RepresentationDimension[Irrep[SO[5], {0, 1}]], 4, TestID->"dim-so5-spinor" ];
VerificationTest[ RepresentationDimension[Irrep[Sp[4], {1, 0}]], 4, TestID->"dim-sp4-defining" ];
VerificationTest[ RepresentationDimension[Irrep[Sp[4], {0, 1}]], 5, TestID->"dim-sp4-omega2" ];
VerificationTest[ RepresentationDimension[Irrep[SU[4], {0, 0, 0}]], 1, TestID->"dim-trivial" ];
```

- [ ] **Step 2: run red.**

- [ ] **Step 3: add usage** `RepresentationDimension::usage = "RepresentationDimension[Irrep[g, w]] gives the dimension of the irrep (Weyl dimension formula).";`

- [ ] **Step 4: implement.** In `Weights.wl`:
```wolfram
weylDim[g_, lambda_] := Module[{lam = toEuclidean[g, lambda], rho = weylVector[g], pos = PositiveRoots[g]},
   Times @@ Table[((lam + rho) . a)/(rho . a), {a, pos}]];
```
In `Representations.wl` (needs `weylDim` visible: declare it in `ClassicalLieAlgebra`Weights`` context — i.e. give it a bare declaration above `Begin["`Private`"]` in `Weights.wl`, like `Common` exposes helpers — OR call it via its full private name. Cleanest: declare `weylDim`, `freudenthal`, `weylOrbit`, `casimir` in the `Weights`` context (above its Private) so `Representations.wl` can use them after `Needs["ClassicalLieAlgebra`Weights`"]`):
```wolfram
RepresentationDimension[Irrep[g_LieAlgebra, w_]] := weylDim[g, w];
```

- [ ] **Step 5: run green; commit.**
```bash
git add Kernel/ Tests/Representations.wlt
git commit -m "feat(rep): RepresentationDimension via Weyl dimension formula"
```

---

## Task 4: CasimirEigenvalue (with C_n form halving)

**Files:** Modify main (usage), `Weights.wl` (`casimir`), `Representations.wl` (accessor); Test append.

- [ ] **Step 1: failing tests** (append to `Tests/Representations.wlt`):
```wolfram
VerificationTest[ CasimirEigenvalue[Irrep[SU[2], {1}]], 3/2, TestID->"cas-su2-fund" ];
VerificationTest[ CasimirEigenvalue[Irrep[SU[3], {1, 1}]], 6, TestID->"cas-su3-adjoint-2hv" ];
VerificationTest[ CasimirEigenvalue[Irrep[SU[3], {0, 0}]], 0, TestID->"cas-trivial" ];
```
(su(2) fundamental → (n²−1)/n = 3/2; su(3) adjoint → 2 h^∨ = 6 in long-root²=2 units.)

- [ ] **Step 2: run red.**

- [ ] **Step 3: add usage** `CasimirEigenvalue::usage = "CasimirEigenvalue[Irrep[g, w]] gives the eigenvalue of the quadratic Casimir on the irrep, normalized so long roots have squared length 2.";`

- [ ] **Step 4: implement** in `Weights.wl`:
```wolfram
(* form normalized to long-root^2 = 2: raw Dot is already correct for A,B,D;
   for C the raw coords give long-root^2 = 4, so halve. *)
casimirForm[g_, u_, v_] := If[g[[1]] === "C", (u . v)/2, u . v];
casimir[g_, lambda_] := Module[{lam = toEuclidean[g, lambda], rho = weylVector[g]},
   casimirForm[g, lam, lam + 2 rho]];
```
In `Representations.wl`: `CasimirEigenvalue[Irrep[g_LieAlgebra, w_]] := casimir[g, w];`

- [ ] **Step 5: run green; commit.**
```bash
git add Kernel/ Tests/Representations.wlt
git commit -m "feat(rep): CasimirEigenvalue with C_n normalization"
```

---

## Task 5: Weyl reflections and orbits

**Files:** `Kernel/Weights.wl` (internal `weylReflect`, `weylOrbit`, `dynkinLabelOf`); Test append to `Tests/Weights.wlt`.

- [ ] **Step 1: failing tests** (append to `Tests/Weights.wlt`):
```wolfram
(* A2: Weyl group S3, the orbit of the fundamental weight omega_1 has 3 elements *)
VerificationTest[ Length[ClassicalLieAlgebra`Weights`Private`weylOrbit[SU[3], ClassicalLieAlgebra`Weights`Private`toEuclidean[SU[3], {1,0}]]], 3, TestID->"orbit-A2-omega1" ];
(* B2: orbit of the vector-rep highest weight (1,0) = {±e1,±e2} has 4 elements *)
VerificationTest[ Length[ClassicalLieAlgebra`Weights`Private`weylOrbit[SO[5], ClassicalLieAlgebra`Weights`Private`toEuclidean[SO[5], {1,0}]]], 4, TestID->"orbit-B2-vector" ];
```

- [ ] **Step 2: run red.**

- [ ] **Step 3: implement** in `Weights.wl`:
```wolfram
dynkinLabelOf[g_, mu_, i_] := With[{a = SimpleRoots[g][[i]]}, 2 (mu . a)/(a . a)];
weylReflect[g_, mu_, i_] := mu - dynkinLabelOf[g, mu, i] SimpleRoots[g][[i]];
(* BFS with descent-only dedup (Snow): from a weight w, reflect at node i only when its
   i-th Dynkin label is positive; collect every distinct image. mu0 should be dominant. *)
weylOrbit[g_, mu0_] := Module[{r = Rank[g], orbit = {mu0}, frontier = {mu0}, next, c, nu},
   While[frontier =!= {},
     next = {};
     Do[ Do[ c = dynkinLabelOf[g, w, i];
             If[c > 0, nu = weylReflect[g, w, i];
                If[! MemberQ[orbit, nu], AppendTo[orbit, nu]; AppendTo[next, nu]]],
          {i, r}], {w, frontier}];
     frontier = next];
   orbit];
```

- [ ] **Step 4: run green; commit.**
```bash
git add Kernel/ Tests/Weights.wlt
git commit -m "feat(weights): Weyl reflections and orbit enumeration"
```

---

## Task 6: WeightSystem (Freudenthal multiplicities)

**Files:** Modify main (usage), `Weights.wl` (`freudenthal`, full weight system), `Representations.wl` (accessor); Test append to `Tests/Representations.wlt`.

- [ ] **Step 1: failing tests** (append):
```wolfram
(* fundamental of su(3): 3 weights, each mult 1 *)
VerificationTest[ Total[Values[WeightSystem[Irrep[SU[3], {1, 0}]]]], 3, TestID->"ws-su3-fund-total" ];
VerificationTest[ Max[Values[WeightSystem[Irrep[SU[3], {1, 0}]]]], 1, TestID->"ws-su3-fund-mult1" ];
(* adjoint of su(3): total dim 8, zero weight {0,0} has multiplicity 2 *)
VerificationTest[ Total[Values[WeightSystem[Irrep[SU[3], {1, 1}]]]], 8, TestID->"ws-su3-adj-total" ];
VerificationTest[ WeightSystem[Irrep[SU[3], {1, 1}]][{0, 0}], 2, TestID->"ws-su3-adj-zeroweight-mult2" ];
(* consistency with the dimension formula, for a B-type rep *)
VerificationTest[ Total[Values[WeightSystem[Irrep[SO[5], {1, 0}]]]], RepresentationDimension[Irrep[SO[5], {1, 0}]], TestID->"ws-so5-sum-eq-dim" ];
```

- [ ] **Step 2: run red.**

- [ ] **Step 3: add usage** `WeightSystem::usage = "WeightSystem[Irrep[g, w]] gives an association <|weight -> multiplicity, ...|> over all weights of the irrep, weights expressed as Dynkin-label lists.";`

- [ ] **Step 4: implement** in `Weights.wl`. Algorithm (Freudenthal, generate-all variant):
  - Work in Euclidean coords. `lam = toEuclidean[g,λ]`, `rho = weylVector[g]`, positive roots `pos`.
  - Generate all weights by BFS subtracting simple roots from `lam`, keeping those whose Dynkin labels keep them inside the weight polytope (a weight `μ` is in the system iff its dominant Weyl-representative `μ⁺` satisfies `λ − μ⁺ ∈ nonneg-integer span of simple roots`; simplest robust generator: BFS from `lam` subtracting each simple root, and only keep `μ` for which Freudenthal yields a positive multiplicity — the recursion assigns 0 outside). Use the layered order (by level = number of simple roots subtracted) so all `μ+kα` (higher) are computed first.
  - Freudenthal per weight `μ` (Euclidean):
```wolfram
(* denom = (lam+rho).(lam+rho) - (mu+rho).(mu+rho)  (>0 for mu != lam) *)
(* accum = Sum over positive roots alpha, Sum over k>=1 while (mu + k alpha) is a known weight:
            mult[mu + k alpha] * ((mu + k alpha) . alpha)   *)
(* mult[mu] = 2 accum / denom ; base mult[lam] = 1 *)
```
  Provide `freudenthalMults[g, λ]` returning `<|euclideanWeight -> mult|>`, plus `weightSystemDynkin[g, λ]` that maps keys through `toDynkin` and returns `<|dynkinWeight -> mult|>`. Include `Assert`/check that `Mod[2 accum, denom] == 0`.
  In `Representations.wl`: `WeightSystem[Irrep[g_LieAlgebra, w_]] := weightSystemDynkin[g, w];`

- [ ] **Step 5: run green; commit.**
```bash
git add Kernel/ Tests/Representations.wlt
git commit -m "feat(rep): WeightSystem via Freudenthal multiplicity recursion"
```

---

## Task 7: Shapovalov form on lowering-words

**Files:** `Kernel/Representations.wl` (internal word algebra + `shapovalov`); Test append to `Tests/Representations.wlt`.

Data model (all internal/private):
- A **word** is a list of simple-root indices `{i1,…,ik}` meaning `f_{i1} f_{i2} … f_{ik} v_λ`; `{}` is `v_λ`.
- The **weight** of a word (in Dynkin labels) is `λ − Σ_j (row j of CartanMatrix)` for each index `j` in the word (subtracting a simple root `α_j` changes Dynkin labels by minus the `j`-th **row** of the Cartan matrix).
- A **vector** is an `Association <|word -> coefficient|>`.
- `applyH[g, λ, i, word]` = `dynkinLabel_i(weight(word)) * word`.
- `applyE[g, λ, i, vec]`: linear; on a single word: `applyE[i, {}] = 0` (`e_i v_λ = 0`); `applyE[i, {j, rest}] = applyF[j, applyE[i, {rest}]] + (i==j) * applyH[i, {rest}]` using `e_i f_j = f_j e_i + δ_ij h_i`. (`applyF[j, vec]` prepends `j` to every word.)
- `shapovalov[g, λ, vec1, vec2]`: bilinear; for single words, `⟨{}, {}⟩ = 1`, and `⟨{i,rest1}, w2⟩ = shapovalov[rest1, applyE[i, w2]]` (using `⟨f_i u, w⟩ = ⟨u, e_i w⟩`); `⟨{}, w2⟩ = coefficient of {} in w2`.

- [ ] **Step 1: failing tests** (append; norms verified from the su(2) spin-1 worked example):
```wolfram
(* su(2), highest weight {2} (spin 1). Words {}, {1}, {1,1}. Shapovalov norms: 1, 2, 4. *)
With[{rp = ClassicalLieAlgebra`Representations`Private`}, Null];
VerificationTest[ ClassicalLieAlgebra`Representations`Private`shapovalov[SU[2], {2}, <|{} -> 1|>, <|{} -> 1|>], 1, TestID->"shap-hw-norm" ];
VerificationTest[ ClassicalLieAlgebra`Representations`Private`shapovalov[SU[2], {2}, <|{1} -> 1|>, <|{1} -> 1|>], 2, TestID->"shap-f-norm" ];
VerificationTest[ ClassicalLieAlgebra`Representations`Private`shapovalov[SU[2], {2}, <|{1, 1} -> 1|>, <|{1, 1} -> 1|>], 4, TestID->"shap-ff-norm" ];
```
(`⟨f v_λ, f v_λ⟩ = ⟨v_λ, e f v_λ⟩ = ⟨v_λ, h v_λ⟩ = 2`; `⟨f² v_λ, f² v_λ⟩ = 4` as in the worked example.)

- [ ] **Step 2: run red.**

- [ ] **Step 3: implement** the word algebra + `shapovalov` in `Representations.wl` per the data model above. Use `CartanMatrix[g]` for the weight bookkeeping and `Rank[g]`. Keep everything exact (rationals/integers).

- [ ] **Step 4: run green; commit.**
```bash
git add Kernel/ Tests/Representations.wlt
git commit -m "feat(rep): Shapovalov contravariant form on lowering-words"
```

---

## Task 8: RepresentationMatrices (build the module; H, F, E)

**Files:** Modify main (usage), `Representations.wl` (engine + accessor); Test append.

Algorithm (uses Task 6 multiplicities + Task 7 Shapovalov):
1. BFS over weights by level from `v_λ` (`{}`). At each step apply each `f_i` to current basis words to produce candidate words at weight `ν`.
2. For each weight `ν`, accumulate candidate words from all sources; using the Shapovalov Gram matrix of the accumulated candidates, pick a maximal subset whose Gram rank equals `mult(ν)` (from Task 6). Orthonormalize that weight space (Gram–Schmidt against the Shapovalov form) → an orthonormal basis expressed as vectors (`<|word->coef|>`).
3. Order all orthonormal basis vectors by weight → indices `1..dim`.
4. `H_i` = `DiagonalMatrix[dynkinLabel_i(weight(b))` for each basis vector `b]`.
5. `F_i[[a,b]]` = Shapovalov coordinate of `f_i·(basis b)` on `basis a` = `⟨basis_a, f_i·basis_b⟩` (in the orthonormal basis, the coordinate is the inner product). `f_i·basis_b` = `applyF[i, basis_b]`.
6. `E_i = ConjugateTranspose[F_i]`.
7. Memoize per `(g, λ)`. Option `"Normalized"->False` returns the unnormalized rational basis (skip the √ normalization; then `E_i` via the diagonal metric, not conjugate-transpose).

- [ ] **Step 1: failing tests** (append; robust, basis-ordering-independent):
```wolfram
bracket[a_, b_] := a.b - b.a;
(* su(2) spin-1 ({2}): 3x3, [E,F]=H, H eigenvalues {2,0,-2} *)
VerificationTest[ Module[{m = RepresentationMatrices[Irrep[SU[2], {2}]]},
   bracket[m["Raising"][[1]], m["Lowering"][[1]]] == m["Cartan"][[1]]], True, TestID->"su2-spin1-EF=H" ];
VerificationTest[ Sort[Diagonal[RepresentationMatrices[Irrep[SU[2], {2}]]["Cartan"][[1]]]], {-2, 0, 2}, TestID->"su2-spin1-H-spectrum" ];
(* su(3) fundamental: 3x3, both [E_i,F_i]=H_i *)
VerificationTest[ Module[{m = RepresentationMatrices[Irrep[SU[3], {1, 0}]]},
   And @@ Table[bracket[m["Raising"][[i]], m["Lowering"][[i]]] == m["Cartan"][[i]], {i, 2}]], True, TestID->"su3-fund-EF=H" ];
VerificationTest[ Length[RepresentationMatrices[Irrep[SU[3], {1, 0}]]["Cartan"][[1]]], 3, TestID->"su3-fund-size3" ];
(* off-diagonal Chevalley relation [E_1,F_2]=0 for su(3) *)
VerificationTest[ Module[{m = RepresentationMatrices[Irrep[SU[3], {1, 0}]]},
   bracket[m["Raising"][[1]], m["Lowering"][[2]]] == 0 IdentityMatrix[3]], True, TestID->"su3-fund-EF-offdiag-0" ];
```

- [ ] **Step 2: run red.**

- [ ] **Step 3: add usage** `RepresentationMatrices::usage = "RepresentationMatrices[Irrep[g, w]] gives <|\"Cartan\"->{H_i}, \"Raising\"->{E_i}, \"Lowering\"->{F_i}|>, the Chevalley generators of g (one per simple root) as matrices in the irrep. Option \"Normalized\"->True (default) returns an orthonormal basis with E_i = ConjugateTranspose[F_i]; False returns the exact rational unnormalized basis.";` Add `Options[RepresentationMatrices] = {"Normalized" -> True};`

- [ ] **Step 4: implement** the engine + accessor `RepresentationMatrices[Irrep[g_LieAlgebra, w_], OptionsPattern[]] := …` per the algorithm. Build matrices in the weight-ordered orthonormal basis.

- [ ] **Step 5: run green; commit.**
```bash
git add Kernel/ Tests/Representations.wlt
git commit -m "feat(rep): RepresentationMatrices via the Shapovalov highest-weight construction"
```

---

## Task 9: Validation — degenerate weights, spinors, sp, and the Phase-1 tie-in

**Files:** Test append to `Tests/Representations.wlt` (no new code unless a test surfaces a bug; if it does, fix the engine).

- [ ] **Step 1: add tests** exercising the hard cases:
```wolfram
bracket[a_, b_] := a.b - b.a;
chevChecks[g_, w_] := Module[{m = RepresentationMatrices[Irrep[g, w]], A = CartanMatrix[g], r},
   r = Rank[g];
   (* [E_i,F_j] = delta_ij H_i  AND  [H_i,E_j] = A[j,i] E_j *)
   (And @@ Flatten[Table[ bracket[m["Raising"][[i]], m["Lowering"][[j]]] ==
        If[i == j, m["Cartan"][[i]], 0 m["Cartan"][[1]]], {i, r}, {j, r}]]) &&
   (And @@ Flatten[Table[ bracket[m["Cartan"][[i]], m["Raising"][[j]]] == A[[j, i]] m["Raising"][[j]], {i, r}, {j, r}]])];

(* su(3) adjoint: 8-dim with a 2-dim zero-weight space (the degenerate case) *)
VerificationTest[ Length[RepresentationMatrices[Irrep[SU[3], {1, 1}]]["Cartan"][[1]]], 8, TestID->"su3-adjoint-dim8" ];
VerificationTest[ chevChecks[SU[3], {1, 1}], True, TestID->"su3-adjoint-chevalley-relations" ];
(* so(5) SPINOR — the case a tableau engine cannot do *)
VerificationTest[ Length[RepresentationMatrices[Irrep[SO[5], {0, 1}]]["Cartan"][[1]]], 4, TestID->"so5-spinor-dim4" ];
VerificationTest[ chevChecks[SO[5], {0, 1}], True, TestID->"so5-spinor-chevalley-relations" ];
(* sp(4) defining *)
VerificationTest[ chevChecks[Sp[4], {1, 0}], True, TestID->"sp4-defining-chevalley-relations" ];
(* Phase-1 tie-in: fundamental irrep H-spectra match the Chevalley defining rep's H diagonals *)
VerificationTest[
   Sort /@ (Diagonal /@ RepresentationMatrices[Irrep[SU[3], {1, 0}]]["Cartan"]) ===
   Sort /@ (Diagonal /@ Generators[SU[3], "Chevalley"]["Cartan"]),
   True, TestID->"su3-fund-matches-phase1-Hspectra" ];
```

- [ ] **Step 2: run.** If any fail, debug the engine (Task 7/8) until green — do NOT weaken the tests. Likely suspects: weight bookkeeping sign (Cartan row vs column), Gram-rank capping in degenerate spaces, C-type normalization. Add focused diagnostics with the `WS` kernel.

- [ ] **Step 3: commit.**
```bash
git add Tests/Representations.wlt Kernel/
git commit -m "test(rep): degenerate-weight, spinor, sp, and Phase-1 consistency checks"
```

---

## Task 10: Protect, docs, README

**Files:** Modify main (`Protect` already auto-covers new symbols via `Names["ClassicalLieAlgebra`*"]` — verify), `README.md`.

- [ ] **Step 1: verify** in a fresh kernel that all six new symbols are `Protected` (the existing `Protect[Evaluate[Names["ClassicalLieAlgebra`*"]]]` should already cover them since they're declared in the main context): `"$WS" -code 'PacletDirectoryLoad["'"$PWD"'"]; Needs["ClassicalLieAlgebra`"]; Print[AllTrue[{Irrep,RepresentationDimension,WeightSystem,CasimirEigenvalue,RepresentationMatrices,HighestWeight}, MemberQ[Attributes[#],Protected]&]]'` → `True`. If any aren't protected, ensure they have `::usage` in the main file (so they're created in the public context before `Protect`).

- [ ] **Step 2: add a "Representations" section to `README.md`** with a worked example:
````markdown
## Representations

Build an irrep from its highest weight (Dynkin labels):

```mathematica
ir = Irrep[SU[3], {1, 1}];          (* the adjoint / octet *)
RepresentationDimension[ir]          (* 8 *)
CasimirEigenvalue[ir]                (* 6 *)
WeightSystem[ir]                     (* <|{1,1}->1, {0,0}->2, ...|> *)
m = RepresentationMatrices[ir];      (* <|"Cartan"->{H1,H2}, "Raising"->{E1,E2}, "Lowering"->{F1,F2}|> *)
```

Works for all classical types, including so/sp spinor representations, e.g. `Irrep[SO[5], {0, 1}]` (the 4-dimensional spinor).
````

- [ ] **Step 3: run the full suite** `"$WS" -file scripts/runTests.wls` → ALL PASS, exit 0. Commit.
```bash
git add Kernel/ README.md
git commit -m "docs: Representations section; confirm new symbols protected"
```

---

## Self-Review

**Spec coverage:** §2 API (`Irrep`+5 accessors) → Tasks 2,3,4,6,8 (+`HighestWeight` Task 2). §3 conventions (Euclidean/Dynkin, ρ, C_n halving) → Tasks 1,4,6. §4 Layer 1 (dim/Freudenthal/orbits/Casimir) → Tasks 3,4,5,6. §5 Layer 2 (Shapovalov, BFS, H/F/E) → Tasks 7,8. §6 module structure → Tasks 1,2. §7 validation (known dims, [E,F]=H, spinor, Phase-1 tie) → Tasks 3,4,6,8,9. §8 phasing → task order matches. §9 decisions baked in (orthonormal default w/ option; Dynkin-label weights; collision check Task 1; scope guard — no tensor products/branching/GT). No gaps.

**Placeholder scan:** Task 6 and Task 8 describe algorithms with explicit data models/formulas and complete test code rather than full final source (the Freudenthal recursion and the module-building BFS are given step-by-step with the exact recurrence and the determining test values). This is intentional for the two genuinely algorithmic tasks; every other task has complete code. No "TBD/handle edge cases" placeholders.

**Type consistency:** `Irrep[g_LieAlgebra, w_]` dispatch and the `g[[1]]`/`Rank[g]` access are consistent across tasks. `RepresentationMatrices` returns the `<|"Cartan","Raising","Lowering"|>` association (same keys as Phase 1 `Generators[g,"Chevalley"]`) everywhere it's used (Tasks 8, 9, 10). Internal `Weights`` helpers (`toEuclidean`, `toDynkin`, `weylVector`, `weylDim`, `casimir`, `weylOrbit`, `freudenthal…`) are named consistently and declared in the `Weights`` context so `Representations.wl` can call them. `WeightSystem` returns Dynkin-label keys in every reference.

**Risk note for execution:** Tasks 7–8 are the hard core (the abstract Shapovalov engine). Build strictly bottom-up with the given exact test values (su(2) spin-1 norms/matrices are fully worked); when a higher case (su(3) adjoint degenerate space, so(5) spinor) fails, the bug is almost always in (a) the Dynkin-label weight bookkeeping (Cartan **row** when subtracting a simple root), (b) Gram-rank capping in multiplicity>1 spaces, or (c) the `e_i f_j = f_j e_i + δ_ij h_i` reduction. Verify each against the kernel before moving on.
