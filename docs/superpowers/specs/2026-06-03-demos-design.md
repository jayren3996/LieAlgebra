# ClassicalLieAlgebra demos — design

## Goal

Add a set of runnable, self-verifying example scripts that show how to use the
v1.0 paclet, and retire the stale `Demo.nb`.

## Motivation

The only existing demo, `Demo.nb`, predates the v1.0 API redesign and no longer
runs: it references removed symbols (`StandardChevalley`, the old `StandardBasis`,
the `{H, E, F}` list returns) and the since-corrected `sp(2n)` generators. A new
user therefore has no current, executable example to start from. The README
quick-start and `docs/walkthrough.md` are prose with code blocks; neither is a
script you can run end to end, and the walkthrough only covers SU(3).

## Deliverables

A new `demos/` folder of `.wls` scripts. Each is self-contained, runs headless
(console output only, no front end), and uses the repo's existing preamble
pattern — resolve the repo root, `PacletDirectoryLoad[root]`, then
`Needs["ClassicalLieAlgebra`"]` — the same one `scripts/runTests.wls` and
`pics/MakeFigures.wls` already use.

Each script ends with a small block of assertions that print `ok` / `FAIL`, so
the demos double as smoke tests and cannot silently rot against a future API
change the way `Demo.nb` did.

- `demos/01-getting-started.wls` — algebra construction and canonicalization
  (`SU`/`SO`/`Sp` → `LieAlgebra[...]`); root-system data (`Rank`,
  `LieAlgebraDimension`, `CartanMatrix`, `SimpleRoots`, `PositiveRoots`,
  `FundamentalWeights`); the three generator bases (standard, Cartan–Weyl,
  Chevalley) with a `[E, F] = H` bracket sanity check.
- `demos/02-representations.wls` — the `Irrep` engine: dimension (Weyl), weight
  system (Freudenthal), Casimir eigenvalue, and explicit `RepresentationMatrices`;
  a dimension table across the su/so/sp families; the `so(5)` spinor as the case
  that does not live in tensor powers of the defining representation.
- `demos/03-young-tableaux.wls` — the SU(3) (1,1) tableau computation from the
  walkthrough, made executable: the Young symmetrizer (`TableauPermute`), inner
  products (`TableauDot`), normalization, and orthogonalization of the degenerate
  (0,0) weight space.
- `demos/04-physics-applications.wls` — SU(3) flavour / the eightfold way (octet
  and decuplet from Dynkin labels and weights); SU(2) spin and the Casimir
  relation `C = 2 j (j+1)`; and a numeric check that the built octet matrices
  satisfy the su(3) Chevalley relations.
- `demos/README.md` — index of the four demos and how to run them (the macOS
  app-bundled `wolframscript` binary, plus plain `wolframscript -file` elsewhere).

## Cleanups (from the review)

- Remove `Demo.nb`; repoint its references in `README.md` and
  `docs/walkthrough.md` to `demos/`.
- Add a `.gitignore` (`.DS_Store`, `*.paclet`, common junk); untrack the
  committed `.DS_Store`.
- Widen the `PacletInfo` `"Description"` to mention the representation engine,
  which it currently omits.

## Out of scope

- Notebook (`.nb`) demos. Chosen format is plain-text `.wls`.
- Deduping the repeated "Regenerating the figures" section.
- Any change to package source under `Kernel/`.

## Verification

- Run each demo with the Wolfram binary
  (`/Applications/Wolfram.app/Contents/MacOS/wolframscript -file demos/NN-*.wls`);
  confirm a clean exit, sensible output, and every assertion `ok`.
- Re-run the `Tests/*.wlt` suite via `scripts/runTests.wls`; confirm no
  regressions from the cleanups.

## Risks

- Demo code drifting from the real API. Mitigated by running every demo and by
  the inline assertions.
- Headless rendering: `MatrixForm` prints as a text grid in a script rather than
  typeset output. The demos label what they print and do not depend on a front
  end.
