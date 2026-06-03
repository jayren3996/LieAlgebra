# Demos

Runnable example scripts for **ClassicalLieAlgebra**. Each is a self-contained
`.wls` script that loads the paclet from a checkout, walks through one theme
with printed output, and ends with a block of self-checks (`[ok]` / `[FAIL]`) so
the demos double as a smoke test of the public API.

| Script | What it shows |
| :--- | :--- |
| [`01-getting-started.wls`](https://github.com/jayren3996/LieAlgebra/blob/master/demos/01-getting-started.wls) | Constructing `SU`/`SO`/`Sp` and the canonical `LieAlgebra[...]` form; root-system data (`Rank`, `CartanMatrix`, `SimpleRoots`, `PositiveRoots`, `FundamentalWeights`); the standard, Cartan–Weyl, and Chevalley generator bases, with a `[E,F]=H` check. |
| [`02-representations.wls`](https://github.com/jayren3996/LieAlgebra/blob/master/demos/02-representations.wls) | `Irrep[g, w]` end to end: dimension, weight system, Casimir eigenvalue, and explicit `RepresentationMatrices`; a dimension/Casimir table across su, so, and sp; the `so(5)` spinor that no tensor power of the vector reaches. |
| [`03-young-tableaux.wls`](https://github.com/jayren3996/LieAlgebra/blob/master/demos/03-young-tableaux.wls) | The Young-tableau toolkit for many-body wavefunctions: the symmetrizer, tensor tableaux and their `Psi` expansion, inner products, normalization, and orthogonalizing a degenerate weight space — the SU(3) computation from the [guided tour](walkthrough.md), made executable. |
| [`04-physics-applications.wls`](https://github.com/jayren3996/LieAlgebra/blob/master/demos/04-physics-applications.wls) | SU(3) flavour and the eightfold way (quarks, the meson/baryon octet, the decuplet, and the dimension counts behind `3⊗3̄` and `3⊗3⊗3`); SU(2) spin with `dim = 2j+1` and `Casimir = 2j(j+1)`; and a numeric check that a built representation obeys the algebra. |

## Running them

The demos print to the console and need no front end, so a command-line kernel is
enough:

```sh
wolframscript -file demos/01-getting-started.wls
```

On macOS the system `wolframscript` in `/usr/local/bin` can be misconfigured; if so,
use the binary bundled with the application:

```sh
/Applications/Wolfram.app/Contents/MacOS/wolframscript -file demos/01-getting-started.wls
```

Each script loads the package from the repository via `PacletDirectoryLoad`, so
you can run them straight from a checkout without installing the paclet first. Run
them from anywhere — the path is resolved relative to the script.

A clean run ends with every self-check reading `[ok]`.
