# Tutorials

Four worked walkthroughs of **ClassicalLieAlgebra**, each mirroring a
self-contained, runnable `.wls` script in the repository's
[`demos/`](https://github.com/jayren3996/LieAlgebra/tree/master/demos) directory.
Read them here, or run the scripts yourself — every one ends with a block of
self-checks (`[ok]` / `[FAIL]`), so the demos double as a smoke test of the public
API and the output you see below is exactly what the package produces.

| Tutorial | Source script | What it shows |
| :--- | :--- | :--- |
| [Getting started](getting-started.md) | [`01-getting-started.wls`](https://github.com/jayren3996/LieAlgebra/blob/master/demos/01-getting-started.wls) | Constructing the algebras, root-system data, and the three generator bases. |
| [Representations](representations.md) | [`02-representations.wls`](https://github.com/jayren3996/LieAlgebra/blob/master/demos/02-representations.wls) | `Irrep[g, w]` end to end: dimension, weights, Casimir, generator matrices — across `su`, `so`, `sp`, and the spinors. |
| [Young tableaux](young-tableaux.md) | [`03-young-tableaux.wls`](https://github.com/jayren3996/LieAlgebra/blob/master/demos/03-young-tableaux.wls) | The many-body wavefunction toolkit: the symmetrizer, inner products, normalization, and orthogonalization. |
| [Physics applications](physics-applications.md) | [`04-physics-applications.wls`](https://github.com/jayren3996/LieAlgebra/blob/master/demos/04-physics-applications.wls) | SU(3) flavour and the eightfold way, SU(2) spin, and a numeric check that a built representation obeys the algebra. |

## Running the scripts

The demos print to the console and do not need the notebook front end:

```sh
wolframscript -file demos/01-getting-started.wls
```

Each script loads the package from the repository via `PacletDirectoryLoad`, so you
can run them straight from a checkout without installing the paclet first. If your
shell cannot find `wolframscript`, use the Wolfram application binary or add it to
your `PATH`. A clean run ends with every self-check reading `[ok]`.
