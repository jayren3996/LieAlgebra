<div align="center">

# ClassicalLieAlgebra

**Exact generators, bases, Young tableaux, and irreducible representations of the classical Lie algebras — for the Wolfram Language.**

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Wolfram Language](https://img.shields.io/badge/Wolfram%20Language-13.0%2B-d10000.svg)](https://www.wolfram.com/language/)
![Paclet](https://img.shields.io/badge/paclet-ClassicalLieAlgebra%201.0-f57c00.svg)

</div>

`ClassicalLieAlgebra` is a Wolfram Language paclet for working with the classical simple Lie algebras — the special unitary `su(n)`, special orthogonal `so(n)`, and symplectic `sp(2n)` families. It gives you their generators in several standard bases, a Young-tableau toolkit for many-body wavefunctions, and a representation engine that builds any irreducible representation from its highest weight: its dimension, weight system, Casimir eigenvalue, and the explicit generator matrices. Everything is computed in exact arithmetic, and the construction reaches every irrep, including the orthogonal and symplectic **spinor** representations.

<p align="center">
  <img src="pics/su3-levels.png" width="32%" alt="SU(3) as a three-level system">
  <img src="pics/rep11-first-level.png" width="32%" alt="Weight-lowering tree of the (1,1) representation">
  <img src="pics/rep11-transition.png" width="32%" alt="A transition in the (1,1) representation">
</p>

<div align="center">
<a href="#features">Features</a> ·
<a href="#installation">Installation</a> ·
<a href="#quick-start">Quick start</a> ·
<a href="#what-it-covers">What it covers</a> ·
<a href="docs/walkthrough.md">Guided tour</a>
</div>

## Features

- **Three families, one interface.** Write `SU[n]`, `SO[n]`, `Sp[n]`, or the canonical `LieAlgebra["A"|"B"|"C"|"D", rank]`; the sugar forms normalize to it.
- **Every basis.** The standard (defining-representation) generators, the Cartan–Weyl basis, and the Chevalley basis, as exact matrices, with the change of basis between realizations.
- **Young-tableau machinery.** Tensor tableaux, their wavefunctions, inner products, normalization, and orthogonalization of degenerate weight spaces.
- **A representation engine.** `Irrep[g, λ]` gives the dimension (Weyl formula), the weight system with multiplicities (Freudenthal), the Casimir eigenvalue, and the explicit Chevalley generator matrices in the irrep, for every classical irrep (**spinors included**).
- **Exact and tested.** Exact arithmetic throughout, with a `VerificationTest` suite and continuous integration.

## Installation

Install a released build and load it:

```mathematica
PacletInstall["https://github.com/jayren3996/LieAlgebra/releases/download/v1.0.0/ClassicalLieAlgebra-1.0.0.paclet"];
Needs["ClassicalLieAlgebra`"];
```

Or work from a local checkout:

```mathematica
PacletDirectoryLoad["/path/to/LieAlgebra"];
Needs["ClassicalLieAlgebra`"];
```

## Quick start

```mathematica
Needs["ClassicalLieAlgebra`"];

(* defining-representation generators *)
Generators[SU[3]]                  (* the eight Gell-Mann matrices *)
Generators[SU[3], "Chevalley"]     (* <|"Cartan"->{H1,H2}, "Raising"->.., "Lowering"->..|> *)

(* build a representation from its highest weight (Dynkin labels) *)
ir = Irrep[SU[3], {1, 1}];         (* the adjoint / octet *)
RepresentationDimension[ir]        (* 8 *)
CasimirEigenvalue[ir]              (* 6 *)
WeightSystem[ir]                   (* <|{1,1}->1, ..., {0,0}->2, ...|> *)
RepresentationMatrices[ir]         (* explicit generator matrices in the irrep *)

(* the engine reaches so/sp spinor representations too *)
RepresentationDimension[Irrep[SO[5], {0, 1}]]   (* 4: the so(5) spinor *)
```

## What it covers

| | $A_n=\mathfrak{su}(n{+}1)$ | $B_n=\mathfrak{so}(2n{+}1)$ | $C_n=\mathfrak{sp}(2n)$ | $D_n=\mathfrak{so}(2n)$ |
| :--- | :---: | :---: | :---: | :---: |
| Generators · Cartan–Weyl · Chevalley | ✓ | ✓ | ✓ | ✓ |
| Irreps: dimension · weights · Casimir | ✓ | ✓ | ✓ | ✓ |
| Irreps: explicit generator matrices | ✓ | ✓ | ✓ | ✓ |
| Spinor representations | — | ✓ | — | ✓ |

## Learn more

- **[Guided tour of SU(3)](docs/walkthrough.md)** — the long-form walkthrough: the construction of the SU(3) representations step by step, with figures, and the representation engine in depth.
- **`Demo.nb`** — a notebook of further examples.

## Regenerating the figures

Every figure is generated from the package by [`pics/MakeFigures.wls`](pics/MakeFigures.wls), so they stay in step with the code. Rasterizing the typeset matrices and tableaux needs a Wolfram front end:

```sh
wolframscript -file pics/MakeFigures.wls
```

## License

Released under the [MIT License](LICENSE).
