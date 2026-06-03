<div align="center">

# ClassicalLieAlgebra

**Exact generators, bases, Young tableaux, and irreducible representations of the classical Lie algebras — for the Wolfram Language.**

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Wolfram Language](https://img.shields.io/badge/Wolfram%20Language-13.0%2B-d10000.svg)](https://www.wolfram.com/language/)
![Paclet](https://img.shields.io/badge/paclet-ClassicalLieAlgebra%201.0-f57c00.svg)

</div>

`ClassicalLieAlgebra` is a Wolfram Language paclet for working with the classical simple Lie algebras — the special unitary `su(n)`, special orthogonal `so(n)`, and symplectic `sp(2n)` families. It gives you their generators in several standard bases, a Young-tableau toolkit for many-body wavefunctions, and a representation engine that builds any irreducible representation from its highest weight: its dimension, weight system, Casimir eigenvalue, and the explicit generator matrices. Everything is computed in exact arithmetic, and the construction reaches every irrep, including the orthogonal and symplectic **spinor** representations.

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
- **Exact and tested.** Exact arithmetic throughout, with a `VerificationTest` suite.

## Installation

Work from a local checkout:

```mathematica
PacletDirectoryLoad["/path/to/LieAlgebra"];
Needs["ClassicalLieAlgebra`"];
```

Once a tagged release is published, you can install the built paclet from its release asset instead:

```mathematica
PacletInstall["https://github.com/jayren3996/LieAlgebra/releases/download/vX.Y.Z/ClassicalLieAlgebra-X.Y.Z.paclet"];
Needs["ClassicalLieAlgebra`"];
```

## Quick start

```mathematica
Needs["ClassicalLieAlgebra`"];

(* defining-representation generators *)
Generators[SU[3]]                  (* the 8 su(3) generators, T_a = lambda_a/2 *)
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
- **[Runnable demos](demos/)** — four self-contained scripts you can run straight from a checkout: a getting-started tour, the representation engine, the Young-tableau toolkit, and physics applications.

## License

Released under the [MIT License](LICENSE).
