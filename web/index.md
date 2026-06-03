# ClassicalLieAlgebra

**Exact generators, bases, Young tableaux, and irreducible representations of the classical Lie algebras — for the Wolfram Language.**

`ClassicalLieAlgebra` is a Wolfram Language paclet for working with the classical simple Lie algebras — the special unitary `su(n)`, special orthogonal `so(n)`, and symplectic `sp(2n)` families. It gives you their generators in several standard bases, a Young-tableau toolkit for many-body wavefunctions, and a representation engine that builds any irreducible representation from its highest weight: its dimension, weight system, Casimir eigenvalue, and the explicit generator matrices. Everything is computed in exact arithmetic, and the construction reaches every irrep, including the orthogonal and symplectic **spinor** representations.

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

- **[Guided tour of SU(3)](walkthrough.md)** — the long-form walkthrough: the construction of the SU(3) representations step by step, with figures, and the representation engine in depth.
- **[Runnable demos](demos.md)** — four self-contained scripts you can run straight from a checkout: a getting-started tour, the representation engine, the Young-tableau toolkit, and physics applications.
- **[API reference](reference/index.md)** — every public symbol, grouped by theme.

## License

Released under the [MIT License](https://github.com/jayren3996/LieAlgebra/blob/master/LICENSE).
