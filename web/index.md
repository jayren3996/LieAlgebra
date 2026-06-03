<div class="cla-hero" markdown>

# ClassicalLieAlgebra

Exact generators, bases, Young tableaux, and irreducible representations of the
classical Lie algebras — **su(n)**, **so(n)**, **sp(2n)** — for the Wolfram Language.

[Guided tour of SU(3)](walkthrough.md){ .md-button .md-button--primary }
[API reference](reference/index.md){ .md-button }

</div>

![The (1,1) "octet" of SU(3), pictured as a three-level system](pics/su3-levels.png){ .center width="520" }

`ClassicalLieAlgebra` is a Wolfram Language paclet for working with the classical
simple Lie algebras — the special unitary `su(n)`, special orthogonal `so(n)`, and
symplectic `sp(2n)` families. It gives you their generators in several standard
bases, a Young-tableau toolkit for many-body wavefunctions, and a representation
engine that builds any irreducible representation from its highest weight: its
dimension, weight system, Casimir eigenvalue, and the explicit generator matrices.

## Why ClassicalLieAlgebra

<div class="grid cards" markdown>

-   :material-check-decagram-outline:{ .lg .middle } __Exact, all the way down__

    ---

    Every root, weight, and matrix element is computed over the rationals — no
    floating point. Results you can trust and verify.

    [:octicons-arrow-right-24: How it's validated](validation.md)

-   :material-set-all:{ .lg .middle } __All four families, one interface__

    ---

    `su`, `so` (both `B` and `D`), and `sp` behind a single API: `SU[n]`,
    `SO[n]`, `Sp[n]`, or the canonical `LieAlgebra[type, rank]`.

    [:octicons-arrow-right-24: The concepts](concepts.md)

-   :material-atom-variant:{ .lg .middle } __Reaches the spinors__

    ---

    Explicit generator matrices for **every** irrep — including the orthogonal and
    symplectic **spinor** representations that no tensor power of the defining
    representation reaches.

    [:octicons-arrow-right-24: Spinors of so(N)](tutorials/spinors.md)

</div>

Coming from another tool? See [how it compares](comparison.md).

## Installation

=== "From a checkout"

    ```mathematica
    PacletDirectoryLoad["/path/to/LieAlgebra"];
    Needs["ClassicalLieAlgebra`"];
    ```

=== "From a release"

    Once a tagged release is published, install the built paclet from its release asset:

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

## Where to next

<div class="grid cards" markdown>

-   :material-school-outline:{ .lg .middle } __Concepts__

    ---

    New to roots, weights, and Dynkin labels? Start with the vocabulary the rest
    of the docs assume — type-agnostic, not just SU(3).

    [:octicons-arrow-right-24: Read the concepts](concepts.md)

-   :material-compass-outline:{ .lg .middle } __Guided tour__

    ---

    Build the SU(3) representations step by step — the Cartan–Weyl and Chevalley
    bases, Young tableaux, and the weight diagram, with figures.

    [:octicons-arrow-right-24: Take the tour](walkthrough.md)

-   :material-flask-outline:{ .lg .middle } __Tutorials__

    ---

    Four worked, runnable walkthroughs: getting started, the representation
    engine, the tableau toolkit, and physics applications.

    [:octicons-arrow-right-24: Browse tutorials](tutorials/index.md)

-   :material-book-open-variant:{ .lg .middle } __API reference__

    ---

    Every public symbol, grouped by theme, with signatures and
    build-evaluated examples.

    [:octicons-arrow-right-24: Open the reference](reference/index.md)

</div>

## License

Released under the [MIT License](https://github.com/jayren3996/LieAlgebra/blob/master/LICENSE).
