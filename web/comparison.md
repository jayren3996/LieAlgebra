# Why this package

If you work with Lie algebras in Mathematica you may already use
[LieART](https://lieart.hepforge.org/) or
[GroupMath](https://renatofonseca.net/groupmath), or reach for
[LiE](http://wwwmathlabo.univ-poitiers.fr/~maavl/LiE/) or
[Sage](https://www.sagemath.org/). They are excellent at what they do —
tabulating dimensions, weights, tensor-product decompositions, and branching
rules, across all types including the exceptional algebras. `ClassicalLieAlgebra`
is not a replacement for them. It fills a different, narrower need.

## What this package is for

Its distinctive feature is **explicit generator matrices**: for any irrep of any
classical algebra, it returns the actual $H_i$, $E_i$, $F_i$ matrices, in an
orthonormal Chevalley basis, in exact arithmetic — and it reaches **every** irrep,
including the orthogonal and symplectic **spinors** that do not live in any tensor
power of the defining representation. Alongside that, it has a **Young-tableau
toolkit** for building and orthogonalizing the many-body wavefunctions that label
representation states.

So if you need the matrices themselves — to build a Hamiltonian, a spin model, or
an explicit symmetry action — rather than a table of dimensions, this is the gap it
fills.

## Capability comparison

| Capability | This package | LieART | GroupMath | LiE / Sage |
| :--- | :---: | :---: | :---: | :---: |
| Classical types A–D | ✓ | ✓ | ✓ | ✓ |
| Exceptional algebras | — | ✓ | ✓ | ✓ |
| Dimension, Casimir | ✓ | ✓ | ✓ | ✓ |
| Weight multiplicities | ✓ | ✓ | ✓ | ✓ |
| Tensor products / branching | — | ✓ | ✓ | ✓ |
| **Explicit generator matrices (orthonormal Chevalley)** | **✓** | — | — | partial |
| **Spinor representation matrices** | **✓** | — | — | — |
| **Young-tableau many-body wavefunctions** | **✓** | — | — | — |
| Exact arithmetic throughout | ✓ | ✓ | ✓ | ✓ |

The comparison is deliberately conservative — the other tools have large feature
sets this table does not enumerate. The point is only where the **bold** rows fall.

## What this package does *not* do

To set expectations honestly:

- **No tensor-product decomposition or branching rules.** It computes dimensions,
  so you can check a decomposition by dimension (see the
  [eightfold-way tutorial](tutorials/physics-applications.md)), but it will not
  hand you $\mathbf 3\otimes\bar{\mathbf 3}=\mathbf 8\oplus\mathbf 1$ as a
  decomposition. For that, keep LieART or GroupMath.
- **No exceptional algebras.** Only the classical families $A$, $B$, $C$, $D$.
- **It is a Wolfram Language paclet**, so it needs Mathematica (LiE and Sage do not).

## See also

- [Validation](validation.md) — how the output is checked.
- [Representations](tutorials/representations.md) and [Spinors](tutorials/spinors.md) — the distinctive features, worked.
