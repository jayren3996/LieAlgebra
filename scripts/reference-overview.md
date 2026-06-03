# API reference

`ClassicalLieAlgebra` exposes 31 public symbols across four areas. This page
introduces how they fit together; each symbol's full reference — signature,
arguments, and runnable examples — is on the themed pages below.

## The mental model

Everything starts from an **algebra**, built with `SU`, `SO`, `Sp`, or the
canonical `LieAlgebra[type, rank]` they normalize to:

```mathematica
g = SU[3];
```

From an algebra you can go in three directions:

- **Read its root-system data** — `Rank`, `LieAlgebraDimension`, `CartanMatrix`,
  `SimpleRoots`, `PositiveRoots`, `FundamentalWeights`.
- **Get its generators** in any of three bases — the defining-representation
  `Generators[g]`, the `CartanWeyl[g]` basis, or the `Chevalley[g]` basis;
  `BasisTransform` relates the two so/sp realizations.
- **Build an irreducible representation** from a highest weight with
  `Irrep[g, w]`, then read its `RepresentationDimension`, `WeightSystem`,
  `CasimirEigenvalue`, and explicit `RepresentationMatrices`. `HighestWeight`
  recovers `w`.

```mathematica
ir = Irrep[SU[3], {1, 1}];          (* the adjoint / octet *)
RepresentationDimension[ir]          (* 8 *)
```

A separate **Young-tableau toolkit** (`Tableau`, `TensorTableau`, `Psi`, and the
`Tableau…` / `Tensor…` operations) builds and manipulates the many-body wave
functions that label representation states — the symmetrizer, inner products,
normalization, and orthogonalization of degenerate weight spaces.

## Conventions

- **Weights** are given as **Dynkin labels** — a list of `Rank[g]` integers.
  Roots and weights returned by the root-system functions are in the Euclidean
  (orthonormal) basis.
- Everything is computed in **exact arithmetic**.
- `SU[n]`, `SO[n]`, and `Sp[n]` all normalize to `LieAlgebra[type, rank]`, so the
  two forms are interchangeable.
