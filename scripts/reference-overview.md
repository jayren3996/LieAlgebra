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

### Types and families

`SU[n]`, `SO[n]`, and `Sp[n]` all normalize to `LieAlgebra[type, rank]`, so the two
forms are interchangeable:

| Type | Algebra | Constructor |
| :-- | :-- | :-- |
| `"A"` | $\mathfrak{su}(n+1)$ | `SU[n+1]` |
| `"B"` | $\mathfrak{so}(2n+1)$ | `SO[2n+1]` |
| `"C"` | $\mathfrak{sp}(2n)$ | `Sp[2n]` |
| `"D"` | $\mathfrak{so}(2n)$ | `SO[2n]` |

### Weights and Dynkin labels

- **Weights** are given as **Dynkin labels** — a list of `Rank[g]` integers. The
  highest weight passed to `Irrep` must be **non-negative**; weights returned by
  `WeightSystem` may have negative entries.
- **Simple roots follow the Bourbaki labelling** of the Dynkin diagram. For
  $B_n=\mathfrak{so}(2n{+}1)$ the short simple root is the last node, so the spinor
  is `{0, ..., 0, 1}`; for $C_n=\mathfrak{sp}(2n)$ the long root is the last node;
  for $D_n=\mathfrak{so}(2n)$ the two half-spinor nodes are the last two.
- Roots and weights returned by the root-system functions are in the Euclidean
  (orthonormal) basis.

### Normalization

- **Root length:** long roots have squared length 2 (the mathematicians'
  normalization).
- **Casimir:** `CasimirEigenvalue` is $(\lambda, \lambda + 2\rho)$ in that
  normalization — twice the physics value, e.g. the su(2) spin-½ value is `3/2`
  against the physics `3/4`.
- **Defining generators:** `Generators[g]` is normalized to
  $\operatorname{Tr}(T_a T_b)=\tfrac12\delta_{ab}$; for su(3) these are the
  Gell-Mann matrices $\lambda_a/2$.
- **Representation matrices:** `RepresentationMatrices` returns an orthonormal basis
  with $E_i=\operatorname{ConjugateTranspose}(F_i)$ and $H_i$ diagonal.

### Associations and exact arithmetic

- `CartanWeyl`, `Chevalley`, and `RepresentationMatrices` all return an association
  keyed `"Cartan"`, `"Raising"`, `"Lowering"`.
- Everything is computed in **exact arithmetic** over the rationals; radicals enter
  only in the final orthonormal rescaling of representation matrices.
- Invalid arguments emit a message and return `$Failed`.

> **The Cartan matrix is transposed relative to some references.**
>
> `CartanMatrix[g]` returns $A_{ij}=2(\alpha_i,\alpha_j)/(\alpha_j,\alpha_j)$,
> equivalently $[H_i,E_j]=A_{ji}E_j$. For the non-simply-laced types ($B$, $C$)
> this is the **transpose** of the matrix some references (e.g. LieART) print —
> for instance `CartanMatrix[SO[5]]` is `{{2, -2}, {-1, 2}}`. The package is
> internally consistent in this convention.
