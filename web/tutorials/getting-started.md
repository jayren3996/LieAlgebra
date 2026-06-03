# Getting started

A tour of the basics: constructing the classical algebras, reading off their
root-system data, and looking at the defining generators in the three standard
bases. This mirrors
[`demos/01-getting-started.wls`](https://github.com/jayren3996/LieAlgebra/blob/master/demos/01-getting-started.wls);
every output below is exactly what the script prints.

```mathematica
Needs["ClassicalLieAlgebra`"];
```

## Constructing the algebras

The classical simple Lie algebras come in four families. Each has a canonical form
`LieAlgebra[type, rank]`, and a familiar-name shorthand the constructors normalize
to:

```mathematica
{SU[3], SO[5], SO[6], Sp[4]}
```

```
{LieAlgebra["A", 2], LieAlgebra["B", 2], LieAlgebra["D", 3], LieAlgebra["C", 2]}
```

So `SO[5]` is $B_2$, `SO[6]` is $D_3$, `Sp[4]` is $C_2$. The shorthand and the
canonical form are the same object:

```mathematica
SU[3] === LieAlgebra["A", 2]
```

```
True
```

## Root-system data

Everything about an algebra's structure comes from its root system. Take
`SO[5]` $=B_2$:

```mathematica
g = SO[5];
Rank[g]                  (* 2 *)
LieAlgebraDimension[g]   (* 10 *)
CartanMatrix[g]          (* {{2, -2}, {-1, 2}} *)
SimpleRoots[g]           (* {{1, -1}, {0, 1}} *)
PositiveRoots[g]         (* {{1, -1}, {1, 1}, {1, 0}, {0, 1}} *)
FundamentalWeights[g]    (* {{1, 0}, {1/2, 1/2}} *)
```

The dimension is `(# positive roots) × 2 + rank` $= 4\times2 + 2 = 10$, and the
half-integer fundamental weight `{1/2, 1/2}` is the one that produces the
4-dimensional **spinor** — see the [Representations tutorial](representations.md).

!!! note "The Cartan matrix is transposed relative to LieART"

    `CartanMatrix[SO[5]]` is `{{2, -2}, {-1, 2}}` — the package uses
    $A_{ij}=2(\alpha_i,\alpha_j)/(\alpha_j,\alpha_j)$ with
    $[H_i,E_j]=A_{ji}E_j$, which for the non-simply-laced types ($B$, $C$) is the
    **transpose** of the matrix some references print. It is internally
    consistent; see the [conventions](../concepts.md#conventions).

## Generators in three bases

`Generators[g]` returns the defining-representation generators. For `su(3)` these
are the eight generators $T_a=\lambda_a/2$ (the normalized Gell-Mann matrices),
each $3\times 3$:

```mathematica
gens = Generators[SU[3]];
Length[gens]          (* 8 *)
Dimensions[gens[[1]]] (* {3, 3} *)
gens[[1]]             (* lambda_1 / 2 *)
```

```
{{0, 1/2, 0}, {1/2, 0, 0}, {0, 0, 0}}
```

$$
\frac{\lambda_1}{2}=\begin{pmatrix}0&\tfrac12&0\\[2pt]\tfrac12&0&0\\[2pt]0&0&0\end{pmatrix}
$$

`CartanWeyl[g]` and `Chevalley[g]` return an association keyed by `"Cartan"`,
`"Raising"`, `"Lowering"`:

```mathematica
cw = CartanWeyl[SU[3]];
Keys[cw]            (* {"Cartan", "Raising", "Lowering"} *)
Length /@ Values[cw] (* {2, 3, 3}: 2 Cartan generators, 3 raising, 3 lowering *)
```

The Chevalley basis is the convenient one for building representations — one
$H_i$, $E_i$, $F_i$ per simple root:

```mathematica
ch = Chevalley[SU[3]];
ch["Cartan"][[1]]   (* H1 *)
ch["Raising"][[1]]  (* E1 *)
```

$$
H_1=\begin{pmatrix}1&0&0\\0&-1&0\\0&0&0\end{pmatrix},\qquad
E_1=\begin{pmatrix}0&1&0\\0&0&0\\0&0&0\end{pmatrix}
$$

These really are generators of the same algebra — for instance $[E_1,F_1]=H_1$:

```mathematica
ch["Raising"][[1]] . ch["Lowering"][[1]] - ch["Lowering"][[1]] . ch["Raising"][[1]] == ch["Cartan"][[1]]
```

```
True
```

For `so`/`sp`, the invariant bilinear form has two realizations, and
`BasisTransform[g]` is the matrix relating them. For `SO[5]` it is $5\times 5$:

```mathematica
Dimensions[BasisTransform[SO[5]]]   (* {5, 5} *)
```

## Next

- The vocabulary behind all of this: [Concepts](../concepts.md).
- Build representations from these generators: [Representations](representations.md).
- Reference entries: [`Generators`](../reference/algebras.md#generators) ·
  [`CartanMatrix`](../reference/root-system.md#cartanmatrix) ·
  [`SimpleRoots`](../reference/root-system.md#simpleroots).
