# Physics applications

Three vignettes showing the algebras at work in physics: SU(3) flavour and the
eightfold way, SU(2) and angular momentum, and a direct numeric check that a built
representation really obeys the algebra. Mirrors
[`demos/04-physics-applications.wls`](https://github.com/jayren3996/LieAlgebra/blob/master/demos/04-physics-applications.wls).

```mathematica
Needs["ClassicalLieAlgebra`"];
dim[g_, w_] := RepresentationDimension[Irrep[g, w]];
```

## SU(3) flavour and the eightfold way

Gell-Mann and Ne'eman organized the hadrons into SU(3) flavour multiplets. Each
multiplet is an irrep, labelled by its Dynkin labels:

```mathematica
dim[SU[3], {1, 0}]   (* quark triplet      3   = (1,0) *)
dim[SU[3], {0, 1}]   (* antiquark triplet  3bar = (0,1) *)
dim[SU[3], {1, 1}]   (* meson/baryon octet 8   = (1,1) *)
dim[SU[3], {3, 0}]   (* baryon decuplet    10  = (3,0) *)
```

```
3
3
8
10
```

The octet's weights are the eight states' (isospin, hypercharge) quantum numbers.
The doubly-degenerate central weight `{0, 0}` is the $\pi^0$/$\eta$ (or
$\Sigma^0$/$\Lambda$) pair:

```mathematica
WeightSystem[Irrep[SU[3], {1, 1}]]
```

```
<|{1, 1} -> 1, {-1, 2} -> 1, {2, -1} -> 1, {0, 0} -> 2,
  {-2, 1} -> 1, {1, -2} -> 1, {-1, -1} -> 1|>
```

The hadron multiplets follow from tensoring quarks. The package computes
dimensions (it does not decompose tensor products), and the dimensions are
consistent with the classic decompositions:

$$
\mathbf{3}\otimes\bar{\mathbf 3}=\mathbf 8\oplus\mathbf 1,\qquad
\mathbf 3\otimes\mathbf 3\otimes\mathbf 3=\mathbf{10}\oplus\mathbf 8\oplus\mathbf 8\oplus\mathbf 1
$$

```mathematica
dim[SU[3], {1, 0}] dim[SU[3], {0, 1}]   (* 9  = 8 + 1 *)
dim[SU[3], {1, 0}]^3                     (* 27 = 10 + 8 + 8 + 1 *)
```

## SU(2) and angular momentum

The `su(2)` irrep with Dynkin label `{2j}` is the spin-$j$ multiplet: dimension
$2j+1$, and quadratic Casimir $2j(j+1)$ in this package's normalization.

| $j$ | Dynkin | `dim` | `Casimir` | $2j(j+1)$ |
| :---: | :---: | :---: | :---: | :---: |
| 1/2 | `{1}` | 2 | 3/2 | 3/2 |
| 1 | `{2}` | 3 | 4 | 4 |
| 3/2 | `{3}` | 4 | 15/2 | 15/2 |
| 2 | `{4}` | 5 | 12 | 12 |

```mathematica
Table[{#/2, dim[SU[2], {#}], CasimirEigenvalue[Irrep[SU[2], {#}]]} &[n], {n, 1, 4}]
```

(The Casimir is twice the physics value $j(j+1)$, because of the long-root²=2
normalization.)

## A built representation really obeys the algebra

`RepresentationMatrices` returns explicit numbers, so we can verify directly that
they satisfy the defining Chevalley relations of `su(3)`,
$[E_i,F_j]=\delta_{ij}H_i$ and $[H_i,E_j]=A_{ji}E_j$:

```mathematica
bracket[a_, b_] := a . b - b . a;

chevalleyHolds[g_, w_] := Module[
   {m = RepresentationMatrices[Irrep[g, w]], a = CartanMatrix[g], r = Rank[g], z},
   z = 0 m["Cartan"][[1]];
   (And @@ Flatten@Table[
       bracket[m["Raising"][[i]], m["Lowering"][[j]]] == If[i == j, m["Cartan"][[i]], z],
       {i, r}, {j, r}]) &&
   (And @@ Flatten@Table[
       bracket[m["Cartan"][[i]], m["Raising"][[j]]] == a[[j, i]] m["Raising"][[j]],
       {i, r}, {j, r}])];

chevalleyHolds[SU[3], {1, 1}]   (* True: the octet matrices form su(3) *)
```

The same engine and the same check work for spinors and for symplectic algebras:

```mathematica
chevalleyHolds[SO[5], {0, 1}]   (* True: the so(5) spinor *)
chevalleyHolds[Sp[4], {1, 0}]   (* True: the sp(4) defining rep *)
```

## Next

- How the irreps are built: [Representations](representations.md).
- The state labels behind the octet: [Young tableaux](young-tableaux.md).
- Reference: [`RepresentationDimension`](../reference/representations.md#representationdimension) ·
  [`CasimirEigenvalue`](../reference/representations.md#casimireigenvalue) ·
  [`RepresentationMatrices`](../reference/representations.md#representationmatrices).
