# Validation

`ClassicalLieAlgebra` is built to be trusted: every quantity is computed in exact
arithmetic, and a `VerificationTest` suite checks the output against known results
and against the defining relations of the algebras themselves. This page says what
is checked, and how to re-run it.

## Exact arithmetic

All root-system data, weight multiplicities, dimensions, and Casimir eigenvalues
are computed over the rationals $\mathbb{Q}$ — no floating point, no `N`. The
explicit representation matrices are exact too; the only radicals that appear are
the $1/\sqrt{\cdot}$ factors of the final orthonormal rescaling, and they are
carried symbolically. So a returned `8`, `5/2`, or `1/Sqrt[6]` is the exact value,
not a numerical approximation.

## What the suite checks

The suite lives in
[`Tests/`](https://github.com/jayren3996/LieAlgebra/tree/master/Tests) — one
`.wlt` file per subsystem (loading, algebras, su/so/sp generators, weights,
representations, Young tableaux). The representation tests
([`Tests/Representations.wlt`](https://github.com/jayren3996/LieAlgebra/blob/master/Tests/Representations.wlt))
are the core of the correctness story:

| Checked | Against | Examples |
| :--- | :--- | :--- |
| **Dimensions** | known irrep dimensions | su(3) **3**, **8**; so(5) vector **5**, spinor **4**; sp(4) **4**, **5**; so(7) spinor **8** |
| **Casimir eigenvalues** | known values | su(2) spin-½ `3/2`; su(3) adjoint `6`; so(5) spinor `5/2`; so(7) spinor `21/4` |
| **Weight multiplicities** | the dimension | $\sum$ multiplicities $=$ `RepresentationDimension` for su/so/sp, incl. spinors |
| **Interior multiplicities** | known degeneracies | su(3) adjoint zero-weight mult **2**; su(4) adjoint and su(3) **27**-plet mult **3** |
| **Chevalley relations** | $[E_i,F_j]=\delta_{ij}H_i$, $[H_i,E_j]=A_{ji}E_j$ | hold for su(3) adjoint, so(5)/so(7) spinors, so(6)/so(8) half-spinors, sp(4) |
| **Serre relations** | $\operatorname{ad}(E_i)^{1-A_{ji}}(E_j)=0$ | su(3), so(5) spinor, sp(4) |
| **Orthonormal basis** | the documented contract | $E_i=F_i^\dagger$ and $H_i$ Hermitian, for every tested irrep |
| **Spinor dimensions** | $2^n$ (type B), $2^{n-1}$ (type D) | so(5)→4, so(7)→8, so(6)→4, so(8) half-spinors→8 |

That the **built matrices satisfy the Chevalley and Serre relations** is the
strongest single check: it verifies, for each irrep, that the explicit numbers
really form a representation of the algebra — spinors and symplectic cases
included. You can reproduce it yourself with the `chevalleyHolds` one-liner in the
[physics-applications tutorial](tutorials/physics-applications.md).

## Re-running it

```sh
wolframscript -file scripts/runTests.wls
```

It runs every `Tests/*.wlt` suite, prints a `PASS`/`FAIL` line per file, and exits
non-zero if anything fails. The [demos](tutorials/index.md), the
[reference examples](reference/index.md) (evaluated at build time), and the
walkthrough figures are all self-checking too, so the documentation cannot drift
from the code without a build failing.

## See also

- [Conventions](concepts.md#conventions) — the normalizations the test values assume.
- [Why this package](comparison.md) — where this fits among the alternatives.
- [Physics applications](tutorials/physics-applications.md) — the algebra-closure check, run live.
