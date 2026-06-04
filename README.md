<div align="center">

# ClassicalLieAlgebra

**Exact classical Lie algebras in the Wolfram Language.**

Generators, Cartan-Weyl and Chevalley bases, Young tableaux, and explicit irreducible representations for `su(n)`, `so(n)`, and `sp(2n)`.

[![Docs](https://img.shields.io/badge/docs-latest-9558B2.svg)](https://jayren3996.github.io/LieAlgebra/) [![Wolfram Language](https://img.shields.io/badge/Wolfram%20Language-13.0%2B-d10000.svg)](https://www.wolfram.com/language/) [![Paclet](https://img.shields.io/badge/paclet-ClassicalLieAlgebra%201.0-f57c00.svg)](PacletInfo.wl) [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

</div>

---

`ClassicalLieAlgebra` is a Wolfram Language paclet for computations with the
classical simple Lie algebras. It gives one interface for the `A`, `B`, `C`, and
`D` families, exact generator matrices in standard bases, a Young-tableau toolkit
for many-body wavefunctions, and a representation engine that builds irreducible
representations from highest weights.

The package is written for calculations where conventions matter: root-system
data, Cartan matrices, Chevalley generators, weight multiplicities, Casimir
eigenvalues, and explicit representation matrices are all exposed directly and
computed in exact arithmetic.

## ✨ Features

|  |  |
| --- | --- |
| 🧩 **One interface for A/B/C/D** | Use `SU[n]`, `SO[n]`, `Sp[n]`, or the canonical `LieAlgebra["A"|"B"|"C"|"D", rank]`; shorthand constructors normalize to the same algebra object. |
| 🧱 **Generator bases** | `Generators`, `CartanWeyl`, and `Chevalley` return exact defining-representation matrices, with consistent conventions across all classical families. |
| 📐 **Root-system data** | `Rank`, `LieAlgebraDimension`, `CartanMatrix`, `SimpleRoots`, `PositiveRoots`, and `FundamentalWeights` expose the structural data used by the rest of the package. |
| 🧮 **Irreducible representations** | `Irrep[g, λ]` gives dimensions, weight systems with multiplicities, Casimir eigenvalues, and explicit Chevalley generator matrices. |
| 🧬 **Spinors included** | Orthogonal spinor representations are part of the same highest-weight workflow, not a separate special case. |
| 🔳 **Young tableaux** | Tensor tableaux, wavefunctions, inner products, normalization, symmetrization, and orthogonalization for degenerate weight spaces. |
| ✅ **Exact and tested** | Arithmetic is exact throughout, with a Wolfram `VerificationTest` suite and self-checking demo scripts. |

## 📦 Installation

From a checkout:

```mathematica
PacletDirectoryLoad["/path/to/LieAlgebra"];
Needs["ClassicalLieAlgebra`"];
```

From a release asset, once a tagged paclet release is available:

```mathematica
PacletInstall["https://github.com/jayren3996/LieAlgebra/releases/download/vX.Y.Z/ClassicalLieAlgebra-X.Y.Z.paclet"];
Needs["ClassicalLieAlgebra`"];
```

## 🚀 Quick Start

Build the `su(3)` adjoint representation and inspect its basic invariants:

```mathematica
Needs["ClassicalLieAlgebra`"];

g = SU[3];

Generators[g]                 (* defining-representation generators *)
Chevalley[g]                  (* <|"Cartan" -> ..., "Raising" -> ..., "Lowering" -> ...|> *)

ir = Irrep[g, {1, 1}];        (* the adjoint / octet *)

RepresentationDimension[ir]   (* 8 *)
CasimirEigenvalue[ir]         (* 6 *)
WeightSystem[ir]              (* weights with multiplicities *)
RepresentationMatrices[ir]    (* explicit matrices in the irrep *)
```

The same representation workflow reaches orthogonal spinors:

```mathematica
RepresentationDimension[Irrep[SO[5], {0, 1}]]   (* 4 *)
```

## 🧭 Choosing a Workflow

| If you want to ... | Start with |
| --- | --- |
| Construct a classical algebra | `SU[n]`, `SO[n]`, `Sp[n]`, or `LieAlgebra[type, rank]` |
| Read root-system conventions | [`Concepts`](https://jayren3996.github.io/LieAlgebra/concepts/) |
| Work with defining generators | [`Generators`](https://jayren3996.github.io/LieAlgebra/reference/algebras/#generators), [`CartanWeyl`](https://jayren3996.github.io/LieAlgebra/reference/algebras/#cartanweyl), [`Chevalley`](https://jayren3996.github.io/LieAlgebra/reference/algebras/#chevalley) |
| Build an irrep from a highest weight | [`Irrep`](https://jayren3996.github.io/LieAlgebra/reference/representations/#irrep) |
| Inspect dimensions, weights, and Casimirs | [`RepresentationDimension`](https://jayren3996.github.io/LieAlgebra/reference/representations/#representationdimension), [`WeightSystem`](https://jayren3996.github.io/LieAlgebra/reference/representations/#weightsystem), [`CasimirEigenvalue`](https://jayren3996.github.io/LieAlgebra/reference/representations/#casimireigenvalue) |
| Use Young-tableau states | [`Young tableaux`](https://jayren3996.github.io/LieAlgebra/tutorials/young-tableaux/) |
| See the SU(3) construction end to end | [`Guided tour of SU(3)`](https://jayren3996.github.io/LieAlgebra/walkthrough/) |

## 📚 Documentation

Full documentation lives at
**[jayren3996.github.io/LieAlgebra](https://jayren3996.github.io/LieAlgebra/)**.

- [Getting started](https://jayren3996.github.io/LieAlgebra/tutorials/getting-started/) — constructors, root-system data, and generator bases.
- [Representations](https://jayren3996.github.io/LieAlgebra/tutorials/representations/) — highest weights, weight systems, Casimirs, and explicit matrices.
- [Spinors of `so(N)`](https://jayren3996.github.io/LieAlgebra/tutorials/spinors/) — orthogonal spinor representations.
- [Young tableaux](https://jayren3996.github.io/LieAlgebra/tutorials/young-tableaux/) — tensor tableaux and wavefunction operations.
- [Physics applications](https://jayren3996.github.io/LieAlgebra/tutorials/physics-applications/) — SU(3) flavor, SU(2) spin, and algebra checks.
- [API reference](https://jayren3996.github.io/LieAlgebra/reference/) — every public symbol with evaluated examples.

## 🧪 Runnable Demos

The [`demos/`](demos/) directory contains self-checking `.wls` scripts. They can
be run directly from a checkout:

```sh
wolframscript -file demos/01-getting-started.wls
wolframscript -file demos/02-representations.wls
wolframscript -file demos/03-young-tableaux.wls
wolframscript -file demos/04-physics-applications.wls
```

Each script loads the local paclet with `PacletDirectoryLoad` and finishes with
`[ok]` / `[FAIL]` checks, so the demos double as smoke tests for the public API.

## ✅ Coverage

|  | $A_n=\mathfrak{su}(n{+}1)$ | $B_n=\mathfrak{so}(2n{+}1)$ | $C_n=\mathfrak{sp}(2n)$ | $D_n=\mathfrak{so}(2n)$ |
| :--- | :---: | :---: | :---: | :---: |
| Generators · Cartan-Weyl · Chevalley | ✓ | ✓ | ✓ | ✓ |
| Root-system data | ✓ | ✓ | ✓ | ✓ |
| Irreps: dimension · weights · Casimir | ✓ | ✓ | ✓ | ✓ |
| Irreps: explicit generator matrices | ✓ | ✓ | ✓ | ✓ |
| Spinor representations | — | ✓ | — | ✓ |

## 🗂 Repository Map

| Path | What lives there |
| --- | --- |
| [`Kernel/`](Kernel/) | Wolfram package implementation |
| [`Tests/`](Tests/) | Wolfram and documentation-staging tests |
| [`demos/`](demos/) | Runnable tutorial scripts |
| [`docs/`](docs/) | Long-form notes and walkthrough material |
| [`web/`](web/) | Source pages for the web documentation |
| [`scripts/`](scripts/) | Test, reference-generation, and documentation build scripts |
| [`pics/`](pics/) | Generated figures used by the README and docs |

## 🛠 Building the Documentation

The web documentation is built with Sphinx and Furo:

```sh
python3 -m pip install -r requirements-docs.txt
scripts/build-docs.sh
```

For a local preview, run `scripts/build-docs.sh serve` and open
`http://localhost:8000/`.

## 📄 License

[MIT](LICENSE) © Jie Ren and contributors.
