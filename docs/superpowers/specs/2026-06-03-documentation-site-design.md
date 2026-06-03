# ClassicalLieAlgebra documentation site — design

## Goal

Publish a browsable documentation website for the `ClassicalLieAlgebra` paclet on
GitHub Pages, built from the Markdown the repo already has plus a generated,
theme-grouped API reference. This is **Phase 1** of a two-phase documentation
effort; Phase 2 (native Wolfram documentation notebooks for the Documentation
Center / `?Irrep` / F1) is sketched at the end and gets its own spec.

## Motivation

The package is already well documented as source files — a landing `README.md`,
a long-form `docs/walkthrough.md`, four runnable `demos/`, and figures generated
from the package by `pics/MakeFigures.wls`. What it lacks is a single web
destination that ties those together, is discoverable on the open web, is
readable without a Wolfram installation, and gives each public symbol a reference
entry. A static site assembled from the existing Markdown delivers that with no
new prose to maintain and no duplication of the source-of-truth files.

## Decisions

Settled during brainstorming:

- **Scope:** both a web site and native Wolfram docs, done in phases — **web site
  first** (this spec), notebooks second (future spec).
- **Web toolchain:** MkDocs with the Material theme.
- **API reference granularity:** grouped by theme (~5 pages), not one page per
  symbol and not a single long page. The generator still emits a per-symbol block
  within each page so Phase 2 can lift each symbol into its own notebook.
- **Automation:** local build the author runs and deploys; **no GitHub Actions**,
  consistent with the recent removal of CI in favour of local testing.
- **Figures:** a staged build that leaves every source file untouched (see
  *Asset handling*), rather than relocating `pics/` into the docs tree.

## Deliverables

A MkDocs site, its source pages, a reference generator, and a build/deploy
script. The site is served at `https://jayren3996.github.io/LieAlgebra/`.

### Site configuration — `mkdocs.yml` (repo root)

- `site_name`, `site_description`, `repo_url` → the GitHub repo, `site_url` → the
  Pages URL.
- `docs_dir: site-src` — a generated, git-ignored staging directory (see *Asset
  handling*). Pointing MkDocs at a staging dir, rather than the existing `docs/`,
  also keeps `docs/superpowers/` (specs and plans) out of the site automatically.
- `theme: material` with `search`, a light/dark palette toggle, `content.code.copy`,
  and `navigation.sections`/`navigation.top`.
- Markdown extensions: `pymdownx.arithmatex` (generic mode, for `$…$` / `$$…$$`),
  `pymdownx.highlight` + `pymdownx.superfences`, `admonition`, `toc` with
  permalinks.
- `extra_javascript` loads MathJax (CDN) configured for arithmatex.
- `mkdocs build --strict` must pass: no broken links, no warnings.

### Navigation

- **Home** — `index.md`
- **Guided tour** — the existing `docs/walkthrough.md`
- **Demos** — adapted from `demos/README.md`
- **Reference**
  - Overview — `reference/index.md`
  - Algebras & generators — `reference/algebras.md`
  - Root system — `reference/root-system.md`
  - Representations — `reference/representations.md`
  - Young tableaux — `reference/young-tableaux.md`

### Authored site pages (tracked source)

Site-only Markdown lives in a tracked `web/` directory and is copied into the
staging tree at build time:

- `web/index.md` — the landing page, adapted from the README. Authored separately
  (not a copy) because the README's links are repo-relative (`docs/walkthrough.md`,
  `demos/`) and must become site-relative (`walkthrough.md`, `demos.md`) on the
  site. Carries the one-paragraph overview, an install snippet, and links into the
  tour, demos, and reference.
- `web/demos.md` — adapted from `demos/README.md`: the table of the four demos and
  how to run them, linking the `.wls` scripts on GitHub.

### API reference generator — `scripts/generate-reference.wls`

A headless `.wls` using the repo's existing preamble (resolve repo root →
`PacletDirectoryLoad[root]` → `Needs["ClassicalLieAlgebra`"]`).

- Enumerates the public symbols as the names in the `ClassicalLieAlgebra`` context
  that carry a `::usage` string.
- Maps each symbol to one of four theme groups (mapping below).
- For each group, writes `site-src/reference/<group>.md`: a short group intro,
  then one block per symbol — heading (the symbol name), the call signature(s)
  and usage prose taken from the `::usage` string (unescaping `\"`; usage strings
  are confirmed plain text with no box markup), and, where one fits, a short
  example lifted from the corresponding demo.
- Writes `site-src/reference/index.md` — a one-line-per-symbol index linking into
  the group pages.
- **Self-check:** asserts every public symbol lands in exactly one group and every
  group symbol exists, printing `ok` / `FAIL`, so the reference cannot silently
  drift from the API — the same discipline as the demos and `MakeFigures.wls`.

Theme mapping (31 symbols):

| Group page | Symbols |
| :-- | :-- |
| `algebras.md` | `LieAlgebra` `SU` `SO` `Sp` `Generators` `CartanWeyl` `Chevalley` `BasisTransform` |
| `root-system.md` | `Rank` `LieAlgebraDimension` `CartanMatrix` `SimpleRoots` `PositiveRoots` `FundamentalWeights` |
| `representations.md` | `Irrep` `HighestWeight` `RepresentationDimension` `RepresentationMatrices` `CasimirEigenvalue` `WeightSystem` |
| `young-tableaux.md` | `Tableau` `TensorTableau` `TableauForm` `TableauPermute` `TableauNormalization` `TableauOrthogonalization` `TableauDot` `ToTensor` `TensorDot` `TensorNorm` `Psi` |

### Asset handling — staged build

`docs/walkthrough.md` references figures as `../pics/foo.png`, correct for GitHub
(the file is in `docs/`, figures in repo-root `pics/`). MkDocs cannot include
files from outside its `docs_dir`, and a relative path that resolves on GitHub
will not resolve on the site. Rather than duplicate the PNGs or disturb the
figure pipeline, the build **stages** a copy:

`scripts/build-site.sh` assembles `site-src/` (git-ignored) by:

1. running `generate-reference.wls` → `site-src/reference/*.md`;
2. copying `web/index.md` → `site-src/index.md` and `web/demos.md` →
   `site-src/demos.md`;
3. copying `docs/walkthrough.md` → `site-src/walkthrough.md`, rewriting `../pics/`
   → `pics/`;
4. copying `pics/*.png` → `site-src/pics/`;

then runs `mkdocs build --strict` (or `mkdocs serve` for preview). Every source
file is left untouched: the figure pipeline and GitHub's own rendering keep
working unchanged.

### Build & deploy

- `scripts/build-site.sh` — assemble + `mkdocs build --strict`; a `serve` flag for
  local preview at `localhost:8000`.
- Deploy: `mkdocs gh-deploy --clean`, which builds and pushes the HTML to the
  `gh-pages` branch.
- **One-time manual step (author):** enable GitHub Pages in repo settings with
  source = `gh-pages` branch. Printed by the build script's usage output; not
  automatable from here.
- `.gitignore`: add `site-src/` and `site/`.

### Dependencies (one-time, local)

- `pip install mkdocs-material` (brings in `pymdown-extensions`).
- MathJax loaded from CDN; nothing added to the paclet.

## Verification

- `mkdocs build --strict` exits clean — no broken internal links, no warnings.
- The generator's self-check prints `ok` for the symbol-coverage assertions.
- `mkdocs serve` spot-check: math renders, `mathematica` code blocks are
  highlighted, walkthrough figures load, dark/light toggle and search work.

## Out of scope (this phase)

- **Native Wolfram documentation notebooks** — Phase 2. A `Documentation/English/`
  tree (a reference page per symbol, a guide page, tech notes) built with
  `PacletDocumentationBuild` so `?Irrep` / F1 work and the paclet is ready for the
  Wolfram Paclet Repository, optionally exporting HTML back into this site. The
  reference generator here emits per-symbol blocks specifically so that content
  can seed the notebook skeletons. Phase 2 gets its own spec and plan.
- GitHub Actions / CI to auto-build or auto-deploy the site.
- Versioned docs (e.g. `mike`), custom domain, analytics.
- Any change to the paclet's public API or to the figure pipeline.
