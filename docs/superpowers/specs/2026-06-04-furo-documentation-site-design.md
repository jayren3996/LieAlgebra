# ClassicalLieAlgebra Furo documentation site - design

## Goal

Replace the current MkDocs Material web documentation with a cleaner, more
elegant Sphinx + Furo site while preserving the current documentation content
and the build-time-evaluated Mathematica examples.

This phase is focused on the public web documentation pages. Native Wolfram
Documentation Center notebooks remain out of scope.

## Motivation

The current documentation already has useful content: a homepage, concepts,
tutorials, a guided tour, validation notes, and a generated grouped API
reference. The main issue is presentation. The site still feels like a default
MkDocs Material site with a custom purple skin, and Mathematica examples render
as disconnected raw input and raw output blocks rather than as Wolfram-style
notebook examples.

Sphinx with Furo is a better fit for the desired web presentation:

- Furo has a quieter, more refined documentation template.
- Sphinx gives more direct control over page templates, roles, directives, and
  reusable components.
- MyST Markdown lets the existing Markdown content migrate without rewriting it
  as reStructuredText.
- The existing Wolfram reference generator can be adapted rather than replaced.

## Decisions

- Use **Sphinx + Furo** as the web documentation stack.
- Use **MyST Parser** so pages can stay in Markdown.
- Keep examples evaluated by `wolframscript` during the documentation build.
- Render examples as notebook-like input/output components, not as adjacent
  plain code blocks.
- Keep the grouped reference-page model for the first migration.
- Keep GitHub Pages deployment from a static build output.
- Do not adopt Documenter.jl for this phase.
- Do not build native Wolfram Documentation Center notebooks in this phase.

## Site architecture

Create a dedicated generated Sphinx source tree:

```text
docs-site/
  conf.py
  index.md
  concepts.md
  tutorials/
  reference/
  _static/
    cla.css
    mathjax.js
  _templates/
```

The existing source material remains the source of truth:

- `web/*.md` and `web/tutorials/*.md` provide authored web pages.
- `docs/walkthrough.md` provides the guided tour.
- `pics/*.png` provides generated figures.
- `scripts/reference-details.wl` and `scripts/reference-overview.md` provide the
  curated reference content.
- `scripts/generate-reference.wls` remains responsible for evaluating examples
  and producing reference pages, but it targets the Sphinx source tree.

The build script stages or generates the Sphinx tree without mutating authored
content.

## Visual direction

Use Furo's default structure and customize lightly:

- no gradient hero;
- restrained color palette in light and dark modes;
- a Wolfram/math accent color used only for links, active navigation, and small
  UI highlights;
- clean homepage with concise package positioning, a real package figure, and
  compact entry points;
- wider, calmer reference pages with clearer headings and symbol sections;
- improved spacing around tables, math, figures, admonitions, and examples.

The site should feel like technical/scientific package documentation, not a
marketing page.

## Mathematica example rendering

Generated examples should render as a single semantic block with paired cells:

```text
In[1]:= RepresentationDimension[Irrep[SU[3], {1, 1}]]

Out[1]= 8
```

The generator should produce MyST-friendly HTML or directive markup that Sphinx
can pass through safely. The rendered block should include:

- Mathematica syntax highlighting for input;
- a copy button for input when available through Sphinx tooling or a small local
  enhancement;
- a visually distinct output region;
- sequential input/output labels within each page or symbol entry;
- horizontal scrolling for wide matrices, associations, and long lists;
- image support for examples whose useful output is a committed or generated
  figure.

Output formatting tiers:

- small scalar/list/association output: compact preformatted Wolfram input form;
- large output: scrollable preformatted block;
- matrices: preformatted by default, with optional TeX rendering later if it is
  clearly more readable;
- front-end-only display forms such as tableau renderings: validated by the
  kernel and displayed as an image or short note.

## Content migration

Migrate the existing web documentation pages first:

- Home
- Concepts
- Guided tour
- Tutorials
- Reference overview and grouped reference pages
- Validation
- Comparison
- Changelog
- Citing
- Glossary

The first migration should avoid rewriting the prose except where links, figure
paths, or heading levels need to change for Sphinx.

## Build and deployment

Add a Sphinx build script:

```text
scripts/build-docs.sh
```

The existing `scripts/build-site.sh` remains in place during the migration for
comparison, but the Furo site is built through `scripts/build-docs.sh`.

The script should:

1. create or refresh the Sphinx staging/source tree;
2. copy authored Markdown pages;
3. copy figures and static assets;
4. run the reference generator;
5. run `sphinx-build -W -b html docs-site site`;
6. support a local preview mode.

Deployment uses the existing GitHub Pages `gh-pages` branch model: build static
HTML into `site/`, then publish that output to `gh-pages` after local
verification passes.

## Dependencies

Expected Python dependencies:

- `sphinx`
- `furo`
- `myst-parser`
- `sphinx-copybutton` if it behaves well with Mathematica blocks

Expected Wolfram dependency:

- a working `wolframscript`, using the existing app-bundled fallback behavior if
  needed.

## Verification

The migration is complete when:

- `sphinx-build -W -b html docs-site site` exits cleanly;
- every existing top-level documentation page has a Sphinx equivalent;
- generated reference pages cover every public symbol currently covered by the
  MkDocs reference;
- every generated example still evaluates cleanly;
- Mathematica input/output examples render as paired notebook-style blocks;
- MathJax equations render correctly;
- figures in the guided tour and tutorials load correctly;
- light and dark modes are readable;
- local spot checks confirm that homepage, tutorial, reference, and long
  walkthrough pages look coherent.

## Out of scope

- Native Wolfram Documentation Center notebooks.
- Wolfram Cloud documentation deployment.
- Major rewrites of tutorial/reference prose.
- Changing the public Wolfram API.
- Adding GitHub Actions unless explicitly requested later.
- Versioned documentation.

## References

- Furo customization: https://furo.readthedocs.io/customisation/index.html
- Sphinx theme customization: https://www.sphinx-doc.org/en/master/tutorial/more-sphinx-customization.html
- MyST Parser: https://myst-parser.readthedocs.io/
