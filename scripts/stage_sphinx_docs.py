#!/usr/bin/env python3
"""Stage Sphinx/Furo documentation sources for ClassicalLieAlgebra."""

from __future__ import annotations

import re
import shutil
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DOCS_SITE = ROOT / "docs-site"
WEB = ROOT / "web"
GH_BLOB = "https://github.com/jayren3996/LieAlgebra/blob/master"
SPHINX_STATIC_ASSETS = ("cla.css",)


def convert_admonitions(text: str) -> str:
    lines = text.splitlines()
    out: list[str] = []
    i = 0
    while i < len(lines):
        match = re.match(r'^!!!\s+([A-Za-z0-9_-]+)(?:\s+"([^"]+)")?\s*$', lines[i])
        if not match:
            out.append(lines[i])
            i += 1
            continue

        kind, title = match.group(1), match.group(2) or ""
        header = f"```{{{kind}}}" + (f" {title}" if title else "")
        out.append(header)
        i += 1

        if i < len(lines) and lines[i] == "":
            i += 1

        while i < len(lines):
            line = lines[i]
            if line == "":
                out.append("")
                i += 1
                continue
            if line.startswith("    "):
                out.append(line[4:])
                i += 1
                continue
            break

        while out and out[-1] == "":
            out.pop()
        out.append("```")

    return "\n".join(out) + ("\n" if text.endswith("\n") else "")


def convert_tabs(text: str) -> str:
    lines = text.splitlines()
    out: list[str] = []
    i = 0
    while i < len(lines):
        match = re.match(r'^===\s+"([^"]+)"\s*$', lines[i])
        if not match:
            out.append(lines[i])
            i += 1
            continue

        out.append(f"#### {match.group(1)}")
        i += 1
        if i < len(lines) and lines[i] == "":
            i += 1

        while i < len(lines):
            line = lines[i]
            if re.match(r'^===\s+"', line):
                break
            if line.startswith("    "):
                out.append(line[4:])
            else:
                out.append(line)
            i += 1

    return "\n".join(out) + ("\n" if text.endswith("\n") else "")


def strip_material_blocks(text: str) -> str:
    text = re.sub(r'^<div[^>\n]*\bmarkdown[^>\n]*>\n?', "", text, flags=re.MULTILINE)
    text = re.sub(r'^</div>\n?', "", text, flags=re.MULTILINE)
    text = re.sub(r':(?:material|octicons)-[A-Za-z0-9_-]+:(?:\{[^}\n]*\})?\s*', "", text)
    return text


def strip_attribute_lists(text: str) -> str:
    text = re.sub(r'(\[[^\]\n]*\]\([^)]+\))\{[^}\n]*\}', r"\1", text)
    return text


def rewrite_markdown_links(text: str) -> str:
    def repl(match: re.Match[str]) -> str:
        target, anchor = match.group(1), match.group(2) or ""
        return f"({target}{anchor})"

    return re.sub(r'\((?!https?:|mailto:)([^)\s#]+)\.md(#[^)]+)?\)', repl, text)


def convert_markdown(text: str, source_name: str) -> str:
    text = convert_admonitions(text)
    text = convert_tabs(text)
    text = strip_material_blocks(text)
    text = strip_attribute_lists(text)
    text = rewrite_markdown_links(text)
    return text


def sphinx_home() -> str:
    return """# ClassicalLieAlgebra

```{raw} html
<section class="cla-home-hero">
  <p class="cla-kicker">Wolfram Language paclet</p>
  <p class="cla-lede">Exact generators, bases, Young tableaux, and irreducible representations of the classical Lie algebras: su(n), so(n), and sp(2n).</p>
  <p class="cla-actions">
    <a class="cla-button cla-button-primary" href="walkthrough/">Guided tour of SU(3)</a>
    <a class="cla-button" href="reference/">API reference</a>
  </p>
</section>
```

![The (1,1) octet of SU(3), pictured as a three-level system](pics/su3-levels.png)

`ClassicalLieAlgebra` is a Wolfram Language paclet for working with the classical
simple Lie algebras: the special unitary `su(n)`, special orthogonal `so(n)`, and
symplectic `sp(2n)` families. It gives you their generators in several standard
bases, a Young-tableau toolkit for many-body wavefunctions, and a representation
engine that builds any irreducible representation from its highest weight.

## Why ClassicalLieAlgebra

```{raw} html
<div class="cla-home-grid">
  <section>
    <h3>Exact, all the way down</h3>
    <p>Every root, weight, and matrix element is computed over exact arithmetic. Results stay inspectable and reproducible.</p>
    <a href="validation/">How it is validated</a>
  </section>
  <section>
    <h3>All four families, one interface</h3>
    <p>Use SU[n], SO[n], Sp[n], or the canonical LieAlgebra[type, rank] across the A, B, C, and D families.</p>
    <a href="concepts/">Read the concepts</a>
  </section>
  <section>
    <h3>Reaches the spinors</h3>
    <p>Build explicit generator matrices for every classical irrep, including orthogonal spinors outside tensor powers of the defining representation.</p>
    <a href="tutorials/spinors/">Spinors of so(N)</a>
  </section>
</div>
```

## Installation

#### From a checkout

```mathematica
PacletDirectoryLoad["/path/to/LieAlgebra"];
Needs["ClassicalLieAlgebra`"];
```

#### From a release

Once a tagged release is published, install the built paclet from its release
asset:

```mathematica
PacletInstall["https://github.com/jayren3996/LieAlgebra/releases/download/vX.Y.Z/ClassicalLieAlgebra-X.Y.Z.paclet"];
Needs["ClassicalLieAlgebra`"];
```

## Quick start

```mathematica
Needs["ClassicalLieAlgebra`"];

Generators[SU[3]]
Generators[SU[3], "Chevalley"]

ir = Irrep[SU[3], {1, 1}];
RepresentationDimension[ir]
CasimirEigenvalue[ir]
WeightSystem[ir]
RepresentationMatrices[ir]

RepresentationDimension[Irrep[SO[5], {0, 1}]]
```

## What it covers

| | $A_n=\\mathfrak{su}(n{+}1)$ | $B_n=\\mathfrak{so}(2n{+}1)$ | $C_n=\\mathfrak{sp}(2n)$ | $D_n=\\mathfrak{so}(2n)$ |
| :--- | :---: | :---: | :---: | :---: |
| Generators, Cartan-Weyl, Chevalley | yes | yes | yes | yes |
| Irreps: dimension, weights, Casimir | yes | yes | yes | yes |
| Irreps: explicit generator matrices | yes | yes | yes | yes |
| Spinor representations | - | yes | - | yes |

```{toctree}
:hidden:
:maxdepth: 2

concepts
walkthrough
tutorials/index
tutorials/getting-started
tutorials/representations
tutorials/spinors
tutorials/young-tableaux
tutorials/physics-applications
reference/index
glossary
comparison
validation
changelog
cite
```
"""


def sphinx_conf() -> str:
    return '''project = "ClassicalLieAlgebra"
author = "Jie Ren"
copyright = "2026, Jie Ren"

extensions = [
    "myst_parser",
    "sphinx.ext.mathjax",
    "sphinx_copybutton",
]

source_suffix = {".md": "markdown"}
master_doc = "index"
exclude_patterns = []

html_theme = "furo"
html_title = "ClassicalLieAlgebra"
html_static_path = ["_static"]
html_css_files = ["cla.css"]
html_theme_options = {
    "light_css_variables": {
        "color-brand-primary": "#b23a2e",
        "color-brand-content": "#9f3128",
    },
    "dark_css_variables": {
        "color-brand-primary": "#f08a7f",
        "color-brand-content": "#f08a7f",
    },
}

myst_enable_extensions = ["amsmath", "dollarmath"]
myst_heading_anchors = 3
myst_links_external_new_tab = True
suppress_warnings = ["myst.header", "myst.xref_missing", "misc.highlighting_failure"]

mathjax3_config = {
    "tex": {
        "inlineMath": [["$", "$"], ["\\\\(", "\\\\)"]],
        "displayMath": [["$$", "$$"], ["\\\\[", "\\\\]"]],
    }
}

copybutton_selector = "div.highlight pre, .cla-example-input pre"
copybutton_prompt_text = r"^In\\[[0-9]+\\]:=\\s*"
copybutton_prompt_is_regexp = True
'''


def write_conf() -> None:
    conf = sphinx_conf()
    (DOCS_SITE / "conf.py").write_text(conf, encoding="utf-8")


def copy_markdown_sources() -> None:
    for src in sorted(WEB.glob("*.md")):
        if src.name == "index.md":
            (DOCS_SITE / "index.md").write_text(sphinx_home(), encoding="utf-8")
            continue
        dest = DOCS_SITE / src.name
        dest.write_text(convert_markdown(src.read_text(encoding="utf-8"), str(src.relative_to(ROOT))), encoding="utf-8")

    tutorials_out = DOCS_SITE / "tutorials"
    tutorials_out.mkdir(parents=True, exist_ok=True)
    for src in sorted((WEB / "tutorials").glob("*.md")):
        dest = tutorials_out / src.name
        dest.write_text(convert_markdown(src.read_text(encoding="utf-8"), str(src.relative_to(ROOT))), encoding="utf-8")


def copy_walkthrough() -> None:
    text = (ROOT / "docs" / "walkthrough.md").read_text(encoding="utf-8")
    text = text.replace("../pics/MakeFigures.wls", f"{GH_BLOB}/pics/MakeFigures.wls")
    text = text.replace("../pics/", "pics/")
    text = text.replace("../README.md", "index")
    text = text.replace("../demos/", "tutorials/index")
    text = convert_markdown(text, "docs/walkthrough.md")
    (DOCS_SITE / "walkthrough.md").write_text(text, encoding="utf-8")


def copy_pics() -> None:
    pic_out = DOCS_SITE / "pics"
    pic_out.mkdir(parents=True, exist_ok=True)
    for src in sorted((ROOT / "pics").glob("*.png")):
        shutil.copy2(src, pic_out / src.name)


def copy_static_assets() -> None:
    static_out = DOCS_SITE / "_static"
    for name in SPHINX_STATIC_ASSETS:
        shutil.copy2(WEB / "sphinx" / name, static_out / name)


def stage() -> None:
    if DOCS_SITE.exists():
        shutil.rmtree(DOCS_SITE)
    (DOCS_SITE / "_static").mkdir(parents=True)
    write_conf()
    copy_markdown_sources()
    copy_walkthrough()
    copy_pics()
    copy_static_assets()


def main() -> None:
    stage()
    print(f"staged Sphinx docs in {DOCS_SITE}")


if __name__ == "__main__":
    main()
