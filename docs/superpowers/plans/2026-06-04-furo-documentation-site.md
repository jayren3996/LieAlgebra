# Furo Documentation Site Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace the MkDocs Material web documentation build with a Sphinx + Furo build that preserves existing content and renders Wolfram examples as notebook-like input/output cells.

**Architecture:** Keep authored content in `web/`, `docs/walkthrough.md`, `pics/`, and `scripts/reference-details.wl`. Generate an ignored `docs-site/` Sphinx source tree through a Python staging script, then run `sphinx-build` into the ignored `site/` output directory. Extend the Wolfram reference generator with a `--sphinx` target so the old MkDocs pipeline can remain available while the new Furo site gets notebook-style examples.

**Tech Stack:** Wolfram Language/wolframscript, Python 3 standard library, Sphinx, Furo, MyST Parser, sphinx-copybutton, Pygments Mathematica lexer, MathJax.

---

### Task 1: Staging Transform Tests

**Files:**
- Create: `Tests/DocsStaging.py`
- Create: `scripts/stage-sphinx-docs.py`

- [ ] **Step 1: Write the failing tests**

Create `Tests/DocsStaging.py` with tests for the Markdown transformations:

```python
import unittest

from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import stage_sphinx_docs as stage


class DocsStagingTests(unittest.TestCase):
    def test_rewrites_mkdocs_admonition_to_myst_fence(self):
        source = '!!! warning "Large irreps are expensive"\n\n    Build small irreps.\n'
        expected = '```{warning} Large irreps are expensive\nBuild small irreps.\n```\n'
        self.assertEqual(stage.convert_markdown(source, "web/tutorials/representations.md"), expected)

    def test_rewrites_tabbed_installation_blocks(self):
        source = '=== "From a checkout"\n\n    ```mathematica\n    Needs["ClassicalLieAlgebra`"];\n    ```\n'
        rendered = stage.convert_markdown(source, "web/index.md")
        self.assertIn("#### From a checkout", rendered)
        self.assertIn('```mathematica\nNeeds["ClassicalLieAlgebra`"];\n```', rendered)

    def test_rewrites_material_image_attributes_to_plain_image(self):
        source = '![Highest weight](../pics/rep11-highest-weight.png){ .center width="120" }\n'
        self.assertEqual(
            stage.convert_markdown(source, "web/tutorials/young-tableaux.md"),
            '![Highest weight](../pics/rep11-highest-weight.png)\n',
        )

    def test_rewrites_markdown_links_for_sphinx_pages(self):
        source = "[Concepts](../concepts.md) and [`Irrep`](../reference/representations.md#irrep)"
        self.assertEqual(
            stage.convert_markdown(source, "web/tutorials/getting-started.md"),
            "[Concepts](../concepts) and [`Irrep`](../reference/representations#irrep)",
        )


if __name__ == "__main__":
    unittest.main()
```

- [ ] **Step 2: Run tests to verify they fail**

Run:

```bash
python3 -m unittest Tests.DocsStaging -v
```

Expected: FAIL or ERROR because `scripts/stage-sphinx-docs.py` does not exist yet.

- [ ] **Step 3: Create the minimal staging module**

Create `scripts/stage-sphinx-docs.py` with:

```python
#!/usr/bin/env python3
from pathlib import Path
import re


ROOT = Path(__file__).resolve().parents[1]
DOCS_SITE = ROOT / "docs-site"


def convert_markdown(text: str, source_name: str) -> str:
    text = convert_admonitions(text)
    text = convert_tabs(text)
    text = strip_material_blocks(text)
    text = strip_attribute_lists(text)
    text = rewrite_markdown_links(text)
    return text
```

Implement the helper functions until the tests pass.

- [ ] **Step 4: Run tests to verify they pass**

Run:

```bash
python3 -m unittest Tests.DocsStaging -v
```

Expected: all four tests pass.

### Task 2: Sphinx Staging Tree

**Files:**
- Modify: `scripts/stage-sphinx-docs.py`
- Modify: `.gitignore`
- Create: `web/sphinx/cla.css`
- Create: `web/sphinx/mathjax.js`

- [ ] **Step 1: Extend the staging script**

Add a `stage()` function that:

```python
def stage() -> None:
    if DOCS_SITE.exists():
        shutil.rmtree(DOCS_SITE)
    (DOCS_SITE / "_static").mkdir(parents=True)
    write_conf()
    copy_markdown_sources()
    copy_walkthrough()
    copy_pics()
    copy_static_assets()
```

It should write `docs-site/conf.py`, copy and convert `web/*.md`, copy and convert `web/tutorials/*.md`, copy `docs/walkthrough.md` to `docs-site/walkthrough.md` with figure links rewritten from `../pics/` to `pics/`, copy `pics/*.png` to `docs-site/pics/`, and copy Sphinx CSS/JS from `web/sphinx/`.

- [ ] **Step 2: Add Furo CSS and MathJax config**

Create `web/sphinx/cla.css` with Furo customizations for:

```css
.cla-example { ... }
.cla-example-row { ... }
.cla-example-label { ... }
.cla-example-input { ... }
.cla-example-output { ... }
.cla-home-hero { ... }
.cla-home-grid { ... }
```

Create `web/sphinx/mathjax.js` with MathJax inline and display math delimiters.

- [ ] **Step 3: Ignore generated Sphinx source**

Add `docs-site/` to `.gitignore`.

- [ ] **Step 4: Run staging**

Run:

```bash
python3 scripts/stage-sphinx-docs.py
```

Expected: `docs-site/conf.py`, converted Markdown pages, `_static/cla.css`, `_static/mathjax.js`, and `pics/*.png` exist.

### Task 3: Sphinx Reference Generator Target

**Files:**
- Modify: `scripts/generate-reference.wls`

- [ ] **Step 1: Write the failing smoke check**

Run the Sphinx staging script followed by the reference generator in Sphinx mode:

```bash
python3 scripts/stage-sphinx-docs.py
WOLFRAM="/Applications/Wolfram.app/Contents/MacOS/wolframscript"; [ -x "$WOLFRAM" ] || WOLFRAM="wolframscript"; "$WOLFRAM" -file scripts/generate-reference.wls --sphinx
rg -n 'cla-example|In\\[1\\]:=|Out\\[1\\]=' docs-site/reference/representations.md
```

Expected before implementation: the generator either ignores `--sphinx` or writes raw code/output blocks, so the `rg` check fails.

- [ ] **Step 2: Add target selection**

Update `scripts/generate-reference.wls` so:

```wolfram
target = If[MemberQ[$ScriptCommandLine, "--sphinx"], "Sphinx", "MkDocs"];
outRoot = If[target === "Sphinx", "docs-site", "site-src"];
outDir = FileNameJoin[{root, outRoot, "reference"}];
```

- [ ] **Step 3: Add notebook-style example rendering for Sphinx**

For the Sphinx target, render bare examples as raw HTML:

```html
<div class="cla-example">
  <div class="cla-example-row cla-example-in">
    <div class="cla-example-label">In[1]:=</div>
    <div class="cla-example-cell"><div class="highlight-mathematica notranslate"><pre>...</pre></div></div>
  </div>
  <div class="cla-example-row cla-example-out">
    <div class="cla-example-label">Out[1]=</div>
    <div class="cla-example-cell cla-example-output"><pre>...</pre></div>
  </div>
</div>
```

HTML-escape input and output. Keep MkDocs rendering unchanged for the default target.

- [ ] **Step 4: Run the smoke check**

Run:

```bash
python3 scripts/stage-sphinx-docs.py
WOLFRAM="/Applications/Wolfram.app/Contents/MacOS/wolframscript"; [ -x "$WOLFRAM" ] || WOLFRAM="wolframscript"; "$WOLFRAM" -file scripts/generate-reference.wls --sphinx
rg -n 'cla-example|In\\[1\\]:=|Out\\[1\\]=' docs-site/reference/representations.md
```

Expected: `rg` finds notebook-style example markup.

### Task 4: Build Script and Dependencies

**Files:**
- Create: `requirements-docs.txt`
- Create: `scripts/build-docs.sh`
- Modify: `README.md`

- [ ] **Step 1: Add docs dependencies**

Create `requirements-docs.txt`:

```text
sphinx>=7.0
furo>=2024.8.6
myst-parser>=2.0
sphinx-copybutton>=0.5
```

- [ ] **Step 2: Add Sphinx build script**

Create `scripts/build-docs.sh` that supports:

```bash
scripts/build-docs.sh          # stage + generate + sphinx build
scripts/build-docs.sh serve    # build + python http.server on :8000
scripts/build-docs.sh clean    # remove docs-site/ and site/
```

The build command should run:

```bash
python3 scripts/stage-sphinx-docs.py
"$WOLFRAM" -file scripts/generate-reference.wls --sphinx
sphinx-build -W -b html docs-site site
```

- [ ] **Step 3: Add README build instructions**

Add a short documentation-build section that tells contributors to run:

```bash
python3 -m pip install -r requirements-docs.txt
scripts/build-docs.sh
```

### Task 5: Build Verification and Polish

**Files:**
- Modify as needed: `scripts/stage-sphinx-docs.py`
- Modify as needed: `scripts/generate-reference.wls`
- Modify as needed: `web/sphinx/cla.css`

- [ ] **Step 1: Install docs dependencies if needed**

Run:

```bash
python3 -m pip install -r requirements-docs.txt
```

Expected: Sphinx, Furo, MyST Parser, and sphinx-copybutton import successfully.

- [ ] **Step 2: Build the site**

Run:

```bash
scripts/build-docs.sh
```

Expected: `sphinx-build` exits with code 0 and no warnings.

- [ ] **Step 3: Verify generated content**

Run:

```bash
rg -n 'cla-example|In\\[1\\]:=|Out\\[1\\]=' site/reference/representations/index.html
test -f site/index.html
test -f site/tutorials/getting-started/index.html
test -f site/walkthrough/index.html
test -f site/pics/su3-levels.png
```

Expected: all commands pass.

- [ ] **Step 4: Run package tests**

Run:

```bash
WOLFRAM="/Applications/Wolfram.app/Contents/MacOS/wolframscript"; [ -x "$WOLFRAM" ] || WOLFRAM="wolframscript"; "$WOLFRAM" -file scripts/runTests.wls
```

Expected: `ALL PASS`.

- [ ] **Step 5: Review the built site locally**

Serve the generated `site/` directory and check home, reference, tutorial, and walkthrough pages in a browser:

```bash
python3 -m http.server 8000 -d site
```

Expected: Furo layout loads, MathJax renders, figures load, and reference examples display as paired Wolfram notebook-style cells.
