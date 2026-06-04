import sys
import unittest
from pathlib import Path


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

    def test_sphinx_mathjax_config_is_not_overridden_by_custom_static_js(self):
        conf = stage.sphinx_conf()
        self.assertIn('"sphinx.ext.mathjax"', conf)
        self.assertIn("mathjax3_config", conf)
        self.assertNotIn("html_js_files", conf)
        self.assertNotIn("mathjax.js", stage.SPHINX_STATIC_ASSETS)


if __name__ == "__main__":
    unittest.main()
