#!/usr/bin/env bash
# build-site.sh -- assemble the documentation site's staging tree and build it.
#
# The source files live in three places and are composed into the git-ignored
# site-src/ directory that mkdocs.yml points at:
#   web/                       authored site pages (index, demos) + mathjax config
#   scripts/generate-reference.wls   the per-symbol API reference (generated)
#   docs/walkthrough.md        the guided tour, with image/links rewritten for the site
#   pics/*.png                 the walkthrough figures
#
# Usage:
#   scripts/build-site.sh           assemble + mkdocs build --strict (output in site/)
#   scripts/build-site.sh serve     assemble + mkdocs serve (live preview at :8000)
#   scripts/build-site.sh deploy    assemble + mkdocs gh-deploy (push to gh-pages)
#
# One-time setup:
#   pip install mkdocs-material
# After the first deploy, enable GitHub Pages in the repo settings with the
# source set to the gh-pages branch.

set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

# Prefer the app-bundled wolframscript; the /usr/local/bin one can be broken on macOS.
WOLFRAM="/Applications/Wolfram.app/Contents/MacOS/wolframscript"
[ -x "$WOLFRAM" ] || WOLFRAM="wolframscript"

GH_BLOB="https://github.com/jayren3996/LieAlgebra/blob/master"

echo "==> assembling site-src/"
rm -rf site-src
mkdir -p site-src
cp -R web/. site-src/                       # index.md, demos.md, javascripts/

echo "==> generating API reference"
"$WOLFRAM" -file scripts/generate-reference.wls   # writes site-src/reference/

echo "==> staging the guided tour and figures"
cp docs/walkthrough.md site-src/walkthrough.md
mkdir -p site-src/pics
cp pics/*.png site-src/pics/
# Rewrite the walkthrough's repo-relative links for the site. The MakeFigures
# link must be rewritten before the general ../pics/ rule (it matches both).
sed -i '' \
  -e "s#\.\./pics/MakeFigures\.wls#${GH_BLOB}/pics/MakeFigures.wls#g" \
  -e 's#\.\./pics/#pics/#g' \
  -e 's#\.\./README\.md#index.md#g' \
  -e 's#\.\./demos/#demos.md#g' \
  site-src/walkthrough.md

case "${1:-build}" in
  serve)  mkdocs serve ;;
  deploy) mkdocs gh-deploy --clean ;;
  build)  mkdocs build --strict ;;
  *) echo "unknown command: $1 (use: build | serve | deploy)" >&2; exit 2 ;;
esac
