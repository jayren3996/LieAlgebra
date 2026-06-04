#!/usr/bin/env bash
# build-docs.sh -- assemble and build the Sphinx/Furo documentation site.
#
# Usage:
#   scripts/build-docs.sh          stage + generate + build into site/
#   scripts/build-docs.sh build    same as above
#   scripts/build-docs.sh serve    build, then serve site/ at localhost:8000
#   scripts/build-docs.sh clean    remove docs-site/ and site/
#
# One-time setup:
#   python3 -m pip install -r requirements-docs.txt

set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

WOLFRAM="/Applications/Wolfram.app/Contents/MacOS/wolframscript"
[ -x "$WOLFRAM" ] || WOLFRAM="wolframscript"

build_docs() {
  echo "==> staging Sphinx sources"
  python3 scripts/stage-sphinx-docs.py

  echo "==> generating Sphinx API reference"
  "$WOLFRAM" -file scripts/generate-reference.wls --sphinx

  echo "==> building Furo site"
  rm -rf site
  python3 -m sphinx.cmd.build -W -b dirhtml docs-site site
}

case "${1:-build}" in
  build)
    build_docs
    ;;
  serve)
    build_docs
    echo "==> serving http://localhost:${PORT:-8000}/"
    python3 -m http.server "${PORT:-8000}" -d site
    ;;
  clean)
    rm -rf docs-site site
    ;;
  *)
    echo "unknown command: $1 (use: build | serve | clean)" >&2
    exit 2
    ;;
esac
