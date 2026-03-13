#!/bin/bash
# ==============================================================================
# Render Manuscript and Supplementary Material to PDF
# ==============================================================================
# Uses pandoc with pdflatex to convert markdown manuscripts to publication-
# quality PDFs with embedded figures and formatted tables.
#
# Requirements:
#   - pandoc (>= 2.0)
#   - pdflatex (TeX Live or similar)
#
# Usage:
#   cd manuscript/
#   bash render_manuscript.sh
#
# Outputs:
#   manuscript/manuscript_main.pdf
#   manuscript/supplementary_material.pdf
# ==============================================================================

set -e

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$SCRIPT_DIR"

echo "============================================================"
echo "Rendering Manuscript and Supplementary Material"
echo "============================================================"
echo ""

# Check dependencies
if ! command -v pandoc &> /dev/null; then
    echo "ERROR: pandoc is not installed."
    echo "  Install with: brew install pandoc"
    exit 1
fi

if ! command -v pdflatex &> /dev/null; then
    echo "ERROR: pdflatex is not installed."
    echo "  Install with: brew install --cask mactex"
    exit 1
fi

echo "pandoc version: $(pandoc --version | head -1)"
echo "pdflatex: $(which pdflatex)"
echo ""

# ---- Render main manuscript ----
echo "--- Rendering main manuscript ---"
pandoc manuscript_main.md \
    -o manuscript_main.pdf \
    --pdf-engine=pdflatex \
    --variable=colorlinks:true \
    --variable=linkcolor:blue \
    --variable=urlcolor:blue \
    --variable=toccolor:blue \
    2>&1

if [ $? -eq 0 ]; then
    echo "  -> Saved: manuscript_main.pdf"
else
    echo "  ERROR: Failed to render main manuscript."
    exit 1
fi

echo ""

# ---- Render supplementary material ----
echo "--- Rendering supplementary material ---"
pandoc supplementary_material.md \
    -o supplementary_material.pdf \
    --pdf-engine=pdflatex \
    --variable=colorlinks:true \
    --variable=linkcolor:blue \
    --variable=urlcolor:blue \
    --variable=toccolor:blue \
    2>&1

if [ $? -eq 0 ]; then
    echo "  -> Saved: supplementary_material.pdf"
else
    echo "  ERROR: Failed to render supplementary material."
    exit 1
fi

echo ""
echo "============================================================"
echo "Rendering complete"
echo "============================================================"
echo ""
echo "Outputs:"
echo "  $SCRIPT_DIR/manuscript_main.pdf"
echo "  $SCRIPT_DIR/supplementary_material.pdf"
echo ""
echo "To view:"
echo "  open manuscript_main.pdf"
echo "  open supplementary_material.pdf"
