#!/bin/bash
#- tutorials of the documentation: Quarto sources (docs/vignettes/*.qmd) rendered to Markdown with their figures (docs/tutorials)
#- the downstream tutorials need the full runs: render them where the results are (SCALPEL_RUN), then collect the outputs here
#- usage: build_tutorials.sh <directory with the rendered *.md & *_files> (default: docs/vignettes)
set -euo pipefail
DOCS=$(cd "$(dirname "$0")/.." && pwd)
SRC=${1:-$DOCS/vignettes}
OUT=$DOCS/tutorials

# Render (where the data are)
# +++++++++++++++++++++++++++
#   quarto render chromium.qmd --to gfm
#   SCALPEL_RUN=<GSE104556 run> quarto render chromium_downstream.qmd --to gfm
#   quarto render dropseq.qmd --to gfm
#   SCALPEL_RUN=<GSE268222 run> quarto render dropseq_downstream.qmd --to gfm

# Collect
# +++++++
for f in chromium chromium_downstream dropseq dropseq_downstream; do
    #- logos from docs/vignettes/img, author line of the vignette header removed
    sed -e 's#\./img/#../vignettes/img/#g' -e '/^author: Franz AKE/d' "$SRC/$f.md" > "$OUT/$f.md"
    if [ -d "$SRC/${f}_files/figure-commonmark" ]; then
        rm -rf "$OUT/${f}_files" && mkdir -p "$OUT/${f}_files"
        cp -r "$SRC/${f}_files/figure-commonmark" "$OUT/${f}_files/"
    fi
    echo "  $f"
done
