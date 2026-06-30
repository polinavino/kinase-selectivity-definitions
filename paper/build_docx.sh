#!/usr/bin/env bash
# Regenerate the Molecular Informatics Word submission from the modular LaTeX
# source. Run from the paper/ directory after editing any sections/*.tex.
#
#   - numbered superscript citations + numbered reference list (nature.csl,
#     closest available match to the journal's Angewandte style; the journal
#     reformats references in production)
#   - run-in \paragraph headings (definitions, required properties) via the
#     Lua filter, matching the PDF
#   - "References" heading on the bibliography
set -euo pipefail

NATURE="$(kpsewhich nature.csl)"

# references.bib uses achemso's "and {et~al.}" form so the LaTeX/PDF renders
# "et al." (achemso prints a literal "others" otherwise). citeproc wants the
# standard "and others" marker, so derive a citeproc-form bib on the fly.
CSL_BIB="${TMPDIR:-/tmp}/refs_csl_$$.bib"
trap 'rm -f "$CSL_BIB"' EXIT
sed 's/and {et~al\.}}/and others}/' references.bib > "$CSL_BIB"

pandoc main.tex \
  --resource-path=.:sections \
  --citeproc \
  --bibliography="$CSL_BIB" \
  --csl="$NATURE" \
  --metadata reference-section-title=References \
  --lua-filter=runin-headings.lua \
  -o molecular_informatics_submission.docx

# pandoc emits no page numbers; add a centered page-number footer
python3 add_page_numbers.py molecular_informatics_submission.docx

echo "Wrote molecular_informatics_submission.docx"
