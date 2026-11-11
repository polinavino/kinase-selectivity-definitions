#!/usr/bin/env bash
# Build the current Word manuscript and synchronized submission text/tables.
# For both PDFs too: python3 paper/build_submission.py
set -euo pipefail
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
python3 "$SCRIPT_DIR/build_submission.py" --word-only
