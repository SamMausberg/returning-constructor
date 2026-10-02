#!/usr/bin/env bash
set -euo pipefail
ROOT=$(cd "$(dirname "$0")" && pwd)
PYTHON=${PYTHON:-python3}
mkdir -p "$ROOT/_build" "$ROOT/pdf"
if command -v bibtex >/dev/null 2>&1; then BIBTEX=$(command -v bibtex)
elif command -v bibtex.original >/dev/null 2>&1; then BIBTEX=$(command -v bibtex.original)
else echo "BibTeX is required." >&2; exit 1; fi
# build <name> [source dir]: three pdflatex passes; BibTeX only where references.bib is used.
build() {
 local name=$1 dir=${2:-$ROOT/source}
 (cd "$dir" && pdflatex -interaction=nonstopmode -halt-on-error -output-directory="$ROOT/_build" "$name.tex" >"$ROOT/_build/${name}_pass1.txt")
 if [[ "$name" == paper || "$name" == referee_packet ]]; then
  (cd "$ROOT/_build" && BIBINPUTS="$ROOT/source:" "$BIBTEX" "$name" >"${name}_bibtex.txt")
 fi
 for pass in 2 3; do
  (cd "$dir" && pdflatex -interaction=nonstopmode -halt-on-error -output-directory="$ROOT/_build" "$name.tex" >"$ROOT/_build/${name}_pass${pass}.txt")
 done
}
build paper
"$PYTHON" "$ROOT/scripts/make_packet.py"
build comment
build referee_packet
cp "$ROOT/_build/paper.pdf" "$ROOT/pdf/Memory_return.pdf"
cp "$ROOT/_build/comment.pdf" "$ROOT/pdf/Comment_PRL.pdf"
cp "$ROOT/_build/referee_packet.pdf" "$ROOT/pdf/Referee_packet.pdf"
"$PYTHON" "$ROOT/scripts/claim_audit.py"
# Private correspondence (git-ignored); built only where it exists locally.
if [[ -f "$ROOT/private/briefing.tex" ]]; then
 build briefing "$ROOT/private"; cp "$ROOT/_build/briefing.pdf" "$ROOT/private/Briefing.pdf"
fi
if [[ -f "$ROOT/private/email.tex" ]]; then
 build email "$ROOT/private"; cp "$ROOT/_build/email.pdf" "$ROOT/private/Cover_email.pdf"
fi
if [[ -f "$ROOT/private/package_email.py" ]]; then "$PYTHON" "$ROOT/private/package_email.py"; fi
printf '\nBuilt paper, Comment, and statement-only packet in pdf/.\n'
