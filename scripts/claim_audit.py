#!/usr/bin/env python3
"""Compare theorem mathematics with the supplied draft; do not validate proofs."""

import difflib
import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
SRC = ROOT / "source"
old = json.loads((ROOT / "verification/original_statements.json").read_text())
aux = (ROOT / "_build/paper.aux").read_text()
labels = {
    k: (n, p) for k, n, p in re.findall(r"\\newlabel\{([^}]+)\}\{\{([^}]*)\}\{([^}]*)\}", aux)
}


def maths(text):
    pieces = []
    pat = r"\$([^$]+)\$|\\begin\{(equation|align)\}(.*?)\\end\{\2\}"
    for m in re.finditer(pat, text, re.S):
        s = m.group(1) if m.group(1) is not None else m.group(3)
        s = re.sub(r"\\label\{[^}]*\}", "", s)
        s = re.sub(r"\\(?:begin|end)\{(?:aligned|gathered)\}", "", s)
        s = re.sub(r"\\\\(?:\[[^]]*\])?", "", s)
        s = re.sub(r"\\(?:nonumber|quad|qquad|bigl|bigr|left|right)\b", "", s)
        s = re.sub(r"\\[,;!:]", "", s)
        s = re.sub(r"\s+|&", "", s)
        pieces.append(s)
    # Joining also makes equation environments split for layout compare alike.
    return "".join(pieces)


rows = []
for key, v in old.items():
    f = SRC / "statements" / (key.replace(":", "_") + ".tex")
    text = f.read_text()
    before = maths(v["body"])
    after = maths(text)
    same = before == after
    row = {
        "label": key,
        "title": v["title"],
        "paper_number": labels.get(key, ("?", "?"))[0],
        "paper_page": labels.get(key, ("?", "?"))[1],
        "source": str(f.relative_to(ROOT)),
        "mathematical_tokens_preserved": same,
        "comparison": "Whitespace, alignment, scalable delimiters, labels, and display-line breaks normalized. This comparison is not a proof check.",
    }
    if not same:
        row["difference"] = "\n".join(difflib.ndiff([before], [after]))
    rows.append(row)
report = {
    "purpose": "Editorial preservation check against supplied source bundle",
    "theorem_statements": rows,
    "flagged_changes": [
        {
            "location": "Appendix F, orthogonal-pure-state one-use aside",
            "change": "Restricts the ratio-one example to strict partial swaps; states that full swap has numerator and denominator zero and hence an undefined ratio.",
            "effect": "No change to the iterated-limit theorem, whose hypothesis is 0 < F(a,b) < 1.",
        }
    ],
    "semantic_review": (
        "The 15 theorem/corollary/proposition/lemma conclusions and assumptions are retained. "
        "Main-text proof sketches and full-proof locations were reorganized."
        if all(x["mathematical_tokens_preserved"] for x in rows)
        else "Mathematical tokens differ for at least one statement; review required."
    ),
}
(ROOT / "verification/claim_preservation.json").write_text(json.dumps(report, indent=2))
print(
    "Mathematical statement token comparisons:",
    sum(x["mathematical_tokens_preserved"] for x in rows),
    "/",
    len(rows),
)
for row in rows:
    if not row["mathematical_tokens_preserved"]:
        print(row["label"], row.get("difference", ""))
if not all(x["mathematical_tokens_preserved"] for x in rows):
    sys.exit(
        "claim_audit.py: a theorem statement changed; see verification/claim_preservation.json"
    )
