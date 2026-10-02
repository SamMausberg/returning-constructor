#!/usr/bin/env python3
"""Audit preservation and estimate Comment length. It does not validate proofs."""

import hashlib
import json
import re
from pathlib import Path

import fitz

R = Path(__file__).resolve().parents[1]
V = R / "verification"
V.mkdir(exist_ok=True)


def digest(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


# Compare against the committed baseline manifest; it is never regenerated here.
old = json.loads((V / "baseline_hashes.json").read_text())
statements = []
for p in sorted((R / "source/statements").glob("*.tex")):
    k = str(p.relative_to(R))
    statements.append({"file": k, "sha256": digest(p), "unchanged": digest(p) == old[k]})
proofs = []
for p in sorted((R / "source/appendices").glob("*.tex")):
    k = str(p.relative_to(R))
    # The verification record and source audit are prose; the proof appendices stay byte-identical.
    if p.name not in ("verification.tex", "audit.tex"):
        proofs.append({"file": k, "sha256": digest(p), "unchanged": digest(p) == old[k]})
assert all(x["unchanged"] for x in statements + proofs)
aux = (R / "_build/paper.aux").read_text()
labels = {
    "cor:exactreturn": ("corollary", "2"),
    "prop:spectralreturn": ("proposition", "5"),
    "cor:singlecost": ("corollary", "7"),
}
for key, (kind, num) in labels.items():
    line = next(l for l in aux.splitlines() if l.startswith("\\newlabel{" + key + "@cref}"))
    assert "[" + kind + "]" in line, (key, line)
    text = "\n".join(p.get_text() for p in fitz.open(R / "_build/paper.pdf"))
    assert "Proof of " + kind.capitalize() + " " + num + "." in text
(V / "editorial_preservation.json").write_text(
    json.dumps(
        {
            "statements": statements,
            "proof_and_audit_appendices": proofs,
            "corrected_reference_types": labels,
            "meaning": "Byte comparison checks preservation, not mathematical correctness.",
        },
        indent=2,
    )
)
# APS word-equivalent estimate. Inline formulae are counted as one text token;
# a second deliberately generous count uses the formulae's TeX token count.
s = (R / "source/comment.tex").read_text().split("\\maketitle", 1)[1]
s = s.split("\\begin{thebibliography}", 1)[0]
s = re.sub(r"\\begin\{acknowledgments\}.*?\\end\{acknowledgments\}", "", s, flags=re.S)
rows = 0
for m in list(re.finditer(r"\\begin\{(equation|align)\}(.*?)\\end\{\1\}", s, re.S)):
    rows += 1 + len(re.findall(r"\\\\", m.group(2)))
s = re.sub(r"\\begin\{(equation|align)\}.*?\\end\{\1\}", "", s, flags=re.S)
maths = re.findall(r"\$([^$]*)\$", s)
s = re.sub(r"\$[^$]*\$", " INLINEFORMULA ", s)
s = re.sub(r"\\(?:includegraphics)(?:\[[^]]*\])?\{[^}]*\}", "", s)
s = re.sub(r"\\(?:label|cite|ref|eqref)\{[^}]*\}", "", s)
s = re.sub(r"\\(?:begin|end)\{[^}]*\}(?:\[[^]]*\])?", "", s)
s = re.sub(r"\\[A-Za-z]+\*?(?:\[[^]]*\])?", " ", s)
words = re.findall(r"[A-Za-z0-9]+(?:['-][A-Za-z0-9]+)*", s)
body = len(words)
math_tokens = sum(max(1, len(re.findall(r"\\[A-Za-z]+|[A-Za-z0-9]+", m))) for m in maths)
f = fitz.open(R / "source/figures/comment_ratio.pdf")
rect = f[0].rect
aspect = rect.width / rect.height
figure = 150 / aspect + 20
estimate = body + 16 * rows + figure
conservative = body - len(maths) + math_tokens + 16 * rows + figure
report = {
    "policy_verified": "2026-10-02",
    "limit_word_equivalents": 750,
    "policy_urls": [
        "https://journals.aps.org/prl/authors",
        "https://journals.aps.org/authors/length-guide",
    ],
    "included": "Body, caption, data availability; displayed mathematics and figure area by APS weights.",
    "excluded": "Title, author/date, bibliography, acknowledgments.",
    "body_caption_and_inline_math_tokens": body,
    "inline_math_expressions": len(maths),
    "inline_math_TeX_tokens": math_tokens,
    "display_rows": rows,
    "display_word_equivalents": 16 * rows,
    "figure_aspect_ratio": aspect,
    "figure_word_equivalents": figure,
    "estimated_total": estimate,
    "generous_inline_math_total": conservative,
    "interpretation": "An author estimate, not the journal production count. The generous variant charges each inline TeX command or alphanumeric token as a word.",
}
assert conservative < 750, report
(V / "comment_length.json").write_text(json.dumps(report, indent=2))
print("Preserved statement files:", len(statements), "/", len(statements))
print("Preserved proof/audit appendix files:", len(proofs), "/", len(proofs))
print(
    "Comment estimate:",
    round(estimate, 1),
    "word-equivalents; generous inline-math estimate:",
    round(conservative, 1),
)
