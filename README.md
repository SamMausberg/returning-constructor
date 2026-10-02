# Memory return in quantum machines

**Samuel Mausberg** (independent researcher) · samuelmausberg@gmail.com

This repository holds the manuscript sources, full proofs, computational checks, certified figure data, and Lean formalisation for *Memory return in quantum machines*. It also holds the accompanying Comment on PRL 128, 080401.

> We ask whether a quantum machine can repeatedly transform fresh inputs while returning its memory to its initial state. […] These results separate the return of a memory's reduced state from its ability to supply fresh independent randomness.

The full abstract is in [`pdf/Memory_return.pdf`](pdf/Memory_return.pdf).

## Contents

| Path | Contents |
|---|---|
| [`pdf/Memory_return.pdf`](pdf/Memory_return.pdf) | Full article, Physical Review A layout (REVTeX 4.2), 18 pp. |
| [`pdf/Comment_PRL.pdf`](pdf/Comment_PRL.pdf) | Comment on PRL **128**, 080401 (2022), 1 p. |
| [`pdf/Referee_packet.pdf`](pdf/Referee_packet.pdf) | The 15 numbered results with hypotheses and source locators, without proofs |
| [`source/`](source/) | LaTeX sources: `paper.tex`, `comment.tex`, theorem statements (`statements/`), full proofs (`appendices/`), TikZ diagrams and generated plots (`figures/`), and `references.bib` |
| [`code/`](code/) | Python: research checks (`return_checks.py`, `exact_homogenizer.py`, `extended_checks.py`), certified figure data (`certified_figures.py`), plots (`make_plots.py`), and the editorial audit (`editorial_audit.py`) |
| [`data/`](data/) | Certified figure data: rational interval endpoints (JSON) and decimal plotting coordinates (CSV) |
| [`lean/`](lean/) | Lean 4 + Mathlib proof that disjoint-gate exchange certificates are sound |
| [`verification/`](verification/) | Machine-readable records of the checks, figure certificates, the statement-preservation audits, and the Lean build |
| [`scripts/`](scripts/) | Build helpers that generate the statement-only packet and audit the theorem statements |

## Citation

If you use or build on this work, please cite it. A machine-readable [`CITATION.cff`](CITATION.cff) is included; GitHub shows it under "Cite this repository".

```bibtex
@unpublished{Mausberg2026MemoryReturn,
  author = {Mausberg, Samuel},
  title  = {Memory return in quantum machines},
  year   = {2026},
  note   = {Manuscript}
}
```

## Reproducing the results

Requires Python ≥ 3.11.

```bash
python3 -m venv .venv
.venv/bin/pip install -r requirements.txt
```

### Research checks

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .venv/bin/python code/return_checks.py
```

This rewrites `verification/return_checks.json`, which should be identical to `verification/prior_return_checks.json`. The suite uses seed 20261002. It includes exact finite algebraic checks and floating-point diagnostics. Its numerical results do not establish the asymptotic theorems. `code/exact_homogenizer.py` and `code/extended_checks.py` write `checks.json` and `extended_checks.json` next to themselves.

### Figures

```bash
.venv/bin/python code/certified_figures.py   # data/*.json, data/*.csv, verification/figure_certificates.json
.venv/bin/python code/make_plots.py          # source/figures/*.pdf and inspection PNGs
```

The figure data use exact rational arithmetic and outward-rounded interval enclosures, with no floating-point eigensolver. The plots show exact evaluations and rigorous enclosures, not fitted asymptotic laws.

### Lean formalisation

Requires [elan](https://github.com/leanprover/elan). The toolchain and Mathlib version are pinned.

```bash
cd lean
lake exe cache get
lake build
```

See [`lean/README.md`](lean/README.md) for exactly what is and is not formalised.

### Manuscripts

Requires TeX Live with REVTeX 4.2 (Debian/Ubuntu package `texlive-publishers`), AMS packages, TikZ, microtype, hyperref, cleveref, aliascnt, capt-of, BibTeX, and Python 3.

```bash
./build.sh
.venv/bin/python code/editorial_audit.py
```

`build.sh` writes intermediate files to `_build/` and the PDFs to `pdf/`. It regenerates `source/referee_packet.tex` from the compiled paper's labels and reruns the statement-preservation audit (`verification/claim_preservation.json`). `editorial_audit.py` confirms that the 15 statement files and the proof appendices match the hashes in `verification/baseline_hashes.json`, and it estimates the Comment's APS length. LaTeX output is not byte-reproducible. The hashes in `verification/publication_checks.json` refer to the committed PDFs.

### Code style

The Python is formatted and linted with [ruff](https://docs.astral.sh/ruff/), configured in `pyproject.toml`:

```bash
.venv/bin/pip install ruff
.venv/bin/ruff format code scripts && .venv/bin/ruff check code scripts
```

## Status

The manuscripts have not yet been peer reviewed. The acknowledgments and the verification record (Appendix K) describe the research assistance used. Appendix K also says which results rest on written proofs, exact finite checks, certified intervals, floating-point tests, or the Lean formalisation.

## License

Copyright © 2026 Samuel Mausberg. The code (`code/`, `scripts/`, `lean/`) is released under the MIT License, and the manuscripts and data (`source/`, `pdf/`, `data/`) under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/). See [`LICENSE`](LICENSE).
