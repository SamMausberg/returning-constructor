# Formalisation: disjoint-gate exchange certificates

`GridExchange.lean` compiles and is kernel-checked by Lean 4 against Mathlib.

| | |
|---|---|
| Lean toolchain | `leanprover/lean4:v4.34.1` (pinned in `lean-toolchain`) |
| Mathlib | tag `v4.34.1`, commit `d13f23b723b8a846827a245b89c10fc7d3f11612` (pinned in `lake-manifest.json`) |
| Axioms used | `propext` only; no `sorryAx` (checked with `#print axioms`) |
| Build record | `../verification/lean_check.txt` |

## Build

Install [elan](https://github.com/leanprover/elan); it selects the pinned toolchain automatically. Then:

```bash
cd lean
lake exe cache get   # download prebuilt Mathlib instead of compiling it
lake build
```

## What is proved

- `exchange_sound`: in any monoid, possibly noncommutative, a certificate of adjacent exchanges of `Separate` gates preserves the ordered product. Gates are separate when they differ in both input wire and reservoir wire. The only hypothesis is that separate gates commute.
- `rectangle_two_by_two`, `rectangle_two_by_two_sound`: a concrete certificate taking the 2×2 row order to column order, and the resulting product equality.

## What is not proved here

The general existence of a certificate for every rectangle is proved in prose and checked on finite examples in Python (`../code/return_checks.py`), but not in Lean. The tensor-product argument that justifies the commutation hypothesis physically is also outside this file. The quantum-channel, spectral, entropy, and resource theorems are not formalised.

## Provenance

The file was first drafted in an environment that could not install Lean (`../verification/prior_lean_environment.txt`). Compiling it produced a single error, `failed to synthesize Decidable (Separate (0, 1) (1, 0))`, because `Separate` is a plain `def`. The only change was adding `instance : DecidableRel Separate`, which decides the existing definition. No definition, hypothesis, or theorem statement was altered.

Appendix K of the paper summarises this result.
