# Formalisation

Two Lean 4 files, compiled and kernel-checked against Mathlib:

- `GridExchange.lean`: soundness of disjoint-gate exchange certificates (the row–column rearrangement of the repeated circuit).
- `CapableState.lean`: the first-use step of the capable-state counterexample (Proposition 14, Appendix G, and the Comment).

| | |
|---|---|
| Lean toolchain | `leanprover/lean4:v4.34.1` (pinned in `lean-toolchain`) |
| Mathlib | tag `v4.34.1`, commit `d13f23b723b8a846827a245b89c10fc7d3f11612` (pinned in `lake-manifest.json`) |
| Axioms used | `GridExchange`: `propext` only. `CapableState`: Lean's standard `propext`, `Classical.choice`, `Quot.sound`. No `sorryAx`; checked with `#print axioms` on every theorem |
| Build record | `../verification/lean_check.txt` |

## Build

Install [elan](https://github.com/leanprover/elan); it selects the pinned toolchain automatically. Then:

```bash
cd lean
lake exe cache get   # download prebuilt Mathlib instead of compiling it
lake build
```

## `GridExchange.lean`

**Proved.**
- `exchange_sound`: in any monoid, possibly noncommutative, a certificate of adjacent exchanges of `Separate` gates preserves the ordered product. Gates are separate when they differ in both input wire and reservoir wire. The only hypothesis is that separate gates commute.
- `rectangle_two_by_two`, `rectangle_two_by_two_sound`: a concrete certificate taking the 2×2 row order to column order, and the resulting product equality.

**Not proved.** The general existence of a certificate for every rectangle is proved in prose and checked on finite examples in Python (`../code/return_checks.py`), but not in Lean. The tensor-product argument that justifies the commutation hypothesis physically is also outside this file.

## `CapableState.lean`

Conventions follow the paper: the system starts in `p = |0⟩⟨0|` (Bloch coordinate `1`), each collision is `U = cos η · I + i sin η · SWAP`, and `q = cos² η`. A diagonal qubit state with Bloch coordinate `z` is `diag((1 + z)/2, (1 - z)/2)`.

**Layer 1: scalar core.** The recursion `z ↦ q z + (1 - q) z_c` over a sequence of cells.
- `run_designated`: the designated initializer `τ^{⊗M}` gives first-output coordinate `q^M`.
- `run_product`: the product initializer `|1⟩⟨1| ⊗ τ^{⊗(M-1)}` gives `(2q - 1) q^(M-1)`.
- `product_within_iff`: for `q > 0`, that coordinate has absolute value at most `q^M` exactly when `1/3 ≤ q ≤ 1`.
- `zTwo_mem`, `run_adjusted`: for `0 < q < 1/3`, the second-cell coordinate `z₂ = q(1 - 2q)/(1 - q)` lies in `[0, 1]`, and the adjusted initializer gives coordinate `0` after two cells and for every later cell.

**Layer 2: physical update.** `ptrace2_pswap` and `ptrace2_partialSwap` use explicit `4 × 4` complex matrices for `U`, the Kronecker product, and a concretely defined partial trace over the cell. They prove that the system's reduced state after `U (ρ ⊗ σ) U†`, with `ρ` and `σ` diagonal, is the diagonal state with coordinate `q z_ρ + (1 - q) z_σ`.

Sanity lemmas pin down the encoding:
- `kron_apply_order`: the Kronecker index order is (system, cell).
- `ptrace2_kron`: `ptrace2 (A ⊗ B) = tr(B) · A`, so the partial trace removes the cell.
- `swapGate_conj_kron`: `SWAP (A ⊗ B) SWAP = B ⊗ A`.
- `pswap_unitary`: `U` is unitary.
- `diagState_one`, `diagState_neg_one`, `diagState_zero`: identify `p = |0⟩⟨0|`, `|1⟩⟨1|`, and `τ`.

**Layer 3: sequential exactness.**
- `redS_collide`: for any two-qubit gate and any system–environment state, however correlated, the system's reduced state after a collision with a fresh, uncorrelated cell depends only on its reduced state beforehand.
- `redS_circuit`: iterating this, the joint state of the system and all visited cells (`circuit`, in the appending model below) has a reduced system state that follows the scalar recursion. Every system–cell correlation is kept in the joint state.

**First use of each initializer.**
- `firstOutput_designated`, `firstOutput_product`, `firstOutput_adjusted`: the first output for each initializer. The adjusted initializer's output is exactly `τ = I/2`.
- `capable_first_use`: for every strict partial swap `0 < q < 1` and every `M ≥ 2`, there is a product initializer with first cell `|1⟩⟨1|` and all cells valid diagonal qubit states. Its first output has coordinate `z` with `|z| ≤ q^M`, and therefore satisfies `F ≥ 1 - e_M` with the closed-form fidelity below. The condition `M ≥ 2` is needed: for `q < 1/3`, the adjusted initializer requires two cells.

**Modelling choices and assumptions.**
- `circuit` appends each cell when the system reaches it. A cell the system has not yet met is untouched and uncorrelated with everything else, so this gives the same state as having all cells present from the start. That equivalence is the model, not a theorem in the file.
- Fidelity to `τ` is taken in the closed form `(1 + √(1 - z²))/2` for a diagonal qubit state (`fidTau`). Mathlib has no Uhlmann fidelity, so this formula is a definition here, not derived. `fidTau_ge_iff` shows that, with this form, `F ≥ 1 - e_M` holds exactly when `|z| ≤ q^M`, where `e_M = (1 - √(1 - q^{2M}))/2`.

**Not formalised.** The orthogonality of the product initializers to the attractor `p^{⊗M}` (immediate from the first cell being `|1⟩⟨1|`, but not stated in Lean), the near-designated initializer `(I - p^{⊗M})/(2^M - 1)` (not a product state), the reverse-direction comparison, the finite-block relaxation lemma (Lemma 12), the limit `S_{e_M}(n) → 0`, and the conclusion `e_M / S_{e_M}(n) → +∞`. Those steps, and the convergence argument in general, rest on the written proofs. No entropy, spectral, or resource theorem is formalised.

## Provenance

`GridExchange.lean` was first drafted in an environment that could not install Lean (`../verification/prior_lean_environment.txt`). Compiling it produced a single error, `failed to synthesize Decidable (Separate (0, 1) (1, 0))`, because `Separate` is a plain `def`. The only change was adding `instance : DecidableRel Separate`, which decides the existing definition. No definition, hypothesis, or theorem statement was altered.

Appendix K of the paper summarises these results.
