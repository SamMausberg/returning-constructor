import Mathlib

/-!
# The first-use step of the capable-state counterexample

Formal evidence for the first-use step of Proposition 14 of
"Memory return in quantum machines" (Appendix G, and the Comment).

The forward task sends a system qubit in `p = |0⟩⟨0|` (Bloch coordinate `1`)
through reservoir cells, one collision per cell. Each collision is the partial swap
`U = cos η · I + i sin η · SWAP`, and `q = cos² η`.

* Layer 1 (scalar core): the Bloch recursion for the designated initializer `τ^{⊗M}`,
  the product initializer `|1⟩⟨1| ⊗ τ^{⊗(M-1)}`, and the adjusted initializer used
  when `q < 1/3`.
* Layer 2 (physical update): with explicit `4 × 4` complex matrices and a concrete
  partial trace, one collision maps a diagonal system state with coordinate `z` and a
  diagonal cell with coordinate `zc` to the diagonal state with coordinate
  `q z + (1 - q) zc`.
* Layer 3 (sequential exactness): the reduced update stays exact when the system is
  already correlated with earlier cells, provided each new cell starts uncorrelated.
  Iterating gives the exact first output in a circuit model where each cell is appended
  when the system reaches it. That this equals the circuit with all cells present from
  the start (unvisited cells are untouched and uncorrelated) is a modelling step, not a
  theorem in this file.

The relaxation lemma and the limit `S_e(n) → 0` (later uses) are not formalised here.
-/

open Matrix Complex
open scoped Kronecker

namespace CapableState

/-! ## Layer 1: the scalar core -/

/-- Reduced Bloch update of the system after one collision with a fresh
diagonal cell of Bloch coordinate `zc`. -/
def step (q z zc : ℝ) : ℝ := q * z + (1 - q) * zc

/-- System Bloch coordinate after meeting cells `zc 0, …, zc (k - 1)` in order,
starting from coordinate `z`. -/
def run (q : ℝ) (zc : ℕ → ℝ) (z : ℝ) : ℕ → ℝ
  | 0 => z
  | k + 1 => step q (run q zc z k) (zc k)

/-- The designated initializer `τ^{⊗M}`: every cell has Bloch coordinate `0`. -/
def designatedCells : ℕ → ℝ := fun _ => 0

/-- The product initializer `|1⟩⟨1| ⊗ τ^{⊗(M-1)}`. -/
def productCells : ℕ → ℝ := fun k => if k = 0 then -1 else 0

/-- Second-cell coordinate used when `q < 1/3`. -/
noncomputable def zTwo (q : ℝ) : ℝ := q * (1 - 2 * q) / (1 - q)

/-- The adjusted product initializer for `q < 1/3`:
`|1⟩⟨1| ⊗ σ(z₂) ⊗ τ^{⊗(M-2)}`. -/
noncomputable def adjustedCells (q : ℝ) : ℕ → ℝ := fun k =>
  if k = 0 then -1 else if k = 1 then zTwo q else 0

/-- The designated initializer gives first-output coordinate `q^M`. -/
theorem run_designated (q : ℝ) (M : ℕ) : run q designatedCells 1 M = q ^ M := by
  induction M with
  | zero => simp [run]
  | succ k ih => simp [run, step, designatedCells, ih, pow_succ, mul_comm]

/-- The product initializer gives first-output coordinate `(2q - 1) q^(M-1)`,
stated with `M = k + 1`. -/
theorem run_product (q : ℝ) (k : ℕ) :
    run q productCells 1 (k + 1) = (2 * q - 1) * q ^ k := by
  induction k with
  | zero => simp [run, step, productCells]; ring
  | succ k ih =>
      rw [run, ih]
      simp [step, productCells, pow_succ]
      ring

/-- For `q > 0`, the product initializer's first-output coordinate has absolute value
at most `q^M` exactly when `1/3 ≤ q ≤ 1` (stated with `M = k + 1`). -/
theorem product_within_iff {q : ℝ} (hq : 0 < q) (k : ℕ) :
    |(2 * q - 1) * q ^ k| ≤ q ^ (k + 1) ↔ 1 / 3 ≤ q ∧ q ≤ 1 := by
  have hk : 0 < q ^ k := pow_pos hq k
  rw [abs_mul, abs_of_pos hk, pow_succ, mul_comm (q ^ k) q,
    mul_le_mul_iff_of_pos_right hk, abs_le]
  constructor
  · rintro ⟨h1, h2⟩
    constructor <;> linarith
  · rintro ⟨h1, h2⟩
    constructor <;> linarith

/-- For `0 < q < 1/3`, the second-cell coordinate is a valid Bloch coordinate. -/
theorem zTwo_mem {q : ℝ} (hq0 : 0 < q) (hq : q < 1 / 3) : zTwo q ∈ Set.Icc (0 : ℝ) 1 := by
  have h1q : 0 < 1 - q := by linarith
  constructor
  · unfold zTwo
    apply div_nonneg
    · apply mul_nonneg hq0.le; linarith
    · exact h1q.le
  · unfold zTwo
    rw [div_le_one h1q]
    nlinarith [sq_nonneg (1 - q), sq_nonneg q]

/-- For `q < 1`, the adjusted initializer gives coordinate `0` after two cells
and keeps it at `0` for every later cell (stated with `M = k + 2`). -/
theorem run_adjusted {q : ℝ} (hq1 : q < 1) (k : ℕ) :
    run q (adjustedCells q) 1 (k + 2) = 0 := by
  have h1q : (1 - q) ≠ 0 := by linarith
  induction k with
  | zero =>
      simp only [run, step, adjustedCells, zTwo]
      simp
      field_simp
      ring
  | succ k ih =>
      rw [show k + 1 + 2 = (k + 2) + 1 by ring, run, ih]
      simp [step, adjustedCells]

/-! ### Link to the capability set `V_{e_M}`

For a diagonal qubit state with Bloch coordinate `z`, fidelity to `τ = I/2` has the
closed form `(1 + √(1 - z²)) / 2`. That closed form is taken as the definition
`fidTau` here; it is not derived from the Uhlmann fidelity, which Mathlib lacks. -/

/-- Fidelity of the diagonal qubit state with Bloch coordinate `z` to `τ = I/2`
(closed form, see above). -/
noncomputable def fidTau (z : ℝ) : ℝ := (1 + Real.sqrt (1 - z ^ 2)) / 2

/-- The first-use tolerance `e_M = (1 - √(1 - q^{2M})) / 2`. -/
noncomputable def eM (q : ℝ) (M : ℕ) : ℝ := (1 - Real.sqrt (1 - q ^ (2 * M))) / 2

/-- With the closed form above, the first output lies in `V_{e_M}` exactly when its
Bloch coordinate has absolute value at most `q^M`. -/
theorem fidTau_ge_iff {q z : ℝ} (hq0 : 0 ≤ q) (hz : |z| ≤ 1) (M : ℕ) :
    1 - eM q M ≤ fidTau z ↔ |z| ≤ q ^ M := by
  have hz2 : 0 ≤ 1 - z ^ 2 := by
    have : z ^ 2 ≤ 1 := by
      rw [← sq_abs]
      nlinarith [abs_nonneg z]
    linarith
  have hqM : 0 ≤ q ^ M := pow_nonneg hq0 M
  have key : 1 - eM q M ≤ fidTau z ↔ Real.sqrt (1 - q ^ (2 * M)) ≤ Real.sqrt (1 - z ^ 2) := by
    unfold eM fidTau
    constructor <;> intro h <;> linarith
  have hpow : q ^ (2 * M) = (q ^ M) ^ 2 := by rw [pow_mul']
  rw [key, Real.sqrt_le_sqrt_iff hz2, hpow]
  constructor
  · intro h
    have h2 : z ^ 2 ≤ (q ^ M) ^ 2 := by linarith
    have h3 := sq_le_sq.mp h2
    rwa [abs_of_nonneg hqM] at h3
  · intro h
    have h2 : z ^ 2 ≤ (q ^ M) ^ 2 := sq_le_sq.mpr (by rwa [abs_of_nonneg hqM])
    linarith

/-! ## Layer 2: the physical update -/

/-- Single-qubit matrices. -/
abbrev Q2 := Matrix (Fin 2) (Fin 2) ℂ
/-- Two-qubit matrices; the index `(s, c)` is (system, cell). -/
abbrev Q4 := Matrix (Fin 2 × Fin 2) (Fin 2 × Fin 2) ℂ

/-- Diagonal qubit state with Bloch coordinate `z`: `diag((1 + z)/2, (1 - z)/2)`. -/
noncomputable def diagState (z : ℝ) : Q2 :=
  Matrix.diagonal ![((1 + z) / 2 : ℂ), ((1 - z) / 2 : ℂ)]

theorem diagState_one : diagState 1 = Matrix.diagonal ![1, 0] := by
  simp [diagState]

theorem diagState_neg_one : diagState (-1) = Matrix.diagonal ![0, 1] := by
  simp [diagState]

theorem diagState_zero : diagState 0 = (1 / 2 : ℂ) • (1 : Q2) := by
  ext i j
  fin_cases i <;> fin_cases j <;> simp [diagState]

/-- The two-qubit SWAP gate. -/
def swapGate : Q4 := fun a b => if a.1 = b.2 ∧ a.2 = b.1 then 1 else 0

/-- Partial swap `c · I + i s · SWAP`; with `c = cos η` and `s = sin η` this is
`U = cos η · I + i sin η · SWAP`. -/
noncomputable def pswap (c s : ℝ) : Q4 := (c : ℂ) • (1 : Q4) + (I * s) • swapGate

/-- Partial trace over the second (cell) qubit. -/
def ptrace2 (X : Q4) : Q2 := fun i i' => ∑ j, X (i, j) (i', j)

/-! Sanity checks on the encoding: the Kronecker index order is (system, cell),
`ptrace2` traces out the cell, `swapGate` exchanges the two factors, and the partial
swap is unitary. -/

theorem kron_apply_order (A B : Q2) (i₁ i₂ j₁ j₂ : Fin 2) :
    (A ⊗ₖ B) (i₁, i₂) (j₁, j₂) = A i₁ j₁ * B i₂ j₂ := rfl

theorem ptrace2_kron (A B : Q2) : ptrace2 (A ⊗ₖ B) = Matrix.trace B • A := by
  ext i j
  simp [ptrace2, Matrix.trace, Fin.sum_univ_two]
  ring

theorem swapGate_conj_kron (A B : Q2) : swapGate * (A ⊗ₖ B) * swapGate = B ⊗ₖ A := by
  ext ⟨a₁, a₂⟩ ⟨b₁, b₂⟩
  fin_cases a₁ <;> fin_cases a₂ <;> fin_cases b₁ <;> fin_cases b₂ <;>
    simp [swapGate, Matrix.mul_apply, Fintype.sum_prod_type, Fin.sum_univ_two] <;> ring

theorem pswap_unitary (c s : ℝ) (hcs : c ^ 2 + s ^ 2 = 1) :
    pswap c s * (pswap c s)ᴴ = 1 := by
  have hcs' : (c : ℂ) ^ 2 + (s : ℂ) ^ 2 = 1 := by exact_mod_cast hcs
  ext ⟨a₁, a₂⟩ ⟨b₁, b₂⟩
  fin_cases a₁ <;> fin_cases a₂ <;> fin_cases b₁ <;> fin_cases b₂ <;>
    simp [pswap, swapGate, Matrix.mul_apply, Fintype.sum_prod_type, Fin.sum_univ_two,
      Matrix.conjTranspose_apply, Matrix.one_apply] <;>
    first
    | (ring_nf; rw [Complex.I_sq]; linear_combination hcs')
    | ring_nf

/-- One collision of a diagonal system state with a diagonal cell: the reduced system
state is diagonal, with Bloch coordinate `c² z_ρ + (1 - c²) z_σ`. -/
theorem ptrace2_pswap (c s x y : ℝ) (hcs : c ^ 2 + s ^ 2 = 1) :
    ptrace2 (pswap c s * (diagState x ⊗ₖ diagState y) * (pswap c s)ᴴ) =
      diagState (c ^ 2 * x + (1 - c ^ 2) * y) := by
  have hcs' : (c : ℂ) ^ 2 + (s : ℂ) ^ 2 = 1 := by exact_mod_cast hcs
  ext i j
  fin_cases i <;> fin_cases j
  · simp [ptrace2, pswap, swapGate, diagState, Matrix.mul_apply, Fintype.sum_prod_type,
      Fin.sum_univ_two, Matrix.kroneckerMap_apply, Matrix.conjTranspose_apply,
      Matrix.one_apply, Matrix.diagonal_apply]
    linear_combination (-(s : ℂ) ^ 2 * (1 + (y : ℂ)) / 2) * Complex.I_sq +
      ((1 + (y : ℂ)) / 2) * hcs'
  · simp [ptrace2, pswap, swapGate, diagState, Matrix.mul_apply, Fintype.sum_prod_type,
      Fin.sum_univ_two, Matrix.kroneckerMap_apply, Matrix.conjTranspose_apply,
      Matrix.one_apply, Matrix.diagonal_apply]
  · simp [ptrace2, pswap, swapGate, diagState, Matrix.mul_apply, Fintype.sum_prod_type,
      Fin.sum_univ_two, Matrix.kroneckerMap_apply, Matrix.conjTranspose_apply,
      Matrix.one_apply, Matrix.diagonal_apply]
  · simp [ptrace2, pswap, swapGate, diagState, Matrix.mul_apply, Fintype.sum_prod_type,
      Fin.sum_univ_two, Matrix.kroneckerMap_apply, Matrix.conjTranspose_apply,
      Matrix.one_apply, Matrix.diagonal_apply]
    linear_combination (-(s : ℂ) ^ 2 * (1 - (y : ℂ)) / 2) * Complex.I_sq +
      ((1 - (y : ℂ)) / 2) * hcs'

/-- Layer 2 at angle `η`: the reduced update is `z ↦ q z_ρ + (1 - q) z_σ` with
`q = cos² η`. -/
theorem ptrace2_partialSwap (η x y : ℝ) :
    ptrace2 (pswap (Real.cos η) (Real.sin η) * (diagState x ⊗ₖ diagState y) *
        (pswap (Real.cos η) (Real.sin η))ᴴ) =
      diagState (step (Real.cos η ^ 2) x y) :=
  ptrace2_pswap _ _ x y (Real.cos_sq_add_sin_sq η)

/-! ## Layer 3: sequential exactness -/

section Sequential

variable {ι : Type*} [Fintype ι] [DecidableEq ι]

/-- Reduced state of the system (first factor) of a system–environment matrix. -/
def redS (X : Matrix (Fin 2 × ι) (Fin 2 × ι) ℂ) : Q2 := fun s s' => ∑ e, X (s, e) (s', e)

/-- Append a fresh cell in state `σ`, uncorrelated with the system and the environment. -/
def appendCell (X : Matrix (Fin 2 × ι) (Fin 2 × ι) ℂ) (σ : Q2) :
    Matrix (Fin 2 × (ι × Fin 2)) (Fin 2 × (ι × Fin 2)) ℂ :=
  fun a b => X (a.1, a.2.1) (b.1, b.2.1) * σ a.2.2 b.2.2

/-- A two-qubit gate `U` acting on the system and the new cell, and trivially on the
environment. -/
def onSC (U : Q4) : Matrix (Fin 2 × (ι × Fin 2)) (Fin 2 × (ι × Fin 2)) ℂ :=
  fun a b => if a.2.1 = b.2.1 then U (a.1, a.2.2) (b.1, b.2.2) else 0

/-- One collision: append a fresh cell `σ`, then apply `U` to the system and that cell. -/
def collide (U : Q4) (X : Matrix (Fin 2 × ι) (Fin 2 × ι) ℂ) (σ : Q2) :
    Matrix (Fin 2 × (ι × Fin 2)) (Fin 2 × (ι × Fin 2)) ℂ :=
  onSC U * appendCell X σ * (onSC U)ᴴ

/-- Sequential exactness: whatever correlations the system already has with the
environment, its reduced state after a collision with a fresh, uncorrelated cell depends
only on its reduced state before the collision. This holds for every gate `U`. -/
theorem redS_collide (U : Q4) (X : Matrix (Fin 2 × ι) (Fin 2 × ι) ℂ) (σ : Q2) :
    redS (collide U X σ) = ptrace2 (U * (redS X ⊗ₖ σ) * Uᴴ) := by
  ext s s'
  simp only [redS, collide, ptrace2, Matrix.mul_apply, Matrix.conjTranspose_apply, onSC,
    appendCell, Matrix.kroneckerMap_apply, Fintype.sum_prod_type]
  -- Collapse the environment guards onto the diagonal `e₁ = e₂ = e`.
  simp only [ite_mul, zero_mul, Finset.sum_ite_irrel, Finset.sum_const_zero,
    Finset.sum_ite_eq, Finset.mem_univ, ite_true, apply_ite star, star_zero, mul_ite, mul_zero]
  -- Distribute, then move the environment sum innermost on the left.
  simp only [Finset.sum_mul, Finset.mul_sum, mul_assoc]
  conv_lhs => rw [Finset.sum_comm]
  refine Finset.sum_congr rfl fun c _ => ?_
  conv_lhs => rw [Finset.sum_comm]
  refine Finset.sum_congr rfl fun s₂ _ => ?_
  conv_lhs => rw [Finset.sum_comm]
  refine Finset.sum_congr rfl fun c₂ _ => ?_
  conv_lhs => rw [Finset.sum_comm]
  refine Finset.sum_congr rfl fun s₁ _ => ?_
  conv_lhs => rw [Finset.sum_comm]

end Sequential

/-- Index type for the first `k` cells, which have already met the system. -/
def Env : ℕ → Type
  | 0 => Unit
  | k + 1 => Env k × Fin 2

instance instFintypeEnv : (k : ℕ) → Fintype (Env k)
  | 0 => inferInstanceAs (Fintype Unit)
  | k + 1 => @instFintypeProd (Env k) (Fin 2) (instFintypeEnv k) _

instance instDecidableEqEnv : (k : ℕ) → DecidableEq (Env k)
  | 0 => inferInstanceAs (DecidableEq Unit)
  | k + 1 => @instDecidableEqProd (Env k) (Fin 2) (instDecidableEqEnv k) _

theorem card_env_zero : Fintype.card (Env 0) = 1 := rfl

/-- Joint state of the system and the first `k` cells after the system has met each of
them once, starting from system state `ρ₀` and cells in states `σ 0, σ 1, …`. Each cell
is appended in its initial state when the system reaches it. This is the appending
model: a cell the system has not yet met is untouched and uncorrelated, so it matches
the circuit with all cells present from the start, but that equality is a modelling
step, not a theorem in this file. -/
def circuit (U : Q4) (ρ₀ : Q2) (σ : ℕ → Q2) :
    (k : ℕ) → Matrix (Fin 2 × Env k) (Fin 2 × Env k) ℂ
  | 0 => fun a b => ρ₀ a.1 b.1
  | k + 1 => collide U (circuit U ρ₀ σ k) (σ k)

/-- With partial swaps and diagonal cells, the system's reduced state after `k`
collisions is the diagonal state whose Bloch coordinate follows the scalar recursion
`run`; every system–cell correlation is retained in the joint state. -/
theorem redS_circuit (c s : ℝ) (hcs : c ^ 2 + s ^ 2 = 1) (z₀ : ℝ) (zc : ℕ → ℝ) (k : ℕ) :
    redS (circuit (pswap c s) (diagState z₀) (fun j => diagState (zc j)) k) =
      diagState (run (c ^ 2) zc z₀ k) := by
  induction k with
  | zero =>
      ext s s'
      simp [redS, circuit, run, card_env_zero]
  | succ k ih =>
      have h := redS_collide (pswap c s)
        (circuit (pswap c s) (diagState z₀) (fun j => diagState (zc j)) k) (diagState (zc k))
      rw [ih, ptrace2_pswap c s _ _ hcs] at h
      exact h

/-! ## The first use of each initializer -/

/-- The partial swap at angle `η`. -/
noncomputable abbrev U (η : ℝ) : Q4 := pswap (Real.cos η) (Real.sin η)

/-- First output for the designated initializer `τ^{⊗M}`: Bloch coordinate `q^M`. -/
theorem firstOutput_designated (η : ℝ) (M : ℕ) :
    redS (circuit (U η) (diagState 1) (fun j => diagState (designatedCells j)) M) =
      diagState ((Real.cos η ^ 2) ^ M) := by
  rw [redS_circuit _ _ (Real.cos_sq_add_sin_sq η), run_designated]

/-- First output for the product initializer `|1⟩⟨1| ⊗ τ^{⊗(M-1)}` with `M = k + 1`:
Bloch coordinate `(2q - 1) q^(M-1)`. -/
theorem firstOutput_product (η : ℝ) (k : ℕ) :
    redS (circuit (U η) (diagState 1) (fun j => diagState (productCells j)) (k + 1)) =
      diagState ((2 * Real.cos η ^ 2 - 1) * (Real.cos η ^ 2) ^ k) := by
  rw [redS_circuit _ _ (Real.cos_sq_add_sin_sq η), run_product]

/-- First output for the adjusted initializer with `M = k + 2` cells, when `q < 1`:
exactly `τ = I/2`. -/
theorem firstOutput_adjusted (η : ℝ) (hq1 : Real.cos η ^ 2 < 1) (k : ℕ) :
    redS (circuit (U η) (diagState 1)
        (fun j => diagState (adjustedCells (Real.cos η ^ 2) j)) (k + 2)) =
      (1 / 2 : ℂ) • (1 : Q2) := by
  rw [redS_circuit _ _ (Real.cos_sq_add_sin_sq η), run_adjusted hq1, diagState_zero]

/-- The first-use capability step. For every strict partial swap `0 < q < 1` and every
`M ≥ 2`, there is a product initializer whose first cell is `|1⟩⟨1|`, whose cells are all
valid diagonal qubit states, and whose first output (in the appending model) has Bloch
coordinate `z` with `|z| ≤ q^M`. With the closed-form fidelity `fidTau`, this is
membership in `V_{e_M}`. -/
theorem capable_first_use {η : ℝ} (hq0 : 0 < Real.cos η ^ 2) (hq1 : Real.cos η ^ 2 < 1)
    (M : ℕ) (hM : 2 ≤ M) :
    ∃ zc : ℕ → ℝ, zc 0 = -1 ∧ (∀ j, zc j ∈ Set.Icc (-1 : ℝ) 1) ∧
      ∃ z : ℝ,
        redS (circuit (U η) (diagState 1) (fun j => diagState (zc j)) M) = diagState z ∧
        |z| ≤ (Real.cos η ^ 2) ^ M ∧
        1 - eM (Real.cos η ^ 2) M ≤ fidTau z := by
  set q := Real.cos η ^ 2 with hqdef
  have hqM1 : q ^ M ≤ 1 := pow_le_one₀ hq0.le hq1.le
  by_cases h : 1 / 3 ≤ q
  · obtain ⟨k, rfl⟩ : ∃ k, M = k + 1 := ⟨M - 1, by omega⟩
    refine ⟨productCells, by simp [productCells], ?_, (2 * q - 1) * q ^ k, ?_, ?_, ?_⟩
    · intro j
      by_cases hj : j = 0 <;> simp [productCells, hj]
    · exact firstOutput_product η k
    · exact (product_within_iff hq0 k).mpr ⟨h, hq1.le⟩
    · have hz := (product_within_iff hq0 k).mpr ⟨h, hq1.le⟩
      exact (fidTau_ge_iff hq0.le (hz.trans hqM1) _).mpr hz
  · replace h : q < 1 / 3 := lt_of_not_ge h
    obtain ⟨k, rfl⟩ : ∃ k, M = k + 2 := ⟨M - 2, by omega⟩
    have hz2 := zTwo_mem hq0 h
    refine ⟨adjustedCells q, by simp [adjustedCells], ?_, 0, ?_, ?_, ?_⟩
    · intro j
      by_cases hj0 : j = 0
      · simp [adjustedCells, hj0]
      · by_cases hj1 : j = 1
        · subst hj1
          simp only [adjustedCells, one_ne_zero, ite_false, ite_true]
          exact ⟨by linarith [hz2.1], hz2.2⟩
        · simp [adjustedCells, hj0, hj1]
    · rw [redS_circuit _ _ (Real.cos_sq_add_sin_sq η), run_adjusted hq1]
    · simpa using pow_nonneg hq0.le (k + 2)
    · exact (fidTau_ge_iff hq0.le (by simp) _).mpr (by simpa using pow_nonneg hq0.le (k + 2))

end CapableState
