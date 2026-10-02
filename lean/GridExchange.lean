import Mathlib

/-!
# Disjoint-gate exchange certificates

Machine-checked with Lean 4 and Mathlib at the versions pinned in
`lean-toolchain` and `lake-manifest.json`; build with `lake build`.

The scope is the algebraic soundness of adjacent disjoint-wire
exchanges, plus a concrete 2 x 2 certificate. It does not prove the existence
of a certificate for every rectangle, or any density-matrix/entropy theorem.

`U` may take values in a noncommutative monoid. The only commutation hypothesis
is that gates on different input wires AND different reservoir wires commute.
The tensor-product representation establishing that physical hypothesis is
outside this file.
-/

namespace Homogenizer

abbrev Gate := Nat × Nat

def Separate (a b : Gate) : Prop := a.1 ≠ b.1 ∧ a.2 ≠ b.2

instance : DecidableRel Separate := fun a b =>
  inferInstanceAs (Decidable (a.1 ≠ b.1 ∧ a.2 ≠ b.2))

def eval {G : Type*} [Monoid G] (U : Gate → G) : List Gate → G
  | [] => 1
  | a :: xs => U a * eval U xs

theorem eval_append {G : Type*} [Monoid G] (U : Gate → G)
    (xs ys : List Gate) : eval U (xs ++ ys) = eval U xs * eval U ys := by
  induction xs with
  | nil => simp [eval]
  | cons a xs ih => simp [eval, ih, mul_assoc]

/-- A certificate whose elementary move exchanges adjacent disjoint gates. -/
inductive Exchange : List Gate → List Gate → Prop
  | refl (xs : List Gate) : Exchange xs xs
  | swap (pre post : List Gate) (a b : Gate) (h : Separate a b) :
      Exchange (pre ++ a :: b :: post) (pre ++ b :: a :: post)
  | trans {xs ys zs : List Gate} :
      Exchange xs ys → Exchange ys zs → Exchange xs zs

/-- Every valid exchange certificate preserves the evaluated product. -/
theorem exchange_sound {G : Type*} [Monoid G] (U : Gate → G)
    (hcomm : ∀ a b, Separate a b → U a * U b = U b * U a)
    {xs ys : List Gate} (certificate : Exchange xs ys) :
    eval U xs = eval U ys := by
  induction certificate with
  | refl xs => rfl
  | swap pre post a b h =>
      have hab : U a * U b = U b * U a := hcomm a b h
      have hc := congrArg (fun z => eval U pre * (z * eval U post)) hab
      simpa only [eval_append, eval, mul_assoc] using hc
  | trans hxy hyz ihxy ihyz => exact ihxy.trans ihyz

def row2 : List Gate := [(0, 0), (0, 1), (1, 0), (1, 1)]
def column2 : List Gate := [(0, 0), (1, 0), (0, 1), (1, 1)]

theorem rectangle_two_by_two : Exchange row2 column2 := by
  exact Exchange.swap [(0, 0)] [(1, 1)] (0, 1) (1, 0) (by decide)

theorem rectangle_two_by_two_sound {G : Type*} [Monoid G] (U : Gate → G)
    (hcomm : ∀ a b, Separate a b → U a * U b = U b * U a) :
    eval U row2 = eval U column2 :=
  exchange_sound U hcomm rectangle_two_by_two

end Homogenizer
