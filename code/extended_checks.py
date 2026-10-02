"""Reproducible diagnostics accompanying the homogenizer manuscript.

The finite-grid rewrite certificates and the pi/4 example use exact discrete
arithmetic. Density-matrix checks use floating point and are not formal proofs.
No tensor-product approximation is made to the reused memory.
"""

from __future__ import annotations

import json
import math
from fractions import Fraction
from pathlib import Path

import numpy as np
from exact_homogenizer import (
    TAU,
    P,
    collide,
    power,
    reduce,
    reused_memory,
    self_test,
    swap_permutation,
    sweep,
    tensor,
    trace_distance,
)


def exchange_certificate(n: int, m: int) -> list[dict]:
    """Adjacent commutations converting row order to column order.

    Gates (i,j) and (k,l) have disjoint supports exactly when i != k and
    j != l (the input and reservoir wire sets are distinct).
    The function emits and checks every local rewrite, without matrices.
    """
    if n < 0 or m < 0:
        raise ValueError("grid dimensions must be nonnegative")
    word = [(i, j) for i in range(n) for j in range(m)]
    target = [(i, j) for j in range(m) for i in range(n)]
    steps = []
    for pos, gate in enumerate(target):
        old = word.index(gate, pos)
        while old > pos:
            left, right = word[old - 1], word[old]
            assert left[0] != right[0] and left[1] != right[1]
            steps.append(dict(position=old - 1, left=list(left), right=list(right)))
            word[old - 1], word[old] = right, left
            old -= 1
    assert word == target
    assert len(steps) == math.comb(n, 2) * math.comb(m, 2)
    return steps


def exact_pi_over_four_example(m: int = 2, n: int = 3) -> dict:
    """Gaussian-integer numerator / integer denominator density matrices.

    At eta=pi/4 the update is (rho+P rho P+i(P rho-rho P))/2.
    Integer overflow is avoided by numpy object arrays (Python integers).
    """
    answers = {}
    for forward in (True, False):
        d = 1 << m
        real = np.eye(d, dtype=object) if forward else np.zeros((d, d), dtype=object)
        if not forward:
            real[0, 0] = 1
        imag = np.zeros((d, d), dtype=object)
        denominator = d if forward else 1
        for _ in range(n):
            inp = np.diag([1, 0]).astype(object) if forward else np.eye(2, dtype=object)
            re, im = np.kron(inp, real), np.kron(inp, imag)
            denominator *= 1 if forward else 2
            for j in range(m):
                perm = swap_permutation(m + 1, 0, j + 1)
                # i*(X+iY) = -Y+iX
                comm_re = re[perm, :] - re[:, perm]
                comm_im = im[perm, :] - im[:, perm]
                re, im = (
                    re + re[np.ix_(perm, perm)] - comm_im,
                    im + im[np.ix_(perm, perm)] + comm_re,
                )
                denominator *= 2
            # Output excited-state probability; all local coherences vanish.
            excited = Fraction(sum(re[d + k, d + k] for k in range(d)), denominator)
            assert all(
                sum(im[(s * d + k), (t * d + k)] for k in range(d)) == 0
                for s in range(2)
                for t in range(2)
            )
            assert sum(re[k, d + k] for k in range(d)) == 0
            real = re[:d, :d] + re[d:, d:]
            imag = im[:d, :d] + im[d:, d:]
        error = abs(excited - Fraction(1, 2)) if forward else excited
        answers["pure_to_mixed" if forward else "mixed_to_pure"] = str(error)
    if (m, n) == (2, 3):
        assert answers == {"pure_to_mixed": "25/64", "mixed_to_pure": "3/8"}
    return answers


def onewire(rho: np.ndarray, gate: np.ndarray, wire: int) -> np.ndarray:
    q = rho.shape[0].bit_length() - 1
    unitary = tensor(*[gate if j == wire else np.eye(2) for j in range(q)])
    return unitary @ rho @ unitary.conj().T


def clockless_outputs(n: int, m: int) -> tuple[np.ndarray, np.ndarray]:
    """Fixed H(input), CNOT(input -> head), H(head), cyclic-memory device.

    It produces min(n,2m) independent maximally mixed qubits and attains
    T=max(0,1-2**(2m-n)) at every later horizon, for the known input P.
    """
    if n < 1 or m < 0 or n + m > 10:
        raise ValueError("require n>=1, m>=0 and n+m<=10")
    rho = tensor(power(P, n), power(TAU, m))
    h = np.array([[1, 1], [1, -1]], complex) / np.sqrt(2)
    q = n + m
    for i in range(n):
        rho = onewire(rho, h, i)
        if m:
            x = np.arange(1 << q)
            # CNOT is a self-inverse computational basis permutation.
            perm = x ^ (((x >> (q - 1 - i)) & 1) << (q - 1 - n))
            rho = rho[np.ix_(perm, perm)]
            rho = onewire(rho, h, n)
            # Move the visited head to the tail; adjacent SWAPs realize it.
            for j in range(m - 1):
                perm = swap_permutation(q, n + j, n + j + 1)
                rho = rho[np.ix_(perm, perm)]
    return reduce(rho, list(range(n))), reduce(rho, list(range(n, n + m)))


def weak_catalytic_design(n: int, epsilon: float) -> dict:
    """Analytic sufficient resources, no large state simulation.

    Ensures every forward marginal T<=epsilon and final whole-memory
    T<=epsilon/2 by the proved chi-square / telescoping bounds.
    """
    if n < 1 or not 0 < epsilon < 0.5:
        raise ValueError("require n>=1 and 0<epsilon<1/2")
    ell = math.log(1 / epsilon)
    m = math.ceil(max(2 * ell, n * n * ell * ell / (epsilon * epsilon)))
    u = ell / m
    eta = math.asin(math.sqrt(u))
    one = 0.5 * math.exp(m * math.log1p(-u))
    chi = -math.expm1(m * math.log1p(-u * u))
    memory_bound = 0.5 * n * math.sqrt(max(0, chi))
    assert one + memory_bound <= epsilon * (1 + 1e-10)
    return dict(
        n=n,
        epsilon=epsilon,
        m=m,
        eta=eta,
        marginal_error_upper_bound=one + memory_bound,
        whole_memory_error_upper_bound=memory_bound,
    )


def frobenius_mixing_diagnostic(b: np.ndarray, n: int, eta: float, max_l: int = 10000) -> dict:
    """A small-n, floating-point diagnostic for the analytic mixing bound.

    With exact arithmetic, kappa=sqrt(2**n)*||Phi**L-R_b||_F<1 certifies
    T(Phi**M(rho),b**n)<=kappa**floor(M/L). Here the matrices and norm
    are floating point, so the returned values are NOT interval certificates.
    """
    if not 1 <= n <= 3:
        raise ValueError("diagnostic is restricted to 1<=n<=3")
    d = 1 << n
    basis = np.eye(d * d, dtype=complex)
    phi = np.column_stack(
        [sweep(basis[:, j].reshape(d, d), b, eta).reshape(-1) for j in range(d * d)]
    )
    replacer = np.outer(power(b, n).reshape(-1), np.eye(d).reshape(-1))
    mat = np.eye(d * d, dtype=complex)
    for length in range(1, max_l + 1):
        mat = phi @ mat
        kappa = float(np.sqrt(d) * np.linalg.norm(mat - replacer, "fro"))
        if kappa < 0.9:
            return dict(n=n, eta=eta, L=length, kappa=kappa, rigorous_interval_certificate=False)
    raise RuntimeError("no contraction found within max_l")


def capability_counterexamples() -> dict:
    """Check product initializers orthogonal to the pure fixed point.

    The proof is analytic in the manuscript. These finite checks include
    strong and weak coupling and retain the full first-use density matrix.
    """
    count = 0
    worst_excess = 0.0
    for q in (0.05, 0.2, 0.333333, 0.4, 0.5, 0.8, 0.999):
        for m in (2, 3, 4):
            excited = np.diag([0.0, 1.0]).astype(complex)
            if q >= 1 / 3:
                mem = tensor(excited, power(TAU, m - 1))
            else:
                z = q * (1 - 2 * q) / (1 - q)
                second = np.diag([(1 + z) / 2, (1 - z) / 2])
                mem = tensor(excited, second, power(TAU, m - 2))
            assert mem[0, 0] == 0
            state = tensor(P, mem)
            eta = math.acos(math.sqrt(q))
            for j in range(m):
                state = collide(state, 0, j + 1, eta)
            excess = trace_distance(reduce(state, [0]), TAU) - q**m / 2
            assert excess < 1e-12
            worst_excess = max(worst_excess, excess)
            count += 1
    return dict(cases=count, max_positive_capability_excess=worst_excess)


def run_checks() -> dict:
    grid_count = 0
    for n in range(9):
        for m in range(9):
            exchange_certificate(n, m)
            grid_count += 1
    max_chi = 0.0
    for m in range(1, 6):
        for eta in [0.2, 0.5, 0.8, 1.2]:
            _, mem = reused_memory(P, TAU, 1, m, eta)
            actual = (1 << m) * np.trace(mem @ mem).real - 1
            expected = 1 - (1 - np.sin(eta) ** 4) ** m
            max_chi = max(max_chi, abs(actual - expected))
    optimal = []
    for m in range(3):
        for n in range(1, 7):
            out, mem = clockless_outputs(n, m)
            actual = trace_distance(out, power(TAU, n))
            predicted = max(0.0, 1 - 2.0 ** (2 * m - n))
            assert abs(actual - predicted) < 1e-11
            assert trace_distance(mem, power(TAU, m)) < 1e-11
            optimal.append(dict(m=m, n=n, error=actual, predicted=predicted))
    return dict(
        base_checks=self_test(),
        product_capability_counterexamples=capability_counterexamples(),
        exact_rewrite_grids_checked=grid_count,
        exact_rational_example=exact_pi_over_four_example(),
        max_chi_square_identity_error=max_chi,
        clockless_optimality_checks=optimal,
        weak_catalytic_design=weak_catalytic_design(100, 0.01),
        floating_mixing_diagnostics=[
            frobenius_mixing_diagnostic(P, 2, 0.5),
            frobenius_mixing_diagnostic(TAU, 2, 0.5),
        ],
    )


if __name__ == "__main__":
    result = run_checks()
    path = Path(__file__).with_name("extended_checks.json")
    path.write_text(json.dumps(result, indent=2))
    print(json.dumps(result, indent=2))
