"""Adversarial diagnostics for Returning a quantum constructor.

This program is not an independent-agent review or a formal proof. The symbolic
singlet identity and dyadic-shift tests are exact. Matrix diagnostics are floating
point, retain all correlations, and cannot establish asymptotic theorems.
Run: python return_checks.py
"""

from __future__ import annotations

import json
import math
from fractions import Fraction
from pathlib import Path

import numpy as np
import sympy as sp
from exact_homogenizer import (
    TAU,
    P,
    outputs_by_columns,
    power,
    reduce,
    reused_memory,
    self_test,
    sweep,
    tensor,
    trace_distance,
)
from extended_checks import (
    capability_counterexamples,
    clockless_outputs,
    exact_pi_over_four_example,
    exchange_certificate,
)

PAULI = [
    np.array([[0, 1], [1, 0]], complex),
    np.array([[0, -1j], [1j, 0]], complex),
    np.diag([1, -1]).astype(complex),
]


def explicit_swap(qubits: int, first: int, second: int) -> np.ndarray:
    """Independent dense implementation, using bit strings rather than masks."""
    d = 2**qubits
    result = np.zeros((d, d), complex)
    for col in range(d):
        bits = list(format(col, f"0{qubits}b"))
        bits[first], bits[second] = bits[second], bits[first]
        result[int("".join(bits), 2), col] = 1
    return result


def dense_sweep(block: np.ndarray, fresh: np.ndarray, eta: float) -> np.ndarray:
    k = int(round(math.log2(block.shape[0])))
    unitary = np.eye(2 ** (k + 1), dtype=complex)
    for j in range(1, k + 1):
        swap = explicit_swap(k + 1, 0, j)
        gate = np.cos(eta) * np.eye(len(swap)) + 1j * np.sin(eta) * swap
        unitary = gate @ unitary
    total = unitary @ np.kron(fresh, block) @ unitary.conj().T
    # Direct block trace, not the tensor-axis partial trace of the older code.
    d = block.shape[0]
    return total[:d, :d] + total[d:, d:]


def exact_singlet_identity() -> dict:
    c, s = sp.symbols("c s", real=True)
    S1 = sp.Matrix(explicit_swap(3, 0, 1).real.astype(int))
    S2 = sp.Matrix(explicit_swap(3, 0, 2).real.astype(int))
    pair_swap = sp.Matrix(explicit_swap(2, 0, 1).real.astype(int))
    singlet = (sp.eye(4) - pair_swap) / 2
    W = (c * sp.eye(8) + sp.I * s * S2) * (c * sp.eye(8) + sp.I * s * S1)
    lifted = W.conjugate().T * sp.kronecker_product(sp.eye(2), singlet) * W
    adjoint = (lifted[:4, :4] + lifted[4:, 4:]) / 2
    target = s**4 * sp.eye(4) / 4 + (1 - s**4) * singlet
    # Reduce every entry modulo c^2+s^2-1 in the polynomial ring Q(i)[c,s].
    basis = sp.groebner([c * c + s * s - 1], c, s, extension=sp.I)
    residuals = [sp.expand(basis.reduce(sp.expand(z))[1]) for z in adjoint - target]
    assert all(z == 0 for z in residuals)
    return {
        "arithmetic": "symbolic polynomial reduction over Q(i)",
        "identity": "Phi_tau,2^*(Pminus)=s^4 I/4+(1-s^4)Pminus",
        "zero_residual_entries": len(residuals),
    }


def dyadic_spectrum(m: int, n: int) -> list[Fraction]:
    if not (1 <= n and m >= 2 * n):
        raise ValueError("require n >= 1 and m >= 2*n")
    k = m - n
    C = Fraction(2, k + 2)
    p = [Fraction(0)] * (2**m)
    p[0] = C
    for i in range(2, 2**k + 1):
        ceil_log = (i - 1).bit_length()
        p[i - 1] = C / Fraction(2**ceil_log)
    assert sum(p) == 1
    return p


def dyadic_checks() -> list[dict]:
    rows = []
    for n, m in [(1, 3), (1, 7), (2, 6), (2, 9), (3, 8), (3, 10), (4, 10)]:
        p = dyadic_spectrum(m, n)
        distances = []
        for t in range(1, n + 1):
            after = [p[i // (2**t)] / (2**t) for i in range(2**m)]
            assert sum(after) == 1
            distance = sum(abs(x - y) for x, y in zip(p, after)) / 2
            expected = Fraction(t, m - n + 2)
            assert distance == expected
            distances.append(str(distance))
        rows.append(dict(n=n, m=m, prefix_distances=distances))
    return rows


def ng_prefix_counterexample() -> dict:
    # Ng et al.'s endpoint-optimal d=4, D=4^3 catalyst: C=4/10.
    n, a = 2, 3
    d = 2**n
    m = n * a
    C = Fraction(d, 1 + (d - 1) * a)
    p = [Fraction(0)] * (2**m)
    p[0] = C
    for i in range(2, d ** (a - 1) + 1):
        shell = 0
        cap = 1
        while cap < i:
            shell += 1
            cap *= d
        p[i - 1] = C / (d**shell)
    assert sum(p) == 1
    distances = []
    for t in range(1, n + 1):
        after = [p[i // (2**t)] / (2**t) for i in range(2**m)]
        distances.append(sum(abs(x - y) for x, y in zip(p, after)) / 2)
    assert distances == [Fraction(2, 5), Fraction(3, 10)]
    return dict(
        n=n,
        m=m,
        prefix_distances=list(map(str, distances)),
        interpretation="Endpoint-optimal spectrum does not optimize prefix return.",
    )


def entropy(rho: np.ndarray) -> float:
    eig = np.maximum(np.linalg.eigvalsh((rho + rho.conj().T) / 2), 0)
    eig = eig[eig > 1e-14]
    return float(-np.sum(eig * np.log2(eig)))


def entropy_checks(rng: np.random.Generator) -> dict:
    worst = 0.0
    bound_violation = 0.0
    for n, m in [(1, 1), (1, 2), (2, 1), (2, 2)]:
        d = 2**m
        total = 2 ** (n + m)
        for _ in range(6):
            x = rng.normal(size=(d, d)) + 1j * rng.normal(size=(d, d))
            sigma = x @ x.conj().T
            sigma /= np.trace(sigma)
            z = rng.normal(size=(total, total)) + 1j * rng.normal(size=(total, total))
            U, _ = np.linalg.qr(z)
            joint = U @ np.kron(power(TAU, n), sigma) @ U.conj().T
            out = reduce(joint, list(range(n)))
            mem = reduce(joint, list(range(n, n + m)))
            mutual = entropy(out) + entropy(mem) - entropy(joint)
            residual = (n - entropy(out)) - (entropy(mem) - entropy(sigma) - mutual)
            worst = max(worst, abs(residual))
            delta = trace_distance(mem, sigma)
            if 0 < delta <= 0.5:
                h = -delta * math.log2(delta) - (1 - delta) * math.log2(1 - delta)
                bound_violation = max(
                    bound_violation, n - entropy(out) - h - delta * math.log2(d - 1)
                )
    assert worst < 1e-10 and bound_violation < 1e-10
    return dict(
        identity_max_abs_residual=worst,
        continuity_bound_violation=bound_violation,
        random_unitary_samples=24,
    )


def moment_checks() -> dict:
    purity_res = 0.0
    fourth_violation = 0.0
    lower_violation = 0.0
    samples = 0
    dense_difference = 0.0
    for M in range(1, 7):
        for u in [0.005, 0.03, 0.1, 0.25, 0.6]:
            eta = math.asin(math.sqrt(u))
            _, mem = reused_memory(P, TAU, 1, M, eta)
            D = 2**M
            A = D * mem - np.eye(D)
            a = float(np.trace(A @ A).real / D)
            b = float(np.trace(A @ A @ A @ A).real / D)
            expected = 1 - (1 - u * u) ** M
            purity_res = max(purity_res, abs(a - expected))
            fourth_violation = max(fourth_violation, b - (3 * M * M + M) * u**4)
            delta = trace_distance(mem, power(TAU, M))
            if delta <= 0.25:
                lower_violation = max(
                    lower_violation, math.sqrt(M * u * u) / (8 * math.sqrt(2)) - delta
                )
            # Direct full-unitary implementation of the memory sweep.
            dense = dense_sweep(power(TAU, M), P, eta)
            dense_difference = max(dense_difference, float(np.linalg.norm(mem - dense)))
            samples += 1
    assert purity_res < 1e-10 and fourth_violation < 1e-9 and lower_violation < 1e-10
    assert dense_difference < 1e-10
    return dict(
        samples=samples,
        purity_max_abs_residual=purity_res,
        fourth_moment_upper_bound_violation=fourth_violation,
        trace_lower_bound_violation=lower_violation,
        separate_dense_implementation_max_difference=dense_difference,
    )


def collective_generator(X: np.ndarray, n: int) -> np.ndarray:
    L = np.zeros_like(X)
    for pauli in PAULI:
        J = sum(
            (tensor(*[pauli / 2 if k == j else np.eye(2) for k in range(n)]) for j in range(n)),
            start=np.zeros_like(X),
        )
        comm = J @ X - X @ J
        L -= 0.5 * (J @ comm - comm @ J)
    return L


def tnorm(X: np.ndarray) -> float:
    return float(np.abs(np.linalg.eigvalsh((X + X.conj().T) / 2)).sum())


def collective_checks(rng: np.random.Generator) -> dict:
    ratios = []
    rank_one_res = 0.0
    projection_norms = []
    for n in [1, 2, 3, 4]:
        d = 2**n
        X = rng.normal(size=(d, d)) + 1j * rng.normal(size=(d, d))
        X = X + X.conj().T
        X /= tnorm(X)
        L = collective_generator(X, n)
        eta = 0.003
        rem = sweep(X, TAU, eta) - X - eta * eta * L
        ratios.append(tnorm(rem) / ((4 / 3) * n**3 * eta**3))
        local = tensor(PAULI[2], *([np.eye(2)] * (n - 1)))
        rank_one_res = max(
            rank_one_res, float(np.linalg.norm(collective_generator(local, n) + local))
        )
        # Exact symmetric-space rank-one projection, represented in full space.
        Pi = np.zeros((d, d), complex)
        for h in range(n + 1):
            v = np.array([int(i.bit_count() == h) for i in range(d)], complex)
            v /= np.linalg.norm(v)
            Pi += np.outer(v, v.conj())
        Jz = sum(
            (tensor(*[PAULI[2] / 2 if k == j else np.eye(2) for k in range(n)]) for j in range(n)),
            start=np.zeros((d, d), complex),
        )
        rank_one = 6 / ((n + 1) * (n + 2)) * Jz @ Pi
        projection_norms.append(tnorm(rank_one))
        assert tnorm(rank_one) <= 1.5 + 1e-12
        assert abs(np.trace(local @ (power(P, n) - rank_one))) < 1e-12
    assert max(ratios) <= 1 + 1e-8 and rank_one_res < 1e-10
    return dict(
        remainder_to_proved_upper_bound=ratios,
        rank_one_generator_max_residual=rank_one_res,
        projected_initial_trace_norms=projection_norms,
    )


def singlet_numerics() -> dict:
    Pminus = (np.eye(4) - explicit_swap(2, 0, 1)) / 2
    worst = 0.0
    obstruction_violation = 0.0
    rows = []
    for M in [1, 2, 3, 5, 8]:
        for u in [0.01, 0.08, 0.3, 0.7]:
            eta = math.asin(math.sqrt(u))
            out = outputs_by_columns(P, TAU, 2, M, eta)
            z = float(np.trace(Pminus @ out).real)
            a = 1 - (1 - u * u) ** M
            worst = max(worst, abs(z - a / 4))
            # 1-use memory calculation is exponentially sized and kept small.
            _, mem = reused_memory(P, TAU, 1, M, eta)
            delta = trace_distance(mem, power(TAU, M))
            err = trace_distance(out, power(TAU, 2))
            obstruction_violation = max(obstruction_violation, 0.25 - delta / 2 - err)
    eta = 0.05
    M = 1200
    out = outputs_by_columns(P, TAU, 2, M, eta)
    rows.append(
        dict(
            n=2,
            M=M,
            eta=eta,
            marginal_errors=[trace_distance(reduce(out, [j]), TAU) for j in range(2)],
            joint_error=trace_distance(out, power(TAU, 2)),
            singlet_probability=float(np.trace(Pminus @ out).real),
            analytic_prefix_return_upper_bound=math.sin(eta) ** 2 * math.sqrt(M),
        )
    )
    assert worst < 1e-9 and obstruction_violation < 1e-9
    return dict(
        identity_max_abs_residual=worst,
        joint_obstruction_violation=obstruction_violation,
        illustrative_computed_outputs=rows,
    )


def reverse_first_use_checks() -> dict:
    residual = 0.0
    cases = 0
    for M in range(1, 7):
        for eta in [0.03, 0.2, 0.5, 1.1]:
            out, mem = reused_memory(TAU, P, 1, M, eta)
            qM = math.cos(eta) ** (2 * M)
            e = trace_distance(out, P)
            d = trace_distance(mem, power(P, M))
            residual = max(residual, abs(e - qM / 2), abs(d - (1 - qM) / 2))
            cases += 1
    assert residual < 1e-10
    return dict(
        cases=cases,
        maximum_absolute_residual=residual,
        identity="output_error + return_error = 1/2 at one reverse use",
    )


def main() -> dict:
    rng = np.random.default_rng(20261002)
    grids = 0
    for n in range(9):
        for m in range(9):
            exchange_certificate(n, m)
            grids += 1
    caps = capability_counterexamples()
    optimal = []
    for n, m in [(1, 1), (2, 1), (3, 1), (4, 2), (5, 2), (6, 3)]:
        out, mem = clockless_outputs(n, m)
        e = trace_distance(out, power(TAU, n))
        expected = max(0.0, 1 - 2 ** (2 * m - n))
        r = trace_distance(mem, power(TAU, m))
        assert abs(e - expected) < 1e-10 and r < 1e-10
        optimal.append(dict(n=n, m=m, joint_error=e, expected_error=expected, return_error=r))
    report = dict(
        status="Executed diagnostics, not an independent review or a formal proof",
        seed=20261002,
        exact_grid_certificates=grids,
        exact_pi_over_four=exact_pi_over_four_example(),
        legacy_full_correlation_tests=self_test(),
        product_capability_counterexamples=caps,
        symbolic_singlet=exact_singlet_identity(),
        exact_dyadic_shift=dyadic_checks(),
        exact_endpoint_vs_prefix_counterexample=ng_prefix_counterexample(),
        entropy=entropy_checks(rng),
        moments=moment_checks(),
        collective_expansion=collective_checks(rng),
        two_output_checks=singlet_numerics(),
        reverse_first_use=reverse_first_use_checks(),
        optimal_forward_device=optimal,
    )
    out = Path(__file__).parent.parent / "verification" / "return_checks.json"
    out.write_text(json.dumps(report, indent=2))
    print(json.dumps(report, indent=2))
    return report


if __name__ == "__main__":
    main()
