"""Exact density-matrix tests of a reused quantum homogenizer.

All reservoir/output correlations are retained. Exponential memory in the number
of retained qubits: intended for small-system verification, not extrapolation.
Fidelity means squared Uhlmann fidelity; trace distance includes the factor 1/2.
"""

from __future__ import annotations

import json
from functools import lru_cache
from pathlib import Path

import numpy as np
from numpy.typing import NDArray

Array = NDArray[np.complex128]
P = np.diag([1.0, 0.0]).astype(complex)
TAU = np.eye(2, dtype=complex) / 2


def tensor(*args: Array) -> Array:
    out = np.ones((1, 1), dtype=complex)
    for a in args:
        out = np.kron(out, a)
    return out


def power(a: Array, n: int) -> Array:
    if n < 0:
        raise ValueError("n must be nonnegative")
    return tensor(*([a] * n))


@lru_cache(None)
def swap_permutation(q: int, i: int, j: int) -> NDArray[np.int64]:
    if not (0 <= i < q and 0 <= j < q):
        raise ValueError("wire index out of range")
    x = np.arange(1 << q, dtype=np.int64)
    bit = ((x >> (q - 1 - i)) ^ (x >> (q - 1 - j))) & 1
    return x ^ (bit << (q - 1 - i)) ^ (bit << (q - 1 - j))


def collide(rho: Array, i: int, j: int, eta: float) -> Array:
    q = (rho.shape[0]).bit_length() - 1
    if rho.shape != (1 << q, 1 << q):
        raise ValueError("rho must be square with power-of-two dimension")
    perm = swap_permutation(q, i, j)
    c, s = np.cos(eta), np.sin(eta)
    left, right = rho[perm, :], rho[:, perm]
    return c * c * rho + s * s * rho[np.ix_(perm, perm)] + 1j * c * s * (left - right)


def reduce(rho: Array, keep: list[int]) -> Array:
    q = rho.shape[0].bit_length() - 1
    if len(set(keep)) != len(keep) or any(k < 0 or k >= q for k in keep):
        raise ValueError("invalid retained wires")
    drop = [k for k in range(q) if k not in keep]
    axes = keep + drop + [k + q for k in keep] + [k + q for k in drop]
    a = rho.reshape([2] * (2 * q)).transpose(axes)
    dk, dd = 1 << len(keep), 1 << len(drop)
    return np.einsum("aibi->ab", a.reshape(dk, dd, dk, dd))


def trace_distance(a: Array, b: Array) -> float:
    z = (a - b + (a - b).conj().T) / 2
    return float(np.abs(np.linalg.eigvalsh(z)).sum() / 2)


def psd_sqrt(a: Array) -> Array:
    val, vec = np.linalg.eigh((a + a.conj().T) / 2)
    if val.min() < -1e-8:
        raise ValueError("matrix is not positive semidefinite")
    return (vec * np.sqrt(np.maximum(val, 0))) @ vec.conj().T


def fidelity(a: Array, b: Array) -> float:
    h = psd_sqrt(a)
    x = h @ b @ h
    val = np.linalg.eigvalsh((x + x.conj().T) / 2)
    return float(np.minimum(1.0, np.maximum(val, 0).__pow__(0.5).sum() ** 2))


def sweep(block: Array, fresh: Array, eta: float) -> Array:
    """One fresh ancilla, wire 0, traverses every block wire and is traced out."""
    k = block.shape[0].bit_length() - 1
    rho = tensor(fresh, block)
    for j in range(k):
        rho = collide(rho, 0, j + 1, eta)
    return reduce(rho, list(range(1, k + 1)))


def outputs_by_columns(a: Array, b: Array, n: int, m: int, eta: float) -> Array:
    block = power(a, n)
    for _ in range(m):
        block = sweep(block, b, eta)
    return block


def full_grid(a: Array, b: Array, n: int, m: int, eta: float, order: str = "rows") -> Array:
    if n + m > 11:
        raise ValueError("full_grid limited to 11 total qubits to avoid excessive RAM")
    rho = tensor(power(a, n), power(b, m))
    if order == "rows":
        gates = [(i, n + j) for i in range(n) for j in range(m)]
    elif order == "columns":
        gates = [(i, n + j) for j in range(m) for i in range(n)]
    else:
        raise ValueError("order must be rows or columns")
    for i, j in gates:
        rho = collide(rho, i, j, eta)
    return rho


def reused_memory(a: Array, b: Array, n: int, m: int, eta: float):
    reservoir = power(b, m)
    out = None
    for _ in range(n):
        rho = tensor(a, reservoir)
        for j in range(m):
            rho = collide(rho, 0, j + 1, eta)
        out = reduce(rho, [0])
        reservoir = reduce(rho, list(range(1, m + 1)))
    return out, reservoir


def self_test() -> dict:
    rng = np.random.default_rng(20261002)

    def rand_density(d):
        x = rng.normal(size=(d, d)) + 1j * rng.normal(size=(d, d))
        z = x @ x.conj().T
        return z / np.trace(z)

    max_grid = max_columns = max_fidelity_violation = 0.0
    for n, m in [(1, 1), (1, 3), (2, 2), (2, 3), (3, 2)]:
        for eta in [0.01, 0.3, 0.8, np.pi / 2]:
            a, b = rand_density(2), rand_density(2)
            r = full_grid(a, b, n, m, eta, "rows")
            c = full_grid(a, b, n, m, eta, "columns")
            o = outputs_by_columns(a, b, n, m, eta)
            max_grid = max(max_grid, float(np.linalg.norm(r - c)))
            max_columns = max(max_columns, float(np.linalg.norm(o - reduce(r, list(range(n))))))
            mem = reduce(r, list(range(n, n + m)))
            lower = fidelity(a, b) ** n
            max_fidelity_violation = max(max_fidelity_violation, lower - fidelity(mem, power(b, m)))
    checks = []
    for a, b, label in [(P, TAU, "pure_to_mixed"), (TAU, P, "mixed_to_pure")]:
        for n in [1, 2, 3]:
            out = outputs_by_columns(a, b, n, n, np.pi / 2)
            assert trace_distance(out, power(b, n)) < 1e-10
        for m, n in [(2, 2), (2, 3), (3, 2), (3, 3), (4, 3)]:
            eta = 0.5
            out, mem = reused_memory(a, b, n, m, eta)
            checks.append(
                dict(
                    direction=label,
                    m=m,
                    n=n,
                    eta=eta,
                    current_trace_distance=trace_distance(out, b),
                    current_infidelity=1 - fidelity(out, b),
                    reservoir_fidelity=fidelity(mem, power(b, m)),
                )
            )
    assert max_grid < 1e-10 and max_columns < 1e-10
    assert max_fidelity_violation < 1e-7
    return dict(
        max_row_column_frobenius_error=max_grid,
        max_reduced_column_frobenius_error=max_columns,
        max_fidelity_lower_bound_violation=max_fidelity_violation,
        samples=checks,
    )


if __name__ == "__main__":
    result = self_test()
    path = Path(__file__).with_name("checks.json")
    path.write_text(json.dumps(result, indent=2))
    print(json.dumps(result, indent=2))
