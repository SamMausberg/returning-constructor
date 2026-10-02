#!/usr/bin/env python3
"""Exact finite homogenizer calculations and rigorous plotting enclosures.

No product approximation, random sampling, or floating-point linear algebra is
used in the calculation of the data. Rational matrices are propagated with
Gaussian-integer arithmetic. Small spectra are enclosed by exact real-root
isolation (SymPy); two-cell spectra use rational square-root enclosures.

Large powers of the 16-dimensional Pauli transfer matrix use a fixed-point
matrix enclosure. Every entry has an integer centre and a common rigorous
absolute radius, both scaled by 2**512. The rounding and matrix-product error
bound is documented in BallMatrix.__matmul__. Floats appear only on CSV export
for graphics. The JSON file retains rational endpoints.
"""

from __future__ import annotations

import csv
import json
import time
from dataclasses import dataclass
from fractions import Fraction as Q
from itertools import product
from math import ceil, isqrt
from pathlib import Path

import numpy as np
import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
BITS = 512
SCALE = 1 << BITS


def ZERO(shape):
    """Zero matrix of exact Python integers/rationals."""
    return np.zeros(shape, dtype=object)


def qi(x):
    return x if isinstance(x, Q) else Q(int(x))


@dataclass(frozen=True)
class Interval:
    lo: Q
    hi: Q

    def __post_init__(self):
        if self.lo > self.hi:
            raise ValueError("Reversed enclosure")

    @classmethod
    def point(cls, x):
        x = qi(x)
        return cls(x, x)

    def __add__(self, other):
        o = other if isinstance(other, Interval) else Interval.point(other)
        return Interval(self.lo + o.lo, self.hi + o.hi)

    __radd__ = __add__

    def __neg__(self):
        return Interval(-self.hi, -self.lo)

    def __sub__(self, o):
        return self + -as_iv(o)

    def __rsub__(self, o):
        return as_iv(o) + -self

    def __mul__(self, o):
        o = as_iv(o)
        a = [self.lo * o.lo, self.lo * o.hi, self.hi * o.lo, self.hi * o.hi]
        return Interval(min(a), max(a))

    __rmul__ = __mul__

    def __truediv__(self, o):
        o = as_iv(o)
        if o.lo <= 0 <= o.hi:
            raise ZeroDivisionError("Interval contains zero")
        return self * Interval(1 / o.hi, 1 / o.lo)

    def __rtruediv__(self, o):
        return as_iv(o) / self

    def square(self):
        if self.lo <= 0 <= self.hi:
            return Interval(Q(0), max(self.lo * self.lo, self.hi * self.hi))
        return Interval(
            min(self.lo * self.lo, self.hi * self.hi), max(self.lo * self.lo, self.hi * self.hi)
        )

    def sqrt(self):
        # Positivity of density spectra is an exact premise, so an enclosure
        # crossing zero may be intersected with the nonnegative real line.
        if self.hi < 0:
            raise ValueError("Negative radicand")
        lo = max(Q(0), self.lo)
        hi = self.hi
        a = isqrt(lo.numerator * SCALE * SCALE // lo.denominator)
        b = isqrt(hi.numerator * SCALE * SCALE // hi.denominator)
        upper = Q(b, SCALE) if Q(b * b, SCALE * SCALE) == hi else Q(b + 1, SCALE)
        return Interval(Q(a, SCALE), upper)

    def abs(self):
        if self.lo >= 0:
            return self
        if self.hi <= 0:
            return -self
        return Interval(Q(0), max(-self.lo, self.hi))

    def clip01(self):
        return Interval(max(Q(0), self.lo), min(Q(1), self.hi))

    def centre_float(self):
        return float((self.lo + self.hi) / 2)

    def json(self):
        return {"lo": str(self.lo), "hi": str(self.hi)}


def as_iv(x):
    return x if isinstance(x, Interval) else Interval.point(x)


@dataclass
class GM:
    """Matrix (real + i imag)/den with exact integer entries."""

    real: np.ndarray
    imag: np.ndarray
    den: int = 1

    def normalized(self):
        from math import gcd

        g = int(self.den)
        for x in list(self.real.flat) + list(self.imag.flat):
            g = gcd(g, int(x))
        if g > 1:
            return GM(self.real // g, self.imag // g, self.den // g)
        return self

    def adjoint(self):
        return GM(self.real.T, -self.imag.T, self.den)

    def __matmul__(self, o):
        return GM(
            self.real @ o.real - self.imag @ o.imag,
            self.real @ o.imag + self.imag @ o.real,
            self.den * o.den,
        ).normalized()

    def __add__(self, o):
        return GM(
            self.real * o.den + o.real * self.den,
            self.imag * o.den + o.imag * self.den,
            self.den * o.den,
        ).normalized()

    def frac(self, i, j):
        return Q(int(self.real[i, j]), self.den), Q(int(self.imag[i, j]), self.den)


def state(k, kind):
    d = 1 << k
    r = ZERO((d, d))
    im = ZERO((d, d))
    if kind == "p":
        r[0, 0] = 1
        den = 1
    elif kind == "tau":
        r = np.eye(d, dtype=object)
        den = d
    else:
        raise ValueError(kind)
    return GM(r, im, den)


def tensor(a, b):
    return GM(
        np.kron(a.real, b.real) - np.kron(a.imag, b.imag),
        np.kron(a.real, b.imag) + np.kron(a.imag, b.real),
        a.den * b.den,
    )


def permutation(n, i, j):
    idx = np.arange(1 << n, dtype=np.int64)
    b1 = n - 1 - i
    b2 = n - 1 - j
    diff = ((idx >> b1) ^ (idx >> b2)) & 1
    return idx ^ (diff << b1) ^ (diff << b2)


# weights (c^2,s^2,cs) have a common denominator.
def collide(x, perm, w=(1, 1, 1, 2)):
    a, b, t, d = w
    R = x.real
    I = x.imag
    rr = a * R + b * R[np.ix_(perm, perm)] - t * I[perm, :] + t * I[:, perm]
    ii = a * I + b * I[np.ix_(perm, perm)] + t * R[perm, :] - t * R[:, perm]
    return GM(rr, ii, x.den * d).normalized()


def sweep(block, fresh, w=(1, 1, 1, 2), postselect=False):
    """Fresh ancilla is the leading wire. Return updated block and ancilla."""
    k = (block.real.shape[0]).bit_length() - 1
    x = tensor(fresh, block)
    for j in range(k):
        x = collide(x, permutation(k + 1, 0, j + 1), w)
    d = 1 << k
    if postselect:
        out = GM(x.real[:d, :d], x.imag[:d, :d], x.den).normalized()
    else:
        out = GM(
            x.real[:d, :d] + x.real[d:, d:], x.imag[:d, :d] + x.imag[d:, d:], x.den
        ).normalized()
    r = ZERO((2, 2))
    im = ZERO((2, 2))
    for i, j in product(range(2), repeat=2):
        r[i, j] = np.trace(x.real[i * d : (i + 1) * d, j * d : (j + 1) * d])
        im[i, j] = np.trace(x.imag[i * d : (i + 1) * d, j * d : (j + 1) * d])
    return out, GM(r, im, x.den).normalized()


def partial_single(x, which):
    k = x.real.shape[0].bit_length() - 1
    d = 1 << k
    r = ZERO((2, 2))
    im = ZERO((2, 2))
    bit = k - 1 - which
    for a in range(d):
        for b in range(d):
            if (a & ~(1 << bit)) == (b & ~(1 << bit)):
                r[(a >> bit) & 1, (b >> bit) & 1] += x.real[a, b]
                im[(a >> bit) & 1, (b >> bit) & 1] += x.imag[a, b]
    return GM(r, im, x.den).normalized()


def diagonal_error(x, target):
    assert (
        x.imag[0, 0]
        == x.imag[1, 1]
        == x.real[0, 1]
        == x.real[1, 0]
        == x.imag[0, 1]
        == x.imag[1, 0]
        == 0
    )
    z = Q(int(x.real[0, 0] - x.real[1, 1]), x.den)
    if target == "tau":
        return (1 - as_iv(1 - z * z).sqrt()) / 2
    if target == "p":
        return as_iv((1 - z) / 2)
    raise ValueError(target)


def block_eigenvalues(x):
    """Certified eigenvalues for an excitation-preserving two-qubit state."""
    assert x.real.shape == (4, 4)
    for i in range(4):
        for j in range(4):
            if i.bit_count() != j.bit_count():
                assert x.real[i, j] == x.imag[i, j] == 0
    a = as_iv(x.frac(0, 0)[0])
    d = as_iv(x.frac(3, 3)[0])
    b = as_iv(x.frac(1, 1)[0])
    c = as_iv(x.frac(2, 2)[0])
    rr, ii = x.frac(1, 2)
    disc = ((b - c) / 2).square() + rr * rr + ii * ii
    mid = (b + c) / 2
    rad = disc.sqrt()
    return [a, d, mid - rad, mid + rad]


def fidelity_tau_two(x):
    return (sum((v.sqrt() for v in block_eigenvalues(x)), as_iv(0)) / 2).square().clip01()


def fidelity_tau_general(x):
    """Exact characteristic polynomials; rational isolation of all real roots."""
    d = len(x.real)
    k = d.bit_length() - 1
    total = as_iv(0)
    eps = sp.Rational(x.den, 1 << 100)
    for excitations in range(k + 1):
        ids = [i for i in range(d) if i.bit_count() == excitations]
        A = sp.Matrix(
            [[sp.Integer(x.real[i, j]) + sp.I * sp.Integer(x.imag[i, j]) for j in ids] for i in ids]
        )
        poly = A.charpoly().as_poly()
        poly = sp.Poly(poly.as_expr(), poly.gens[0], domain=sp.ZZ)
        assert all(c.is_Integer for c in poly.all_coeffs())
        count = 0
        for (lo, hi), mult in sp.polys.polytools.intervals(poly, eps=eps):
            lo = Q(int(lo.p), int(lo.q)) / x.den
            hi = Q(int(hi.p), int(hi.q)) / x.den
            if lo < 0 and hi < 0:
                raise AssertionError("Density operator lost positivity")
            total += mult * Interval(max(Q(0), lo), hi).sqrt()
            count += mult
        assert count == len(ids)
    return (total.square() / d).clip01()


# Hermitian Pauli basis, tensor order II,IX,IY,IZ,XI,...,ZZ.
I2 = GM(np.eye(2, dtype=object), ZERO((2, 2)))
X = GM(np.array([[0, 1], [1, 0]], dtype=object), ZERO((2, 2)))
Y = GM(ZERO((2, 2)), np.array([[0, -1], [1, 0]], dtype=object))
Z = GM(np.array([[1, 0], [0, -1]], dtype=object), ZERO((2, 2)))
PAULI = [tensor(a, b) for a, b in product([I2, X, Y, Z], repeat=2)]


def pauli_transfer(w):
    vals = np.empty((16, 16), dtype=object)
    for j, P in enumerate(PAULI):
        z, _ = sweep(P, state(1, "tau"), w)
        for i, Qp in enumerate(PAULI):
            re = np.trace(Qp.real @ z.real - Qp.imag @ z.imag)
            im = np.trace(Qp.real @ z.imag + Qp.imag @ z.real)
            assert im == 0
            vals[i, j] = Q(int(re), 4 * z.den)
    assert list(vals[0, :]) == [Q(1)] + [Q(0)] * 15
    return vals


@dataclass
class BallMatrix:
    mid: np.ndarray
    rad: int

    @classmethod
    def fractions(cls, a):
        z = ZERO(a.shape)
        for ij in np.ndindex(a.shape):
            q = qi(a[ij])
            z[ij] = q.numerator * SCALE // q.denominator
        return cls(z, 1)

    @classmethod
    def identity(cls, n):
        return cls(np.eye(n, dtype=object) * SCALE, 0)

    def __matmul__(self, o):
        assert self.mid.shape[1] == o.mid.shape[0]
        d = self.mid.shape[1]
        aa = max(abs(int(x)) for x in self.mid.flat)
        bb = max(abs(int(x)) for x in o.mid.flat)
        raw = self.mid @ o.mid
        # Each summand's error is <= aa*rB + bb*rA + rA*rB
        # in units SCALE**2. Divide with outward rounding and add one unit
        # for flooring the central product. This encloses every real entry.
        numer = d * (aa * o.rad + bb * self.rad + self.rad * o.rad)
        rad = (numer + SCALE - 1) // SCALE + 1
        return BallMatrix(raw // SCALE, int(rad))

    def power(self, n):
        if n < 0:
            raise ValueError("negative exponent")
        a = self
        out = BallMatrix.identity(len(self.mid))
        while n:
            if n & 1:
                out = out @ a
            n //= 2
            if n:
                a = a @ a
        return out

    def entry(self, i, j):
        return Interval(
            Q(int(self.mid[i, j]) - self.rad, SCALE), Q(int(self.mid[i, j]) + self.rad, SCALE)
        )


def pauli_initial():
    x = np.zeros((16, 1), dtype=object)
    for i in [0, 3, 12, 15]:
        x[i, 0] = 1
    return BallMatrix(x * SCALE, 0)


def state_entries(vec):
    r = [[as_iv(0) for j in range(4)] for i in range(4)]
    im = [[as_iv(0) for j in range(4)] for i in range(4)]
    for v, P in enumerate(PAULI):
        for i, j in product(range(4), repeat=2):
            if P.real[i, j]:
                r[i][j] += vec.entry(v, 0) * int(P.real[i, j]) / 4
            if P.imag[i, j]:
                im[i][j] += vec.entry(v, 0) * int(P.imag[i, j]) / 4
    return r, im


def certified_two_errors(vec):
    r, im = state_entries(vec)
    mid = (r[1][1] + r[2][2]) / 2
    rad = (((r[1][1] - r[2][2]) / 2).square() + r[1][2].square() + im[1][2].square()).sqrt()
    eigs = [r[0][0], r[3][3], mid - rad, mid + rad]
    joint = sum(((x - Q(1, 4)).abs() for x in eigs), as_iv(0)) / 2
    m1 = vec.entry(12, 0).abs() / 2
    m2 = vec.entry(3, 0).abs() / 2
    marginal = Interval(max(m1.lo, m2.lo), max(m1.hi, m2.hi))
    singlet = (r[1][1] + r[2][2] - 2 * r[1][2]) / 2
    return marginal, joint, singlet


def scalar_power(q, n):
    arr = np.array([[qi(q)]], dtype=object)
    return BallMatrix.fractions(arr).power(n).entry(0, 0)


def dump_records(path, rows):
    if not rows:
        raise ValueError("empty data")
    serial = []
    flat = []
    for row in rows:
        s = {}
        f = {}
        for key, v in row.items():
            if isinstance(v, Interval):
                s[key] = v.json()
                f[key + "_lo"] = float(v.lo)
                f[key + "_hi"] = float(v.hi)
                f[key] = v.centre_float()
            elif isinstance(v, Q):
                s[key] = str(v)
                f[key] = float(v)
            else:
                s[key] = v
                f[key] = v
        serial.append(s)
        flat.append(f)
    path.with_suffix(".json").write_text(json.dumps(serial, indent=2))
    keys = list(dict.fromkeys(k for r in flat for k in r))
    with path.with_suffix(".csv").open("w", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=keys)
        writer.writeheader()
        writer.writerows(flat)


def make_ratio_data():
    # Reservoir scan: n=2. Output state is four dimensional for every M.
    op = state(2, "p")
    om = state(2, "tau")
    conditional = state(2, "tau")
    rows = []
    for M in range(1, 33):
        op, _ = sweep(op, state(1, "tau"))
        om, _ = sweep(om, state(1, "p"))
        conditional, _ = sweep(conditional, state(1, "p"), postselect=True)
        ep = diagonal_error(partial_single(op, 1), "tau")
        em = diagonal_error(partial_single(om, 1), "p")
        fm = as_iv(Q(int(np.trace(conditional.real)), conditional.den))
        row = {
            "M": M,
            "n": 2,
            "forward_bound": Interval(ep.lo, 4 * ep.hi),
            "forward_error": ep,
            "reverse_ratio": em / fm,
        }
        if M <= 6:
            reservoir = state(M, "tau")
            for t in range(2):
                reservoir, out = sweep(reservoir, state(1, "p"))
            assert diagonal_error(out, "tau") == ep
            fp = fidelity_tau_general(reservoir)
            ratio = ep / fp
            assert ratio.lo >= ep.lo and ratio.hi <= 4 * ep.hi
            row["forward_ratio"] = ratio
        rows.append(row)
    dump_records(ROOT / "data/ratio_reservoir", rows)
    # Use scan: M=2. Both memory spectra reduce to 1x1,2x2,1x1 blocks.
    rp = state(2, "tau")
    rm = state(2, "p")
    rows = []
    for n in range(1, 49):
        rp, op = sweep(rp, state(1, "p"))
        rm, om = sweep(rm, state(1, "tau"))
        ep = diagonal_error(op, "tau")
        em = diagonal_error(om, "p")
        fp = fidelity_tau_two(rp)
        fm = as_iv(rm.frac(0, 0)[0])
        rows.append(
            {"M": 2, "n": n, "forward_ratio": ep / fp, "reverse_ratio": em / fm, "limit": 2}
        )
    dump_records(ROOT / "data/ratio_uses", rows)
    return {"reservoir_points": 32, "forward_certified_points": 6, "uses_points": 48}


def make_singlet_data():
    rows = []
    largest = Q(0)
    for k in range(1, 33):
        # Rational c,s parametrization with M=k**5 and eta=2 arctan(1/(2k**2)).
        # Thus eta~M**(-2/5), inside the proved collective-limit range.
        d = 4 * k**4 + 1
        c = 4 * k**4 - 1
        s = 4 * k * k
        M = k**5
        w = (c * c, s * s, c * s, d * d)
        u = Q(s * s, d * d)
        trans = BallMatrix.fractions(pauli_transfer(w))
        vec = trans.power(M) @ pauli_initial()
        marg, joint, singlet = certified_two_errors(vec)
        witness = scalar_power(1 - u * u, M) / 4
        assert singlet.lo <= Q(1, 4) - witness.lo and singlet.hi >= Q(1, 4) - witness.hi
        bound = u * as_iv(M).sqrt()  # n=2 telescope upper bound, clipped at 1.
        bound = Interval(min(Q(1), bound.lo), min(Q(1), bound.hi))
        largest = max(largest, marg.hi - marg.lo, joint.hi - joint.lo)
        rows.append(
            {
                "k": k,
                "M": M,
                "u": u,
                "c": Q(c, d),
                "s": Q(s, d),
                "marginal_error": marg,
                "joint_error": joint,
                "singlet_lower": witness,
                "prefix_return_upper": bound,
            }
        )
    dump_records(ROOT / "data/singlet", rows)
    return {
        "points": len(rows),
        "bits": BITS,
        "maximum_trace_error_enclosure_width": str(largest),
        "width_float": float(largest),
    }


def spectral_lower(n, delta):
    # Proposition 5, exact integer rounding. Also include rank m>=n.
    return max(n, n * (ceil(1 / delta - Q(1, (1 << n) - 1)) - 1) + 1)


def dyadic_upper(n, delta):
    assert 0 < delta <= Q(1, 4)
    return ceil(Q(n) / delta) + n - 2


def make_cost_data():
    rows = []
    for n in range(1, 129):
        d = Q(1, 10)
        rows.append(
            {
                "n": n,
                "delta": d,
                "forward": ceil(Q(n, 2)),
                "reverse_lower": spectral_lower(n, d),
                "reverse_upper": dyadic_upper(n, d),
            }
        )
    dump_records(ROOT / "data/cost_horizon", rows)
    rows = []
    for t in range(4, 101):
        d = Q(1, t)
        n = 32
        rows.append(
            {
                "n": n,
                "delta": d,
                "inverse_delta": t,
                "forward": ceil(Q(n, 2)),
                "reverse_lower": spectral_lower(n, d),
                "reverse_upper": dyadic_upper(n, d),
            }
        )
    dump_records(ROOT / "data/cost_tolerance", rows)
    return {"horizon_points": 128, "tolerance_points": 97}


def consistency_checks():
    checks = []
    # Rational row-column equality, comparing exact integer fractions.
    for M in range(1, 5):
        for n in [1, 2, 3]:
            for a, b in [("p", "tau"), ("tau", "p")]:
                mem = state(M, b)
                last = None
                for _ in range(n):
                    mem, last = sweep(mem, state(1, a))
                block = state(n, a)
                for _ in range(M):
                    block, _ = sweep(block, state(1, b))
                last2 = partial_single(block, n - 1)
                assert (
                    last.real.tolist() == last2.real.tolist()
                    and last.imag.tolist() == last2.imag.tolist()
                    and last.den == last2.den
                )
                checks.append((M, n, a, b))
    # Test the generic fixed-point enclosure against exact rational transfer.
    w = (9, 16, 12, 25)
    T = pauli_transfer(w)
    bm = BallMatrix.fractions(T)
    exact = np.array(
        [Q(1) if i in [0, 3, 12, 15] else Q(0) for i in range(16)], dtype=object
    ).reshape(16, 1)
    for power in [1, 2, 5, 9]:
        ee = exact.copy()
        for _ in range(power):
            ee = T @ ee
        vv = bm.power(power) @ pauli_initial()
        assert all(vv.entry(i, 0).lo <= ee[i, 0] <= vv.entry(i, 0).hi for i in range(16))
    return {"exact_row_column_cases": len(checks), "interval_vs_exact_powers": [1, 2, 5, 9]}


def main():
    (ROOT / "data").mkdir(exist_ok=True)
    start = time.time()
    report = {
        "precision_bits": BITS,
        "claim": "Rigorous finite-size enclosures; not proofs of asymptotic statements.",
    }
    report["consistency"] = consistency_checks()
    print("Exact consistency checks passed.", flush=True)
    report["ratios"] = make_ratio_data()
    print("Ratio data finished.", flush=True)
    report["singlet"] = make_singlet_data()
    print("Singlet data finished.", flush=True)
    report["costs"] = make_cost_data()
    report["elapsed_seconds"] = time.time() - start
    (ROOT / "verification/figure_certificates.json").write_text(json.dumps(report, indent=2))
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
