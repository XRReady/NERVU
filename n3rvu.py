"""
Attribution License (Free Use, Derivatives Allowed)

Copyright (c) 2025 Dale Spencer

Permission is hereby granted, free of charge, to any person obtaining a copy of
this work and associated materials (the "Work"), to use, reproduce, modify,
publish, distribute, sublicense, and/or sell copies of the Work, for any purpose,
commercial or non-commercial, provided that the following attribution requirement
is met:

    The original author must be credited in any public or distributed use of
    this Work or derivative works, and this attribution notice must be included
    with all copies or substantial portions of the Work.

No restrictions are placed on licensing of derivative works beyond preservation
of the attribution notice. This license imposes no copyleft requirements.

THE WORK IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED.

Author:
    Dale Spencer
    Creator of the NERVU Algebra (N3RVU baric triple formulation)
    TexasDataSafe | texasdatasafe@gmail.com
"""


from dataclasses import dataclass
from typing import Iterable, Tuple, Union

Number = Union[int, float]

class _ZeroLiteral:
    def __repr__(self) -> str:
        return "ZERO"

class _InvZeroLiteral:
    def __repr__(self) -> str:
        return "INV_ZERO"

ZERO = _ZeroLiteral()
INV_ZERO = _InvZeroLiteral()


@dataclass(frozen=True, slots=True)
class N3RVU:
    r: float
    v: float
    u: float

    @staticmethod
    def from_iterable(x: Iterable[Number]) -> "N3RVU":
        r, v, u = x
        return N3RVU(float(r), float(v), float(u))

    def as_tuple(self) -> Tuple[float, float, float]:
        return (self.r, self.v, self.u)

    def __repr__(self) -> str:
        return f"({self.r},{self.v},{self.u})"

    def __neg__(self) -> "N3RVU":
        return N3RVU(-self.r, -self.v, -self.u)

    # ---- Addition / Subtraction ----
    def __add__(self, other: "N3RVU") -> "N3RVU":
        if not isinstance(other, N3RVU):
            return NotImplemented
        return N3RVU(self.r + other.r, self.v + other.v, self.u + other.u)

    def __radd__(self, other: "N3RVU") -> "N3RVU":
        return self.__add__(other)

    def __sub__(self, other: "N3RVU") -> "N3RVU":
        if not isinstance(other, N3RVU):
            return NotImplemented
        return N3RVU(self.r - other.r, self.v - other.v, self.u - other.u)

    def __rsub__(self, other: "N3RVU") -> "N3RVU":
        if not isinstance(other, N3RVU):
            return NotImplemented
        return N3RVU(other.r - self.r, other.v - self.v, other.u - self.u)

    # ---- Scalar Multiplication ----
    def _scalar_mul(self, k: Number) -> "N3RVU":
        return N3RVU(k * self.r, k * self.v, k * self.u)

    # ---- RVU Multiplication ----
    def __mul__(self, other):
        # numeric scalars, but treat 0 specially to enforce RV identity
        if isinstance(other, (int, float)):
            if other == 0:
                # "pure zero" effect for numeric 0:
                # a * 0 = (a2, a0 + a1, 0)
                return N3RVU(self.u, self.r + self.v, 0.0)
            return self._scalar_mul(other)

        if isinstance(other, _ZeroLiteral):
            # ZERO sentinel: same mapping as numeric 0
            return N3RVU(self.u, self.r + self.v, 0.0)

        if isinstance(other, _ZeroLiteral):
            return N3RVU(self.u, self.r + self.v, 0.0)

        if not isinstance(other, N3RVU):
            return NotImplemented

        a0, a1, a2 = self.r, self.v, self.u
        b0, b1, b2 = other.r, other.v, other.u

        c0 = a0*b0 + a2*b1 + a1*b2
        c1 = a0*b1 + a1*b0 + a1*b1
        c2 = a0*b2 + a2*b0 + a2*b2

        return N3RVU(c0, c1, c2)
    
    def __rmul__(self, other):
        if isinstance(other, (int, float)):
            if other == 0:
                return N3RVU(self.u, self.r + self.v, 0.0)
            return self._scalar_mul(other)
        if isinstance(other, _ZeroLiteral):
            return self.__mul__(ZERO)
        return NotImplemented

    # ---- Division (algebraically correct right-division) ----
    @staticmethod
    def _det_right_map(b0: float, b1: float, b2: float) -> float:
        # det M(b) = (b0^2 - b1*b2) * (b0 + b1 + b2)
        return (b0*b0 - b1*b2) * (b0 + b1 + b2)

    @staticmethod
    def _solve3(M, a, eps: float = 1e-12):
        """
        Solve M x = a for 3x3 M using Gaussian elimination with partial pivoting.
        Returns a 3-tuple or None if the system is singular/ill-conditioned.
        """
        # Build augmented matrix [M | a]
        A = [
            [float(M[0][0]), float(M[0][1]), float(M[0][2]), float(a[0])],
            [float(M[1][0]), float(M[1][1]), float(M[1][2]), float(a[1])],
            [float(M[2][0]), float(M[2][1]), float(M[2][2]), float(a[2])],
        ]

        # Forward elimination with partial pivoting
        for col in range(3):
            # Find pivot
            pivot_row = max(range(col, 3), key=lambda r: abs(A[r][col]))
            if abs(A[pivot_row][col]) < eps:
                return None  # singular / ill-conditioned

            # Swap if needed
            if pivot_row != col:
                A[col], A[pivot_row] = A[pivot_row], A[col]

            # Eliminate below
            piv = A[col][col]
            for r in range(col + 1, 3):
                factor = A[r][col] / piv
                for c in range(col, 4):
                    A[r][c] -= factor * A[col][c]

        # Back substitution
        x = [0.0, 0.0, 0.0]
        for i in reversed(range(3)):
            piv = A[i][i]
            if abs(piv) < eps:
                return None
            s = A[i][3] - sum(A[i][j] * x[j] for j in range(i + 1, 3))
            x[i] = s / piv

        return tuple(x)

    def __truediv__(self, other):
        # ---- scalar division (keep the totalized rule for scalar zero) ----
        if isinstance(other, (int, float)):
            if other == 0:
                # RVU rule: a / 0 = (a1, 0, a0 + a2)
                return N3RVU(self.v, 0.0, self.r + self.u)
            return self._scalar_mul(1.0 / float(other))

        # ---- special ZERO/INV_ZERO sentinels (preserve semantics) ----
        if isinstance(other, (_ZeroLiteral, _InvZeroLiteral)):
            return N3RVU(self.v, 0.0, self.r + self.u)

        # ---- algebraic right-division by an N3RVU element ----
        if not isinstance(other, N3RVU):
            return NotImplemented

        a0, a1, a2 = self.r, self.v, self.u
        b0, b1, b2 = other.r, other.v, other.u

        # Right-multiplication matrix M(b) from product x*b:
        # y0 = x0*b0 + x1*b2 + x2*b1
        # y1 = x0*b1 + x1*(b0+b1) + x2*0
        # y2 = x0*b2 + x1*0 + x2*(b0+b2)
        M = [
            [b0,     b2,       b1],
            [b1, b0 + b1,     0.0],
            [b2,     0.0,  b0 + b2],
        ]
        a_vec = [a0, a1, a2]

        # Invertibility check (fast, exact formula)
        det = self._det_right_map(b0, b1, b2)
        if abs(det) < 1e-12:
            raise ZeroDivisionError(
                "Right-division undefined: divisor is non-invertible "
                "(det M(b) = 0 where det = (b0^2 - b1*b2)*(b0 + b1 + b2))."
            )

        x = self._solve3(M, a_vec)
        if x is None:
            # Extremely ill-conditioned numerics fallback (should be rare since det≠0),
            # but we keep the contract clear.
            raise ZeroDivisionError("Right-division failed due to numerical ill-conditioning.")

        return N3RVU(*x)


    def __rtruediv__(self, other):
        return NotImplemented

    # ---- Equality ----
    def __eq__(self, other) -> bool:
        if not isinstance(other, N3RVU):
            return False
        return (self.r, self.v, self.u) == (other.r, other.v, other.u)

    # ---- Basis Constructors ----
    @staticmethod
    def e0(x: Number = 1.0) -> "N3RVU":
        return N3RVU(float(x), 0.0, 0.0)

    @staticmethod
    def e1(x: Number = 1.0) -> "N3RVU":
        return N3RVU(0.0, float(x), 0.0)

    @staticmethod
    def e2(x: Number = 1.0) -> "N3RVU":
        return N3RVU(0.0, 0.0, float(x))

    # ---- Weight ----
    def weight(self) -> float:
        return self.r + self.v + self.u

    # ---- Literal Operator Helpers ----
    def Z(self) -> "N3RVU":
        return self * ZERO

    def U(self) -> "N3RVU":
        return self / ZERO


def _demo() -> None:
    e0, e1, e2 = N3RVU.e0(), N3RVU.e1(), N3RVU.e2()
    assert e1 * e2 == e2 * e1
    a = N3RVU(1, 2, 3)
    b = N3RVU(2, 1, 1)
    print("a = ", a)
    print("b = ", b)
    print("a * b = ", a * b)    # Bilinear product
    print("a / b = ", a / b)    # Solve for M(b),x = a
    print("a / 2 = ", a / 2)    # Scalar division
    print("a * 2 = ", a * 2)    # Scalar multiplication
    print("a * 0 = ", a * 0)    # Multiplication by zero (transforms reals to the V-space and U's to real-space)
    print("0 * a = ", a * 0)    # Multiplication by zero (transforms reals to the V-space and U's to real-space)
    print("a / 0 = ", a / 0)    # Division by zero (transforms reals to the U-space and V's to real-space)
    print("a.Z() = ", a.Z())    # Multiply by pure zero
    print("a.U() = ", a.U())    # Divide by pure zero

    a = N3RVU(1, 0, 0)
    print()
    print("a = ", a)
    print("(a * 0) / 0 = ", (a * 0) / 0)

    #basic case that shows non-associativty
    x = N3RVU(0,1,0)
    y = N3RVU(0,1,0)
    z = N3RVU(0,0,1)

    xy   = x * y
    yz   = y * z
    xy_z = xy * z
    x_yz = x * yz

    print("x     =", x)
    print("y     =", y)
    print("z     =", z)
    print()

    print("xy    = x * y   =", xy)
    print("yz    = y * z   =", yz)
    print()

    print("(xy)z = (x*y)*z =", xy_z)
    print("x(yz) = x*(y*z) =", x_yz)
    print()

    if xy_z == x_yz:
        print("Associative.")
    else:
        print("NOT associative.")

if __name__ == "__main__":
    _demo()
