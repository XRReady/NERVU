# rvu_checks.py — RVU algebra diagnostics (baric + Jordan witness)
# (Formatted printing version)

from __future__ import annotations
import itertools
from typing import Optional, Dict, Tuple, List
import sympy as sp

E = [sp.Matrix([1,0,0]), sp.Matrix([0,1,0]), sp.Matrix([0,0,1])]
def e(i:int): return E[i]

# ----------------------------
# Structure constants and core
# ----------------------------
C = None
def _to_sympy_3x3x3(obj) -> sp.MutableDenseNDimArray:
    A = sp.MutableDenseNDimArray.zeros(3,3,3)
    for i in range(3):
        for j in range(3):
            for k in range(3):
                A[i,j,k] = sp.sympify(obj[i][j][k])
    return A

def set_structure_constants(tensor_3x3x3):
    global C
    C = _to_sympy_3x3x3(tensor_3x3x3)

def mul(u: sp.Matrix, v: sp.Matrix) -> sp.Matrix:
    assert C is not None
    u = sp.Matrix(u); v = sp.Matrix(v)
    out = sp.zeros(3,1)
    for k in range(3):
        acc = 0
        for i in range(3):
            for j in range(3):
                acc += u[i]*v[j]*C[i,j,k]
        out[k] = sp.simplify(acc)
    return out

def left_matrix(x: sp.Matrix) -> sp.Matrix:
    assert C is not None
    x = sp.Matrix(x)
    L = sp.zeros(3)
    for j in range(3):
        col = sp.zeros(3,1)
        for k in range(3):
            s = 0
            for i in range(3):
                s += x[i]*C[i,j,k]
            col[k] = sp.simplify(s)
        L[:,j] = col
    return L

def assoc(u: sp.Matrix, v: sp.Matrix, w: sp.Matrix) -> sp.Matrix:
    return sp.simplify(mul(mul(u,v), w) - mul(u, mul(v,w)))

# -------------
# RVU algebra
# -------------
C_rvu = [[[0,0,0] for _ in range(3)] for __ in range(3)]
C_rvu[0][0][0] = 1
C_rvu[2][1][0] = 1
C_rvu[1][2][0] = 1
C_rvu[0][1][1] = 1
C_rvu[1][0][1] = 1
C_rvu[1][1][1] = 1
C_rvu[0][2][2] = 1
C_rvu[2][0][2] = 1
C_rvu[2][2][2] = 1
set_structure_constants(C_rvu)

# -----------------
# Property checkers
# -----------------
def is_unital() -> bool:
    for i in range(3):
        if mul(e(0), e(i)) != e(i): return False
        if mul(e(i), e(0)) != e(i): return False
    return True

def is_associative() -> bool:
    for i,j,k in itertools.product(range(3), repeat=3):
        if assoc(e(i), e(j), e(k)) != sp.zeros(3,1):
            return False
    return True

def is_commutative() -> bool:
    for i,j in itertools.product(range(3), repeat=2):
        if mul(e(i), e(j)) != mul(e(j), e(i)):
            return False
    return True

def check_flexibility() -> bool:
    for i,j in itertools.product(range(3), repeat=2):
        if mul(e(i), mul(e(j), e(i))) != mul(mul(e(i), e(j)), e(i)):
            return False
    return True

def check_alternativity() -> bool:
    for i,j in itertools.product(range(3), repeat=2):
        if mul(mul(e(i), e(i)), e(j)) != mul(e(i), mul(e(i), e(j))): return False
        if mul(e(j), mul(e(i), e(i))) != mul(mul(e(j), e(i)), e(i)): return False
    return True

def check_power_associativity() -> bool:
    a,b,c = sp.symbols('a b c')
    x = a*e(0) + b*e(1) + c*e(2)
    return assoc(x,x,x) == sp.zeros(3,1)

def bracket(u: sp.Matrix, v: sp.Matrix) -> sp.Matrix:
    return mul(u,v) - mul(v,u)

def jacobi_identity_holds() -> bool:
    for i,j,k in itertools.product(range(3), repeat=3):
        x,y,z = e(i), e(j), e(k)
        J = bracket(x, bracket(y,z)) + bracket(y, bracket(z,x)) + bracket(z, bracket(x,y))
        if J != sp.zeros(3,1): return False
    return True

# ----------------
# Jordan tests
# ----------------
def check_jordan_identity_symbolic() -> bool:
    a,b,c,d,e_,f = sp.symbols('a b c d e f')
    x = a*e(0) + b*e(1) + c*e(2)
    y = d*e(0) + e_*e(1) + f*e(2)
    x2 = mul(x,x)
    return sp.simplify(mul(mul(x2, y), x) - mul(x2, mul(y, x))) == sp.zeros(3,1)

def jordan_counterexample():
    x = sp.Matrix([1,1,-1])
    y = e(1)
    diff = mul(mul(mul(x,x),y), x) - mul(mul(x,x), mul(y,x))
    return x, y, diff

# ----------------
# Baric
# ----------------
def weight(x: sp.Matrix):
    x = sp.Matrix(x)
    return sp.simplify(x[0] + x[1] + x[2])

def check_baric() -> bool:
    x0,x1,x2,y0,y1,y2 = sp.symbols('x0 x1 x2 y0 y1 y2')
    x = sp.Matrix([x0,x1,x2]); y = sp.Matrix([y0,y1,y2])
    lhs = weight(mul(x,y))
    rhs = sp.simplify(weight(x)*weight(y))
    return sp.simplify(lhs - rhs) == 0

# -------------------------
# Center / Idempotents / Units
# -------------------------
def center():
    z0,z1,z2 = sp.symbols('z0 z1 z2')
    z = sp.Matrix([z0,z1,z2])
    eqs = []
    for i in range(3):
        eqs += list(mul(z, e(i)) - mul(e(i), z))
    A, _ = sp.linear_eq_to_matrix(eqs, [z0,z1,z2])
    Ns = A.nullspace()
    return Ns if Ns else [sp.zeros(3,1)]

def idempotents():
    x0,x1,x2 = sp.symbols('x0 x1 x2')
    x = sp.Matrix([x0,x1,x2])
    return sp.solve(list(mul(x,x) - x), [x0,x1,x2], dict=True)

def det_left_matrix_symbolic() -> sp.Expr:
    x0,x1,x2 = sp.symbols('x0 x1 x2')
    x = sp.Matrix([x0,x1,x2])
    return sp.simplify(left_matrix(x).det())

def has_zero_divisors_sample(samples: Optional[int] = 10) -> bool:
    rng_vals = [-2,-1,1,2]; checked = 0
    for a in rng_vals:
        for b in rng_vals:
            for c in rng_vals:
                if a==0 and b==0 and c==0: continue
                x = a*e(0) + b*e(1) + c*e(2)
                if left_matrix(x).det() == 0: return True
                checked += 1
                if samples and checked >= samples: return False
    return False

def units_sample(samples: Optional[int] = 10):
    rng_vals = [-2,-1,1,2]; res = []; checked = 0
    for a in rng_vals:
        for b in rng_vals:
            for c in rng_vals:
                if a==0 and b==0 and c==0: continue
                x = a*e(0) + b*e(1) + c*e(2)
                if left_matrix(x).det() != 0:
                    res.append(sp.Matrix([a,b,c]))
                checked += 1
                if samples and checked >= samples: return res
    return res

# -------------------------
# Literal zero ops
# -------------------------
def Z(a):
    a = sp.Matrix(a); a0,a1,a2 = a
    return sp.Matrix([a2, a0 + a1, 0])

def U(a):
    a = sp.Matrix(a); a0,a1,a2 = a
    return sp.Matrix([a1, 0, a0 + a2])

# ------- Report As Dict -------
def report_all():
    out: Dict[str, object] = {}
    out['unital(e0)'] = is_unital()
    out['associative'] = is_associative()
    out['commutative'] = is_commutative()
    out['flexible'] = check_flexibility()
    out['alternative'] = check_alternativity()
    out['power_associative'] = check_power_associativity()
    out['jordan_identity'] = check_jordan_identity_symbolic()
    out['jacobi(commutator)'] = jacobi_identity_holds()
    out['baric(weight multiplicative)'] = check_baric()
    out['zero_divisors_sample'] = has_zero_divisors_sample()

    out['L_e0'], out['L_e1'], out['L_e2'] = left_matrix(e(0)), left_matrix(e(1)), left_matrix(e(2))
    out['center_basis'] = center()
    out['idempotents'] = idempotents()
    out['units_sample'] = units_sample()

    detL = det_left_matrix_symbolic()
    out['det_L_(x0,x1,x2)'] = detL
    out['det_L factorization'] = sp.factor(detL)

    x,y,diff = jordan_counterexample()
    out['jordan_counterexample'] = {'x': x, 'y': y, 'diff': diff}

    return out

# -------------------------
# Pretty Printer
# -------------------------
def _format_matrix(M: sp.Matrix, indent: str = "    ") -> str:
    return "\n".join(f"{indent}[ " + " ".join(str(M[i,j]) for j in range(M.cols)) + " ]"
                     for i in range(M.rows))

def _vec_tuple(M):
    M = sp.Matrix(M)
    return tuple(M[i,0] for i in range(M.rows))

def _is_identity(M):
    return M == sp.eye(M.rows)

def print_formatted_report(res):
    print("=== RVU Algebraic Property Summary ===\n")
    print("Core Structural Properties:")
    print(f"    • Unital (e₀ is identity):        {res['unital(e0)']}")
    print(f"    • Associative:                    {res['associative']}")
    print(f"    • Commutative:                    {res['commutative']}")
    print(f"    • Flexible:                       {res['flexible']}")
    print(f"    • Alternative:                    {res['alternative']}")
    print(f"    • Power-associative:              {res['power_associative']}\n")

    print("Jordan / Lie Structure:")
    print(f"    • Jordan identity holds:          {res['jordan_identity']}")
    print(f"    • Jacobi identity (commutator):   {res['jacobi(commutator)']}\n")

    print("Baric Structure:")
    print(f"    • Weight map multiplicative:      {res['baric(weight multiplicative)']}\n")

    print("Zero Divisors:")
    print(f"    • Zero divisors exist (sampled):  {res['zero_divisors_sample']}\n")

    print("\n=== Left Multiplication Operators L_{eᵢ} ===\n")
    for name in ['L_e0','L_e1','L_e2']:
        print(f"{name.replace('L_','L(e')} ) =")
        print(_format_matrix(res[name]))
        print()

    print("\n=== Center of the Algebra ===")
    center_basis = res['center_basis']
    if len(center_basis) == 3 and all(v in center_basis for v in E):
        print("    Basis: e₀, e₁, e₂   (the center is the whole algebra)")
    else:
        print("    Basis:")
        for v in center_basis:
            print(f"        {_vec_tuple(v)}")

    print("\n\n=== Idempotents (solutions to x² = x) ===")
    for sol in res['idempotents']:
        print(f"    ({sol[sp.Symbol('x0')]}, {sol[sp.Symbol('x1')]}, {sol[sp.Symbol('x2')]})")

    print("\n\n=== Sample Units (non-zero elements with det(Lₓ) ≠ 0) ===")
    units = res['units_sample']
    line = []
    for u in units:
        line.append(str(_vec_tuple(u)))
        if len(line) == 4:
            print("    " + ", ".join(line))
            line = []
    if line:
        print("    " + ", ".join(line))

    print(f"\n\n=== Determinant of Left Multiplication ===")
    print(f"    det(Lₓ) = {res['det_L_(x0,x1,x2)']}")
    print(f"    Factored: {res['det_L factorization']}")

    x,y,diff = res['jordan_counterexample']['x'], res['jordan_counterexample']['y'], res['jordan_counterexample']['diff']
    print("\n\n=== Jordan Identity Counterexample ===")
    print(f"    x = {_vec_tuple(x)}")
    print(f"    y = {_vec_tuple(y)}")
    print(f"    ( (x²·y)·x ) − ( x²·(y·x) ) = {_vec_tuple(diff)}")

    print("\n\n=== Literal Zero Operators (Spec Check) ===")
    print(f"    Z(e₁) = {_vec_tuple(Z(e(1)))}")
    print(f"    U(e₂) = {_vec_tuple(U(e(2)))}")

if __name__ == "__main__":
    res = report_all()
    print_formatted_report(res)
