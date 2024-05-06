# SageMath imports
from sage.all import (
    EllipticCurve,
    PolynomialRing,
    discrete_log,
    matrix
)

# action of square root of -1
def i_action(P, zeta4):
    F = P.base_ring()
    E = P.curve()
    assert zeta4 in F
    assert E == EllipticCurve(F, [1,0]) # P should be on the curve y^2 = x^3 + x
    X, Y, Z = P
    return E([-X, zeta4*Y, Z])

# Frobenius endomorphism
def Frobenius(P):
    p = P.base_ring().characteristic()
    E = P.curve()
    X, Y, Z = P
    return E([X**p, Y**p, Z**p])

# return retP s.t. 2*retP = P. Note that retP is over an extension field.
def half_point(P, F2):
    F = P.base_ring()
    E = P.curve()
    assert E == EllipticCurve(F, [1,0]) # P should be on the curve y^2 = x^3 + x
    assert F.is_subring(F2)

    E2 = EllipticCurve(F2, [1, 0])
    R = PolynomialRing(F2, name="X")
    X = R.gens()[0]
    if P.is_zero():
        f = X**3 + X
    else:
        f = P[2]*(X**2 - 1)**2 - P[0]*4*(X**3 + X)
    xs = f.roots(multiplicities=False)
    assert len(xs) > 0
    x = xs[0]
    y = (x**3 + x).sqrt()
    retP = E2([x, y])
    if not E(2*retP) == P:
        retP = -retP
    assert E(2*retP) == P
    return retP

# the action of (i + j)/2
def i_j_2_action(P, zeta4, F2):
    F = P.base_ring()
    E = P.curve()
    assert E == EllipticCurve(F, [1,0]) # P should be on the curve y^2 = x^3 + x
    halfP = half_point(P, F2)
    return E(i_action(halfP, zeta4) + Frobenius(halfP))

# the action of (1 + ij)/2
def one_ij_2_action(P, zeta4, F2):
    F = P.base_ring()
    E = P.curve()
    assert E == EllipticCurve(F, [1,0]) # P should be on the curve y^2 = x^3 + x
    halfP = half_point(P, F2)
    return E(halfP + i_action(Frobenius(halfP), zeta4))

# the action of a + bi + c(i + j)/2 + d(1 + ij)/2
def action(alpha, P, zeta4, F2):
    F = P.base_ring()
    E = P.curve()
    assert E == EllipticCurve(F, [1,0]) # P should be on the curve y^2 = x^3 + x
    a, b, c, d = alpha
    ret = a*P
    ret += b*i_action(P, zeta4)
    ret += c*i_j_2_action(P, zeta4, F2)
    ret += d*one_ij_2_action(P, zeta4, F2)
    return ret

# matrix of the multiplication by alpha w.r.t. basis of N-torsion subgroup
def action_matrix(alpha, basis, N, zeta4, F2):
    P, Q = basis
    aP, aQ = [action(alpha, R, zeta4, F2) for R in [P, Q]]
    a, b, c, d = bi_dlp(aP, aQ, P, Q, N)
    assert a*P + b*Q == aP and c*P + d*Q == aQ
    return matrix([[a, c], [b, d]])

# matrices of the multiplications by i, (i + j)/2, (1 + ij)/2
def action_matrices(basis, N, zeta4, F2):
    one = [1,0,0,0]
    Ms = []
    for alpha in [one[-i:]+one[:-i] for i in range(1, 4)]:
        Ms.append(action_matrix(alpha, basis, N, zeta4, F2))
    return Ms

# return a, b, c, d s.t P = [a]R + [b]S, Q = [c]R + [d]S
def bi_dlp(P, Q, R, S, n):
    e0 = R.weil_pairing(S, n)
    e1 = P.weil_pairing(S, n)
    e2 = R.weil_pairing(P, n)
    e3 = Q.weil_pairing(S, n)
    e4 = R.weil_pairing(Q, n)

    a = discrete_log(e1, e0, n)
    b = discrete_log(e2, e0, n)
    c = discrete_log(e3, e0, n)
    d = discrete_log(e4, e0, n)

    return a, b, c, d

