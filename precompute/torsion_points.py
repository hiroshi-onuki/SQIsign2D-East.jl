from sage.all import (
    GF,
    PolynomialRing,
    ZZ,
    is_prime
)

# return GF(p^4), GF(p^2), square root of -1 in GF(p^2)
def calcFields(p):
    assert p % 4 == 3
    R = PolynomialRing(GF(p), name="x")
    x = R.gens()[0]
    Fp2 = GF(p**2, modulus=x**2+1, name="i")
    z = Fp2.random_element()
    while z.is_square():
        z = Fp2.random_element()
    t = ZZ(z + z**p)
    n = ZZ(z**(p+1))
    R = PolynomialRing(ZZ, name="x")
    x = R.gens()[0]
    Fp4 = GF(p**4, modulus=x**4 - t*x**2 + n, name="z")
    Fp2 = Fp4.subfield(2)
    i = Fp2.gen()**((p**2 - 1)//4)
    assert i**2 == -1
    return Fp4, Fp2, i

# return a point whose x-coordinate in Fp2 on E of order l^e, where l is prime.
def point_ord(E, Fp2, is_twist, l, e):
    p = E.base_ring().characteristic()
    if is_twist:
        ord = p - 1
    else:
        ord = p + 1
    assert is_prime(l)
    assert ord % l**e == 0

    while True:
        xP = Fp2.random_element()
        P = E.lift_x(xP)
        y = P.xy()[1]
        if (is_twist and not y in Fp2) or (not is_twist and y in Fp2):
            P = ord//l**e * P
            if not (l**(e-1)*P).is_zero():
                break
    assert (l**e*P).is_zero()

    # this makes this function deterministic when the seed is fixed.
    if P[1][0] >= (p + 1)//2 or (P[1][0] == 0 and P[1][1] >= (p + 1)//2):
        P = -P
    return P

# return a basis of E[l^e], where l is prime.
def basis(E, Fp2, is_twist, l, e):
    P = point_ord(E, Fp2, is_twist, l, e)
    Q = point_ord(E, Fp2, is_twist, l, e)
    while (l**(e-1)*P).weil_pairing(l**(e-1)*Q, l) == 1:
        Q = point_ord(E, Fp2, is_twist, l, e)
    return P, Q

def basis_2e_special(E, Fp2, zeta4, e):
    P = point_ord(E, Fp2, False, 2, e)
    while (2**(e-1)*P).xy()[0] == 0:
        P = point_ord(E, Fp2, False, 2, e)
    Q = E([-P.xy()[0], zeta4*P.xy()[1]])
    assert (2**(e-1)*P).weil_pairing(2**(e-1)*Q, 2) == -1
    return P, Q

# transform an element x in Fp(zeta4) to an element in Fp2(zeta4d)
def Fp2ToFp2d(x, zeta4, zeta4d):
    p = x.base_ring().characteristic()
    return ZZ((x + x**p)/2) + ZZ((x - x**p)/(2*zeta4)) * zeta4d
