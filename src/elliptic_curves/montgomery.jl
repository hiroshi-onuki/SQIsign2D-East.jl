export xDBL, xADD, xDBLADD, xDBLe, ladder, ladder3pt, x_add_sub,
    linear_comb_2_e, random_point, random_point_order_2power,
    Montgomery_coeff, A_to_a24, a24_to_A, jInvariant_a24, jInvariant_A,
    two_e_iso, odd_isogeny, torsion_basis, isomorphism_Montgomery,
    Montgomery_normalize, complete_baisis

# random point on a Montgomery curve: y^2 = x^3 + Ax^2 + x
function random_point(A::T) where T <: RingElem
    F = parent(A)
    while true
        x = rand(F)
        is_square(x^3 + A*x^2 + x) && return Proj1(x)
    end
end

# random point on a Montgomery curve with order 2^e
function random_point_order_2power(A::T, curve_order::Integer, e::Integer) where T <: RingElem
    F = parent(A)
    a24 = Proj1(A + 2, F(4))
    n = div(curve_order, BigInt(2)^e)
    while true
        P = random_point(A)
        P = ladder(n, P, a24)
        if !is_infinity(xDBLe(P, a24, e-1))
            return P
        end
    end
end

# return [2]P on Mont_A with a24 = (A + 2)/4.
function xDBL(P::Proj1{T}, a24::Proj1{T}) where T <: RingElem
    t0 = P.X - P.Z
    t1 = P.X + P.Z
    t0 = t0^2
    t1 = t1^2
    Z = a24.Z*t0
    X = Z*t1
    t1 -= t0
    t0 = a24.X*t1
    Z += t0
    Z *= t1

    return Proj1(X, Z);
end

# return [2^e]P on Mont_A with a24 = (A + 2)/4.
function xDBLe(P::Proj1{T}, a24::Proj1{T}, e::Int) where T <: RingElem
    outP = P
    for _ in 1:e
        outP = xDBL(outP, a24)
    end

    return outP
end

# return P + Q from P, Q, Q-P
function xADD(P::Proj1{T}, Q::Proj1{T}, QmP::Proj1{T}) where T <: RingElem
    a = P.X + P.Z
    b = P.X - P.Z
    c = Q.X + Q.Z
    d = Q.X - Q.Z
    a *= d
    b *= c
    c = a + b
    d = a - b
    c = c^2
    d = d^2
    return Proj1(QmP.Z*c, QmP.X*d)
end

# return 2P and P + Q
function xDBLADD(P::Proj1{T}, Q::Proj1{T}, QmP::Proj1{T}, a24::Proj1{T}) where T <: RingElem
    t0 = P.X + P.Z;
    t1 = P.X - P.Z;
    tPX = t0^2;
    t2 = Q.X - Q.Z;
    PpQX = Q.X + Q.Z;
    t0 *= t2;
    tPZ = t1^2;
    t1 *= PpQX;
    t2 = tPX - tPZ;
    tPZ *= a24.Z
    tPX *= tPZ;
    PpQX = a24.X * t2;
    PpQZ = t0 - t1;
    tPZ = PpQX + tPZ;
    PpQX = t0 + t1;
    tPZ = tPZ * t2;
    PpQZ = PpQZ^2;
    PpQX = PpQX^2;
    PpQZ *= QmP.X;
    PpQX *= QmP.Z;

    return Proj1(tPX, tPZ), Proj1(PpQX, PpQZ);
end

# return x(P + Q) or x(P - Q)
function x_add_sub(P::Proj1{T}, Q::Proj1{T}, a24::Proj1{T}) where T <: RingElem
    is_infinity(P) && return Q
    is_infinity(Q) && return P

    A = a24.X + a24.X
    A -= a24.Z
    A += A
    C = a24.Z
    X1X2 = P.X * Q.X
    Z1Z2 = P.Z * Q.Z
    X1Z2 = P.X * Q.Z
    Z1X2 = P.Z * Q.X
    a = (X1Z2 - Z1X2)^2 * C
    b = (X1X2 + Z1Z2) * (X1Z2 + Z1X2) * C + (A + A) * X1X2 * Z1Z2
    c = (X1X2 - Z1Z2)^2 * C
    d = square_root(b^2 - a*c)
    return Proj1(b + d, a)
end

# return [m]P
function ladder(m::Integer, P::Proj1{T}, a24::Proj1{T}) where T <: RingElem
    m == 0 && return infinity_point(parent(P.X))
    m == 1 && return P
    m == 2 && return xDBL(P, a24)

    t = m >> 1
    b = BigInt(1)
    while t != 1
        t >>= 1
        b <<= 1 
    end

    P0, P1 = P, xDBL(P, a24)
    while b != 0
        if m & b == 0
            P0, P1 = xDBLADD(P0, P1, P, a24)
        else
            P1, P0 = xDBLADD(P1, P0, P, a24)
        end
        b >>= 1
    end
    return P0
end

# return x(P + [m]Q) from x(P), x(Q), x(Q-P)
function ladder3pt(m::Integer, P::Proj1{T}, Q::Proj1{T}, QmP::Proj1{T}, a24::Proj1{T}) where T <: RingElem
    m < 0 && error("m in Ladder3pt must be nonnegative")
    m == 0 && return P
    is_infinity(QmP) && error("Q == P is not allowed")

    P0 = Q;
    P1 = P;
    P2 = QmP;
    t = BigInt(m)
    while (t != 0)
        if (t & 1 == 1)
            P0, P1 = xDBLADD(P0, P1, P2, a24);
        else
            P0, P2 = xDBLADD(P0, P2, P1, a24);
        end
        t >>= 1;
    end
    return P1;
end

# return x([a]P + [b]Q) from x(P), x(Q), x(Q-P), where ord(P) = ord(Q) = 2^e
function linear_comb_2_e(a::Integer, b::Integer, xP::Proj1{T}, xQ::Proj1{T}, xQmP::Proj1{T}, a24::Proj1{T}, e::Int) where T <: RingElem
    a = a % (BigInt(2)^e)
    b = b % (BigInt(2)^e)
    a < 0 && (a += BigInt(2)^e)
    b < 0 && (b += BigInt(2)^e)
    a == 0 && return ladder(b, xQ, a24)
    b == 0 && return ladder(a, xP, a24)
    g = gcd(a, b)
    f = 0
    while g & 1 == 0
        g >>= 1
        f += 1
    end
    a >>= f
    b >>= f
    if a & 1 == 0
        c = invmod(b, BigInt(2)^e)
        a = (a * c) % (BigInt(2)^e)
        xR = ladder3pt(a, xQ, xP, xQmP, a24)
        xR = ladder(b, xR, a24)
    else
        c = invmod(a, BigInt(2)^e)
        b = (b * c) % (BigInt(2)^e)
        xR = ladder3pt(b, xP, xQ, xQmP, a24)
        xR = ladder(a, xR, a24)
    end
    return xDBLe(xR, a24, f)
end

# Montgomery coefficient of a24
function Montgomery_coeff(a24::Proj1)
    a = a24.X + a24.X
    a -= a24.Z
    a += a
    return a//a24.Z
end

function a24_to_A(a24::Proj1)
    a = a24.X + a24.X
    a -= a24.Z
    a += a
    return Proj1(a, a24.Z)
end

function A_to_a24(A::T) where T <: RingElem
    F = parent(A)
    return Proj1(A + 2, F(4))
end

function A_to_a24(A::Proj1)
    Z2 = A.Z + A.Z
    Z4 = Z2 + Z2
    return Proj1(A.X + Z2, Z4)
end

# the j-invariant of a24
function jInvariant_a24(a24::Proj1)
    j = a24.X + a24.X - a24.Z
    j += j
    j = j^2
    t1 = a24.Z^2
    t0 = t1 + t1
    t0 = j - t0
    t0 -= t1
    j = t0 - t1
    t1 = t1^2
    j *= t1
    t0 += t0
    t0 += t0
    t1 = t0^2
    t0 *= t1
    t0 += t0
    t0 += t0
    j = 1//j
    j *= t0

    return j
end

function jInvariant_A(A::Proj1)
    return jInvariant_a24(A_to_a24(A))
end

function jInvariant_A(A::T) where T <: RingElem
    return jInvariant_a24(A_to_a24(A))
end

function two_iso(a24::Proj1{T}, P::Proj1{T}) where T <: RingElem
    if P.X == 0
        a24, _ = two_iso_zero(a24, Proj1{T}[])
    else
        a24 = two_iso_curve(P)
    end
    return a24
end

# 2-isogey. return (A + 2)/4, where Mont_A = E/<P> with ord(P) = 2.
function two_iso_curve(P::Proj1)
    X = P.X^2;
    Z = P.Z^2;
    X = Z - X;

    return Proj1(X, Z);
end

# The image of 2-isogey. return f(Q), where f is an isogeny with kernel <P>.
function two_iso_eval(P::Proj1{T}, Q::Proj1{T}) where T <: RingElem
    t0 = P.X + P.Z;
    t1 = P.X - P.Z;
    t2 = Q.X + Q.Z;
    t3 = Q.X - Q.Z;
    t0 *= t3;
    t1 *= t2;
    t2 = t0 + t1;
    t3 = t0 - t1;
    X = Q.X * t2;
    Z = Q.Z * t3;

    return Proj1(X, Z);
end

# 2-isogeny with kernel (0, 0)
function two_iso_zero(a24::Proj1{T}, Qs::Vector{Proj1{T}}) where T <: RingElem
    a = a24.X + a24.X - a24.Z
    d = sqrt(a24.X*(a24.X - a24.Z))
    d += d
    d2 = d+d
    a24d = Proj1(-a + d, d2)
    a += a

    retQ = Proj1{T}[]
    for Q in Qs
        XZ = Q.X*Q.Z
        X = Q.X^2 + Q.Z^2
        X *= a24.Z
        X += a*XZ
        Z = d2*XZ
        push!(retQ, Proj1(X, Z))
    end

    return a24d, retQ
end

# 4-isogey. return (A + 2)/4 and (K1, K2, K3), where Mont_A = E/<P> with ord(P) = 4.
function four_iso_curve(P::Proj1)
    K2 = P.X - P.Z;
    K3 = P.X + P.Z;
    K1 = P.Z^2;
    K1 += K1;
    Z = K1^2;
    K1 += K1;
    X = P.X^2;
    X += X;
    X = X^2;

    return Proj1(X, Z), K1, K2, K3;
end

# The image of 4-isogey using an output (K1, K2, K3) of four_iso_curve().
function four_iso_eval(K1::T, K2::T, K3::T, Q::Proj1{T}) where T <: RingElem
    t0 = Q.X + Q.Z;
    t1 = Q.X - Q.Z;
    QX = t0 * K2;
    QZ = t1 * K3;
    t0 *= t1;
    t0 *= K1;
    t1 = QX + QZ;
    QZ = QX - QZ;
    t1 = t1^2;
    QZ = QZ^2;
    QX = t0 + t1;
    t0 = QZ - t0;
    X = QX * t1;
    Z = QZ * t0;

    return Proj1(X, Z);
end

# 4-isogeny with kernel <(1, -)>
four_iso_curve_one(a24::Proj1) = Proj1(a24.X, a24.X - a24.Z)

# 4-isogeny with kernel <(-1, -)>
four_iso_curve_neg_one(a24::Proj1) = Proj1(a24.Z, a24.X)

# evaluation under 4-isogeny with kernel <(1, -)>
function four_iso_eval_one(a24::Proj1{T}, P::Proj1{T}) where T <: RingElem
    t0 = (P.X - P.Z)^2
    t1 = P.X * P.Z
    t1 += t1
    t1 += t1
    X = (t0 + t1) * (a24.Z*t0 + a24.X*t1)
    Z = a24.Z - a24.X
    Z *= t0
    Z *= t1
    return Proj1(X, Z)
end

# evaluation under 4-isogeny with kernel <(-1, -)>
function four_iso_eval_neg_one(a24::Proj1{T}, P::Proj1{T}) where T <: RingElem
    Q = four_iso_eval_one(Proj1(a24.Z-a24.X, a24.Z), Proj1(-P.X, P.Z))
    return Proj1(-Q.X, Q.Z)
end

# 4-isogeny with kernel <ker>
function four_iso(a24::Proj1{T}, ker::Proj1{T}, Qs::Vector{Proj1{T}}) where T <: RingElem
    imQs = Vector{Proj1{T}}(undef, length(Qs))
    if ker.X == ker.Z
        for i in 1:length(Qs)
            imQs[i] = four_iso_eval_one(a24, Qs[i])
        end
        a24 = four_iso_curve_one(a24)
    elseif ker.X == -ker.Z
        for i in 1:length(Qs)
            imQs[i] = four_iso_eval_neg_one(a24, Qs[i])
        end
        a24 = four_iso_curve_neg_one(a24)
    else
        a24, K1, K2, K3 = four_iso_curve(ker)
        for i in 1:length(Qs)
            imQs[i] = four_iso_eval(K1, K2, K3, Qs[i])
        end
    end
    return a24, imQs
end

# 2^e-isogeny with kernel <K>. A = (a + 2)/4.
function two_e_iso(a24::Proj1{T}, P::Proj1{T}, e::Int, Qs::Vector{Proj1{T}}) where T <: RingElem
    k = e-2
    while k >= 0
        ker = xDBLe(P, a24, k)
        if ker.X == ker.Z
            k != 0 && (P = four_iso_eval_one(a24, P))
            for i in 1:length(Qs)
                Qs[i] = four_iso_eval_one(a24, Qs[i])
            end
            a24 = four_iso_curve_one(a24)
        elseif ker.X == -ker.Z
            k != 0 && (P = four_iso_eval_neg_one(a24, P))
            for i in 1:length(Qs)
                Qs[i] = four_iso_eval_neg_one(a24, Qs[i])
            end
            a24 = four_iso_curve_neg_one(a24)
        else
            a24, K1, K2, K3 = four_iso_curve(ker)
            k != 0 && (P = four_iso_eval(K1, K2, K3, P))
            for i in 1:length(Qs)
                Qs[i] = four_iso_eval(K1, K2, K3, Qs[i])
            end
        end
        k -= 2
    end
    if k == -1
        if P.X == 0
            a24, Qs = two_iso_zero(a24, Qs)
        else
            a24 = two_iso_curve(P)
            for i in 1:length(Qs)
                Qs[i] = two_iso_eval(P, Qs[i])
            end
        end
    end

    return a24, Qs
end

# 2^e-isogeny using strategy. A = (a + 2)/4.
function two_e_iso(a24::Proj1{T}, P::Proj1{T}, e::Int, Qs::Vector{Proj1{T}}, strategy::Vector{Int}, out_put_neighbor=0) where T <: RingElem
    if out_put_neighbor == 1
        K = xDBLe(P, a24, e-1)
        a24_neighbor = two_iso(a24, K)
    end
    S = [div(e, 2)]
    Ps = vcat(Qs, [P])
    i = 1
    while length(S) > 0
        h = pop!(S)
        K = pop!(Ps)
        if h == 1
            if length(S) == 0 && out_put_neighbor == -1
                a24_neighbor = two_iso(a24, xDBL(K, a24))
            end
            a24, Ps = four_iso(a24, K, Ps)
            S = [h - 1 for h in S]
        else
            push!(S, h)
            push!(Ps, K)
            K = xDBLe(K, a24, 2*strategy[i])
            push!(S, h - strategy[i])
            push!(Ps, K)
            i += 1
        end
    end
    if out_put_neighbor == 0
        return a24, Ps
    else
        return a24, Ps, a24_neighbor
    end
end

# isogeny of odd degree d
function odd_isogeny(a24::Proj1{T}, ker::Proj1{T}, d::Integer, Qs::Vector{Proj1{T}}) where T <: RingElem
    F = parent(a24.X)
    s = div(d, 2)
    K = ker
    s >= 2 && (R = xDBL(ker, a24))
    A = 1
    C = 1
    imQs = [[F(1), F(1)] for _ in Qs]
    for i in 1:s 
        tp = K.X + K.Z
        tm = K.X - K.Z
        A *= tp
        C *= tm
        for i in 1:length(Qs)
            mp = (Qs[i].X - Qs[i].Z) * tp
            pm = (Qs[i].X + Qs[i].Z) * tm
            imQs[i][1] *= mp + pm
            imQs[i][2] *= mp - pm
        end
        if i < s
            K, R = R, xADD(R, ker, K)
        end
    end
    A = A^8
    C = C^8
    A *= a24.X^d
    C *= (a24.X - a24.Z)^d
    retQs = Proj1{T}[Proj1(Qs[i].X*imQs[i][1]^2, Qs[i].Z*imQs[i][2]^2) for i in 1:length(Qs)]
    return Proj1(A, A - C), retQs
end

# Algorithm 2 in SQIsign documentation
# return a fixed basis (P, Q) of E[2^e] from P
function complete_baisis(a24::Proj1{T}, P::Proj1{T}, Pd::Proj1{T}, x::T, e::Int) where T <: RingElem
    F = parent(a24.X)
    p = Integer(characteristic(F))
    N = (p + 1) >> e
    A = a24_to_A(a24)
    i = gen(F)
    Q = Proj1(x)
    while true
        x += i
        if is_square(A.Z * x * (A.Z * (x^2 + 1) + A.X * x))
            Q = Proj1(x)
            Q = ladder(N, Q, a24)
            Qd = xDBLe(Q, a24, e-1)
            if !is_infinity(Qd) && Qd != Pd
                break
            end
        end
    end
    PQ = x_add_sub(P, Q, a24)
    return P, Q, PQ
end

# Algorithm 3 in SQIsign documentation
# return a fixed basis of E[2^e]
function torsion_basis(a24::Proj1{T}, e::Int) where T <: RingElem
    F = parent(a24.X)
    p = Integer(characteristic(F))
    N = (p + 1) >> e
    A = a24_to_A(a24)
    i = gen(F)
    x = F(1)
    P = Proj1(x)
    Pd = Proj1(x)
    while true
        x += i
        if is_square(A.Z * x * (A.Z * (x^2 + 1) + A.X * x))
            P = Proj1(x)
            P = ladder(N, P, a24)
            Pd = xDBLe(P, a24, e-1)
            if !is_infinity(Pd)
                break
            end
        end
    end
    return complete_baisis(a24, P, Pd, x, e)
end

# isomorphism from Montgomery curnve with a24 to Montgomery curve mapping P4 to (1, *)
function isomorphism_Montgomery(a24::Proj1{T}, P4::Proj1{T}, Ps::Vector{Proj1{T}}) where T <: RingElem
    P2 = xDBL(P4, a24)
    A = a24_to_A(a24)
    u = P4.X * P2.Z - P2.X * P4.Z
    Ad = Proj1((A.X * P2.Z + 3*P2.X*A.Z) * P4.Z, u * A.Z)
    imPs = [Proj1(P4.Z * (P.X * P2.Z - P2.X * P.Z), u * P.Z) for P in Ps]
    return A_to_a24(Ad), imPs
end

# Algorithm 1 in SQIsign documentation
function Montgomery_normalize(a24::Proj1{T}, Ps::Vector{Proj1{T}}) where T <: RingElem
    A = a24_to_A(a24)
    A2 = A.X^2
    C2 = A.Z^2
    twoC2 = C2 + C2
    threeC2 = twoC2 + C2
    d = A2 - twoC2 - twoC2
    d = square_root(d)
    inv_2C2_d = 1/(twoC2 * d)
    Z0 = (A2 + A2) * inv_2C2_d * d
    t1 = (threeC2 + threeC2 + threeC2 - A2) * d
    t2 = (A2 - threeC2) * A.X
    Z1 = (t1 + t2) * inv_2C2_d
    Z2 = (t1 - t2) * inv_2C2_d

    Z = Z0
    lex_order(Z1, Z) && (Z = Z1)
    lex_order(Z2, Z) && (Z = Z2)

    Ad = Proj1(square_root(Z))

    if A == Ad
        R1, R2 = 0, 1
        U1, U2 = 1, 1
    elseif A == -Ad
        R1, R2 = 0, 1
        U1, U2 = -1, 1
    else
        Ad2 = Ad.X^2
        Cd2 = Ad.Z^2
        ACd2 = A2 * Cd2
        AdC2 = Ad2 * C2
        threeCCd2 = Cd2 * threeC2
        nineCCd2 = threeCCd2 + threeCCd2 + threeCCd2
        R1 = (ACd2 + AdC2 - threeCCd2 - threeCCd2) * A.X
        R2 = (ACd2 + AdC2 + AdC2 - nineCCd2) * A.Z
        U1 = Ad.X * A.Z * R2
        R1C = R1 * A.Z
        U2 = (A.X * R2 - R1C - R1C - R1C) * Ad.Z
    end

    images = Vector{Proj1{T}}(undef, length(Ps))
    for i in 1:length(Ps)
        P = Ps[i]
        X = P.X
        Z = P.Z
        images[i] = Proj1(U1*(X*R2 + Z*R1), U2*Z*R2)
    end

    return A_to_a24(Ad), images
end