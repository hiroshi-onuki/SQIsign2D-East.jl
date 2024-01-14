# random point on a Montgomery curve: y^2 = x^3 + Ax^2 + x
function random_point(A::T) where T <: RingElem
    F = parent(A)
    while true
        x = rand(F)
        is_square(x^3 + A*x^2 + x) && return Proj1(x)
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

# return [m]P
function Ladder(m::Integer, P::Proj1{T}, a24::Proj1{T}) where T <: RingElem
    if P.X == 0
        m % 2 == 0 ? (return InfPoint(T)) : return P
    end
    m == 0 && return InfPoint(T)
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

# return P + [m]Q
function Ladder3pt(m::Integer, P::Proj1{T}, Q::Proj1{T}, QmP::Proj1{T}, a24::Proj1{T}) where T <: RingElem
    m < 0 && error("m in Ladder3pt must be nonnegative")
    IsInfinity(QmP) && error("Q == P is not allowed")

    P0 = Q;
    P1 = P;
    P2 = QmP;
    t = m
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

# Montgomery coefficient of a24
function Montgomery_coeff(a24::Proj1)
    a = a24.X + a24.X
    a -= a24.Z
    a += a
    return a//a24.Z
end

# the j-invariant of a24
function jInvariant(a24::Proj1)
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

# 2^e-isogey with kernel <K>. A = (a + 2)/4.
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
            println("hoge")
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
