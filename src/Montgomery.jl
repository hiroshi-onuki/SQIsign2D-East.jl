
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
