# double of P. Alorithm 4 in DMPR2023
function double(tnull::ThetaNullLv2{T}, P::ThetaPtLv2{T}) where T <: RingElem
    lam1, lam2, lam3, lamd1, lamd2, lamd3 = precomputation!(tnull)
    x, y, z, w = Hadamard(square(P))
    x2 = x^2
    y2 = lamd1*y^2
    z2 = lamd2*z^2
    w2 = lamd3*w^2
    xd, yd, zd, wd = Hadamard(x2, y2, z2, w2)
    return ThetaPtLv2(xd, lam1*yd, lam2*zd, lam3*wd)
end

# differential addition of P and Q. Algorithm 3 in DMPR2023
function diff_add(tnull::ThetaNullLv2{T}, P::ThetaPtLv2{T}, Q::ThetaPtLv2{T}, PmQ::ThetaPtLv2{T}) where T <: RingElem
    _, _, _, lamd1, lamd2, lamd3 = precomputation!(tnull)
    xP, yP, zP, wP = Hadamard(square(P))
    xQ, yQ, zQ, wQ = Hadamard(square(Q))
    xPQ = xP*xQ
    yPQ = lamd1*yP*yQ
    zPQ = lamd2*zP*zQ
    wPQ = lamd3*wP*wQ
    xPQ, yPQ, zPQ, wPQ = Hadamard(xPQ, yPQ, zPQ, wPQ)
    xyPmQ = PmQ[1]*PmQ[2]
    zwPmQ = PmQ[3]*PmQ[4]
    x = xPQ*zwPmQ*PmQ[2]
    y = yPQ*zwPmQ*PmQ[1]
    z = zPQ*xyPmQ*PmQ[4]
    w = wPQ*xyPmQ*PmQ[3]
    return ThetaPtLv2(x, y, z, w)
end

# return [2^e]P by double
function double_iter(tnull::ThetaNullLv2{T}, P::ThetaPtLv2{T}, e::Integer) where T <: RingElem
    for _ in 1:e
        P = double(tnull, P)
    end
    return P
end

# return [m]P by Montgomey ladder
function ladder(tnull::ThetaNullLv2{T}, m::Integer, P::ThetaPtLv2{T}) where T <: RingElem
    m == 0 && return ThetaPtLv2([tnull[i] for i in 1:4])
    m == 1 && return P
    m == 2 && return double(tnull, P)

    t = m >> 1
    b = BigInt(1)
    while t != 1
        t >>= 1
        b <<= 1 
    end

    P0, P1 = P, double(tnull, P)
    while b != 0
        if m & b == 0
            P0, P1 = double(tnull, P0), diff_add(tnull, P0, P1, P)
        else
            P1, P0 = double(tnull, P1), diff_add(tnull, P0, P1, P)
        end
        b >>= 1
    end
    return P0
end