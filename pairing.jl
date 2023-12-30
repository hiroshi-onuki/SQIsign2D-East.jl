using Nemo
using Pkg
using BenchmarkTools
Pkg.activate(".")
import SQIsign2D: Montgomery, PointOrder, Ladder, DefintionField, IsRational, Proj1, IsInfinity,
        OddIsogeny

struct Point{T <: RingElem}
    X::T
    Y::T
end

function RandomSSCurve(F)
    p = characteristic(F)
    @assert p % 12 == 11
    E = Montgomery(F(0))
    for _ in 1:100
        P = PointOrder(E, 3, BigInt(p + 1)^2)
        E = OddIsogeny(E, P, 3)
    end
    return E
end

function Dbl(Mont::Montgomery{T}, P::Point{T}, lam::T) where T <: RingElem
    A = Mont.A
    X = lam^2 - A - 2*P.X
    Y = lam * (P.X - X) - P.Y
    return Point(X, Y)
end

function hPPQ(Mont::Montgomery{T}, P::Point{T}, Q::Point{T}, lam) where T <: RingElem
    P.Y == 0 && return Q.X - P.X
    A = Mont.A
    return (Q.Y - P.Y - lam * (Q.X - P.X)) / (Q.X + 2*P.X + A - lam^2)
end

function MillerFuncPowerTwoInv(Mont::Montgomery{T}, P::Point{T}, Q::Point{T}, e::Integer) where T <: RingElem
    F = DefintionField(Mont)
    A = Mont.A
    R = P
    f = F(1)
    for i in 1:e-1
        lam = (3*R.X^2 + 2*A*R.X + 1) / (2*R.Y)
        f = f^2 * hPPQ(Mont, R, Q, lam)
        R = Dbl(Mont, R, lam)
    end
    f = f^2 * (Q.X - R.X)
    return f
end

function WeilPairingPowerTwoInv(Mont::Montgomery{T}, P::Point{T}, Q::Point{T}, e::Integer) where T <: RingElem
    fPQ = MillerFuncPowerTwoInv(Mont, P, Q, e)
    fQP = MillerFuncPowerTwoInv(Mont, Q, P, e)
    return fPQ / fQP
end

function MillerFuncPowerTwo(Mont::Montgomery{T}, P::Point{T}, Q::Point{T}, e::Integer) where T <: RingElem
    F = DefintionField(Mont)
    A = Mont.A
    X, Y, Z = P.X, P.Y, F(1)
    f1, f2 = F(1), F(1)
    for i in 1:e-1
        AZ = A*Z
        QXZ = Q.X*Z
        lam1 = (3*X + 2*AZ) * X + Z^2
        lam12Z = lam1^2*Z
        lam2 = 2*Y*Z
        lam22 = lam2^2
        h1 = lam22 * (Q.Y*Z - Y) - lam1*lam2 * (QXZ - X)
        h2 = lam22 * (QXZ + 2*X + AZ) - lam12Z
        f1 = f1^2 * h1
        f2 = f2^2 * h2
        if i < e-1
            lam23 = lam2^3
            X, Y, Z = lam12Z*lam2 - lam23*(AZ + 2*X),
                lam1 * (lam22 * (3*X + AZ) - lam12Z) - Y * lam23,
                Z * lam23
        else
            X, Z = lam1^2*Z - lam22*(AZ + 2*X), Z * lam22
        end
    end
    f1 = f1^2 * (Q.X*Z - X)
    f2 = f2^2 * Z
    return f1, f2
end

function WeilPairingPowerTwo(Mont::Montgomery{T}, P::Point{T}, Q::Point{T}, e::Integer) where T <: RingElem
    fPQ1, fPQ2 = MillerFuncPowerTwo(Mont, P, Q, e)
    fQP1, fQP2 = MillerFuncPowerTwo(Mont, Q, P, e)
    return (fPQ1*fQP2) / (fPQ2*fQP1)
end

function PointOrderPowerTwo(E::Montgomery, e::Integer, ord::Integer)
    n = div(ord, BigInt(2)^e)
    xP = rand(Fp2)
    while true
        xP = rand(Fp2)
        if IsRational(E, Proj1(xP))
            P = Ladder(E, Proj1(xP), n)
            if !IsInfinity(Ladder(E, P, BigInt(2)^(e-1)))
                return P
            end
        end
    end
end

p = 13 * ZZ(2)^126 * ZZ(3)^78 - 1
println("p = ", p)
R, T = PolynomialRing(GF(p), "T")
Fp2, i = FiniteField(T^2 + 1, "i")
e = 126

E = RandomSSCurve(Fp2)
A = E.A
println("E = ", E)

P = PointOrderPowerTwo(E, e, BigInt(p + 1))
Q = PointOrderPowerTwo(E, e, BigInt(p + 1))
while Ladder(E, P, BigInt(2)^(e-1)) == Ladder(E, Q, BigInt(2)^(e-1))
    global Q = PointOrderPowerTwo(E, e, BigInt(p + 1))
end
@assert !IsInfinity(Ladder(E, P, BigInt(2)^(e-1)))
@assert !IsInfinity(Ladder(E, Q, BigInt(2)^(e-1)))
@assert IsInfinity(Ladder(E, P, BigInt(2)^e))
@assert IsInfinity(Ladder(E, Q, BigInt(2)^e))

xP = P.X / P.Z
xQ = Q.X / Q.Z
println("xP = ", xP)
println("xQ = ", xQ)
P = Point(xP, sqrt(xP^3 + A*xP^2 + xP))
Q = Point(xQ, sqrt(xQ^3 + A*xQ^2 + xQ))
w = WeilPairingPowerTwo(E, P, Q, e)
wd = WeilPairingPowerTwoInv(E, P, Q, e)
println(w^(ZZ(2)^(e-1)) == Fp2(-1))
println(w == wd)
r = @benchmark WeilPairingPowerTwo(E, P, Q, e)
println(r)
r = @benchmark WeilPairingPowerTwoInv(E, P, Q, e)
println(r)
