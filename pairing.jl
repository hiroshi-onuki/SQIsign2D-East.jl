using Nemo
using Pkg
using BenchmarkTools
Pkg.activate(".")
import SQIsign2D: Montgomery, PointOrder, Ladder, DefintionField, IsRational, Proj1, IsInfinity,
        OddIsogeny

struct Point{T <: RingElem}
    X::T
    Y::T
    Z::T
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

function Dbl(Mont::Montgomery{T}, P::Point{T}) where T <: RingElem
    A = Mont.A
    lam1 = 3*P.X^2 + 2*A*P.X*P.Z + P.Z^2
    lam2 = 2*P.Y*P.Z
    X = lam1^2*lam2*P.Z - lam2^3*(A*P.Z + 2*P.X)
    Y = lam1 * (lam2^2 * (3*P.X + A*P.Z) - lam1^2*P.Z) - P.Y * lam2^3
    Z = P.Z * lam2^3
    return Point(X, Y, Z)
end

function xDbl(Mont::Montgomery{T}, P::Point{T}, lam1::T, lam2::T) where T <: RingElem
    A = Mont.A
    X = lam1^2*lam2*P.Z - lam2^3*(A*P.Z + 2*P.X)
    Z = P.Z * lam2^3
    return X, Z
end

function hPPQ(Mont::Montgomery{T}, P::Point{T}, Q::Point{T}, lam1::T, lam2::T) where T <: RingElem
    P.Y == 0 && return Q.X*P.Z - P.X*Q.Z, P.Z*Q.Z
    A = Mont.A
    return lam2^2 * (Q.Y*P.Z - P.Y*Q.Z) - lam1*lam2 * (Q.X*P.Z - P.X*Q.Z), lam2^2 * (Q.X*P.Z + 2*P.X*Q.Z + A*P.Z*Q.Z) - lam1^2*P.Z*Q.Z
end

function MillerFuncPowerTwo(Mont::Montgomery{T}, P::Point{T}, Q::Point{T}, e::Integer) where T <: RingElem
    F = DefintionField(Mont)
    A = Mont.A
    X, Y, Z = P.X, P.Y, P.Z
    f1, f2 = F(1), F(1)
    for i in 1:e-1
        lam1 = 3*X^2 + 2*A*X*Z + Z^2
        lam2 = 2*Y*Z
        lam22 = lam2^2
        h1 = lam22 * (Q.Y*Z - Y*Q.Z) - lam1*lam2 * (Q.X*Z - X*Q.Z)
        h2 = lam22 * (Q.X*Z + 2*X*Q.Z + A*Z*Q.Z) - lam1^2*Z*Q.Z
        f1 = f1^2 * h1
        f2 = f2^2 * h2
        if i < e-1
            lam12Z = lam1^2*Z
            lam23 = lam2^3
            X, Y, Z = lam12Z*lam2 - lam23*(A*Z + 2*X),
                lam1 * (lam22 * (3*X + A*Z) - lam12Z) - Y * lam23,
                Z * lam23
        else
            X, Z = lam1^2*Z - lam22*(A*Z + 2*X), Z * lam22
        end
    end
    f1 = f1^2 * (Q.X*Z - X*Q.Z)
    f2 = f2^2 * Z*Q.Z
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
P = Point(xP, sqrt(xP^3 + A*xP^2 + xP), Fp2(1))
Q = Point(xQ, sqrt(xQ^3 + A*xQ^2 + xQ), Fp2(1))
T = P
for _ in 1:e-1
    global T = Dbl(E, T)
end
println("T = ", T)
w = WeilPairingPowerTwo(E, P, Q, e)
println(w^(ZZ(2)^(e-1)) == Fp2(-1))
r = @benchmark WeilPairingPowerTwo(E, P, Q, e)
println(r)
