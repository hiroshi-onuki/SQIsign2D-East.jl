# Montgomery cofficient
struct Montgomery{T <: RingElem}
    A::T
end

coefficient(Mont::Montgomery) = Mont.A
Base.:(==)(Mont1::Montgomery, Mont2::Montgomery) = Mont1.A == Mont2.A

function IsRational(Mont::Montgomery{T}, x::Proj1{T}) where T <: RingElem
    return is_square(x.X^3*x.Z + Mont.A*x.X^2*x.Z^2 + x.X*x.Z^3)
end

# return x(P + Q) from x(P), x(Q), x(P - Q)
function DiffAdd(Mont::Montgomery{T}, P::Proj1{T}, Q::Proj1{T}, PmQ::Proj1{T}) where T <: RingElem
    IsInfinity(P) && return copy(Q)
    IsInfinity(Q) && return copy(P)
    X = PmQ.Z * (P.X*Q.X - P.Z*Q.Z)^2
    Z = PmQ.X * (P.X*Q.Z - P.Z*Q.X)^2
    return Proj1(X, Z)
end

# return x(2P) from A, x(P)
function Dbl(Mont::Montgomery{T}, P::Proj1{T}) where T <: RingElem
    A = Mont.A
    X2 = P.X^2
    Z2 = P.Z^2
    XZ = P.X*P.Z
    X = (X2 - Z2)^2
    Z = 4*XZ*(X2 + A*XZ + Z2)
    return Proj1(X, Z)
end

# return x(nP) from E, x(P), n
function Ladder(E::Montgomery{T}, P::Proj1{T}, n::Integer) where T <: RingElem
    IsInfinity(P) && return InfPoint(DefintionField(E))
    n == 0 && return InfPoint(DefintionField(E))
    n < 0 && (n = abs(n))
    n == 1 && return copy(P)
    P0 = copy(P)
    P1 = Dbl(E, P)
    n == 2 && return P1

    t = n >> 1
    b = BigInt(1)
    while t != 1
        t >>= 1
        b <<= 1 
    end

    while b != 0
        if n & b == 0
            P1 = DiffAdd(E, P0, P1, P)
            P0 = Dbl(E, P0)
        else
            P0 = DiffAdd(E, P0, P1, P)
            P1 = Dbl(E, P1)
        end
        b >>= 1
    end
    return P0
end

# return a point of order l in an elliptic curve of order ord
function PointOrder(E::Montgomery{T}, l::Integer, ord::Integer) where T <: RingElem
    n = ord
    while n % l == 0
        n = div(n, l)
    end
    while true
        x = rand(DefintionField(E))
        if IsRational(E, Proj1(x))
            P = Ladder(E, Proj1(x), n)
            if !IsInfinity(P)
                lP = Ladder(E, P, l)
                while !IsInfinity(lP)
                    P = lP
                    lP = Ladder(E, P, l)
                end
                return P
            end
        end
    end
end

# The Montgomery coefficient of the codomain of an isogeny of odd degree
function OddIsogeny(Mont::Montgomery{T}, ker::Proj1{T}, d::Integer, Ps::Vector{Proj1{T}} = Proj1{T}[]) where T <: RingElem
    F = DefintionField(Mont)
    A = Mont.A
    @assert (d % 2 == 1)
    s = div(d, 2)
    K = ker
    s >= 2 && (R = Dbl(Mont, ker))
    dp = F(1)
    dm = F(1)
    imXs = [one(F) for _ in 1:length(Ps)]
    imZs = [one(F) for _ in 1:length(Ps)]
    for i in 1:s
        tp = K.X + K.Z
        tm = K.X - K.Z
        dp *= tp
        dm *= tm
        for j in 1:length(Ps)
            mp = (Ps[j].X - Ps[j].Z) * tp
            pm = (Ps[j].X + Ps[j].Z) * tm
            imXs[j] *= mp + pm
            imZs[j] *= mp - pm
        end
        if i < s
            K, R = R, DiffAdd(Mont, R, ker, K)
        end
    end
    dp = dp^8
    dm = dm^8
    dp *= (A + 2)^d
    dm *= (A - 2)^d
    Ad = 2*(dp + dm)//(dp - dm)
    length(Ps) == 0 && return Montgomery(Ad)
    imPs = [Proj1(Ps[i].X*imXs[i]^2, Ps[i].Z*imZs[i]^2) for i in 1:length(Ps)]
    return Montgomery(Ad), imPs
end

function DefintionField(Mont::Montgomery{T}) where T <: RingElem
    return parent(Mont.A)
end

function j_invarint(Mont::Montgomery{T}) where T <: RingElem
    A = Mont.A
    return 256*(A^2 - 3)^3/(A^2 - 4)
end
