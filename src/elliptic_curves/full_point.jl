export Point, double, add, mult, is_infinity

# full point on a Montgomery curve
struct Point{T <: RingElem}
    X::T
    Y::T
    Z::T
end

function Base.:(==)(P::Point{T}, Q::Point{T}) where T <: RingElem
    return P.X*Q.Z == P.Z*Q.X && P.Y*Q.Z == P.Z*Q.Y
end

function Base.:-(P::Point{T}) where T <: RingElem
    return Point(P.X, -P.Y, P.Z)
end

function Point(X::T, Y::T) where T <: RingElem
    F = parent(X)
    return Point(X, Y, F(1))
end

function Point(A::T, XZ::Proj1{T}) where T <: RingElem
    X, Z = XZ.X, XZ.Z
    XZ = X*Z
    Z2 = Z^2
    Y, Yd = square_root(XZ*(X^2 + A*XZ + Z2))
    X = XZ * Yd
    Z = Z2 * Yd
    return Point(X, Y, Z)
end

function is_infinity(P::Point{T}) where T <: RingElem
    return P.Z == 0
end

function infinity_full_point(F::T) where T <: Ring
    return Point(F(0), F(1), F(0))
end

function double(P::Point{T}, A::Proj1{T}) where T <: RingElem
    X, Y, Z = P.X, P.Y, P.Z
    XY = X * Y
    YZ = Y * Z
    M = X^2
    M = M + M + M
    M += Z^2
    M *= A.Z
    AXZ = A.X * X * Z
    M += AXZ + AXZ
    W = A.Z * YZ
    W += W
    W2 = W^2
    CXY = A.Z * XY
    X2 = W * (A.X * YZ + CXY + CXY)
    X2 += X2
    X2 = M^2 - X2
    CWXY = W * CXY
    CW2Y2 = A.Z * W2 * Y^2
    Y2 = M * (CWXY + CWXY - X2) - CW2Y2 - CW2Y2
    Z2 = W2 * W
    return Point(X2 * W, Y2, Z2)
end

function add(P::Point{T}, Q::Point{T}, A::Proj1{T}) where T <: RingElem
    P.Z == 0 && return Q
    Q.Z == 0 && return P
    X1, Y1, Z1 = P.X, P.Y, P.Z
    X2, Y2, Z2 = Q.X, Q.Y, Q.Z
    X1Z2 = X1 * Z2
    X2Z1 = X2 * Z1
    Y1Z2 = Y1 * Z2
    Y2Z1 = Y2 * Z1
    if X1Z2 == X2Z1
        if Y1Z2 == Y2Z1
            return double(P, A)
        else
            F = parent(A.X)
            return infinity_full_point(F)
        end
    end
    U = X2Z1 - X1Z2
    V = Y2Z1 - Y1Z2
    ZZ = Z1 * Z2
    W = ZZ * A.Z
    U2 = U^2
    U3 = U2 * U
    X3 = V^2 * W - U2 * (A.X * ZZ + X1Z2 * A.Z + X2Z1 * A.Z)
    Y3 = V * (X1Z2 * U2 * A.Z - X3) - Y1Z2 * U3 * A.Z
    Z3 = U3 * W
    return Point(X3 * U, Y3, Z3)
end

function mult(m::Integer, P::Point{T}, A::Proj1{T}) where T <: RingElem
    m < 0 && return -mult(-m, P, A)
    F = parent(A.X)
    m == 0 && return infinity_full_point(F)
    m == 1 && return P
    m == 2 && return double(P, A)

    t = m >> 1
    b = BigInt(1)
    while t != 1
        t >>= 1
        b <<= 1 
    end

    R = P
    while b != 0
        R = double(R, A)
        if m & b != 0
            R = add(R, P, A)
        end
        b >>= 1
    end
    return R
end