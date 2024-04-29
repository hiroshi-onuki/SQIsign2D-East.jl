export Weil_pairing_2power, Tate_pairing_P0, Tate_pairing_iP0, make_pairing_table, Tate_pairings

# Miller function f_{P}(Q)
function Miller_function(A::T, P::Point{T}, Q::Point{T}, e::Integer) where T <: RingElem
    F = parent(A)
    X, Y, Z = P.X, P.Y, P.Z
    f1, f2 = F(1), F(1)
    for i in 1:e-1
        AZ = A*Z
        QXZ = Q.X*Z
        QZX = Q.Z*X
        QZY = Q.Z*Y
        QZZ = Q.Z*Z
        lam1 = (3*X + 2*AZ) * X + Z^2
        lam12Z = lam1^2*Z
        lam2 = 2*Y*Z
        lam22 = lam2^2
        h1 = lam22 * (Q.Y*Z - QZY) - lam1*lam2 * (QXZ - QZX)
        h2 = lam22 * (QXZ + 2*QZX + A*QZZ) - lam12Z*Q.Z
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
    f1 = f1^2 * (Q.X * Z - Q.Z * X)
    f2 = f2^2 * Q.Z * Z
    return f1, f2
end

# Weil pairing e_{2^e}(P, Q)
function Weil_pairing_2power(A::T, P::Point{T}, Q::Point{T}, e::Integer) where T <: RingElem
    fPQ1, fPQ2 = Miller_function(A, P, Q, e)
    if fPQ1 == 0 || fPQ2 == 0
        return parent(A)(1)
    end
    fQP1, fQP2 = Miller_function(A, Q, P, e)
    return (fPQ1*fQP2) / (fPQ2*fQP1)
end

# precomputed table for Tate pairings
function make_pairing_table(A::FqFieldElem, P::Point{FqFieldElem}, e::Integer)
    R = P
    x, y = R.X/R.Z, R.Y/R.Z
    table = [[x, y, parent(A)(0)]]
    for i in 1:e-1
        lam = (3*x^2 + 2*A*x + 1) / (2*y)
        R = double(R, Proj1(A))
        x = R.X/R.Z
        y = R.Y/R.Z
        push!(table, [x, y, lam])
    end
    return table
end

# Tate pairing t_{2^e}(P0, P) using precomputed table for P0
function Tate_pairing_P0(P::Point{FqFieldElem}, table::Vector{Vector{FqFieldElem}}, f::Integer)
    x, y, z = P.X, P.Y, P.Z
    x_frob = Frob(x)
    z_frob = Frob(z)
    x0, y0 = table[1][1], table[1][2]
    f0 = 1
    h0 = 1
    for (xt, yt, lam) in table[2:end]
        t0 = x - x0 * z
        t1 = y - y0 * z
        t0 *= lam
        g = t0 - t1
        h = x_frob - Frob(xt) * z_frob
        g *= h
        f0 = f0^2 * g
        h0 = h0^2 * z * z_frob
        x0, y0 = xt, yt
    end
    g = x - x0 * z
    f0 = f0^2 * g
    h0 = h0^2 * z
    f0 = Frob(f0) * h0 / (f0 * Frob(h0))
    f0 = f0^f
    return f0
end

# Tate pairing t_{2^e}(iP0, P) using precomputed table for P0
function Tate_pairing_iP0(P::Point{FqFieldElem}, table::Vector{Vector{FqFieldElem}}, f::Integer)
    x, y, z = P.X, P.Y, P.Z
    x_frob = Frob(x)
    z_frob = Frob(z)
    x0, y0 = table[1][1], table[1][2]
    f0 = 1
    h0 = 1
    for (xt, yt, lam) in table[2:end]
        t0 = x + x0 * z
        t1 = y - mult_by_i(y0) * z
        t0 *= -mult_by_i(lam)
        g = t0 - t1
        h = x_frob + Frob(xt) * z_frob
        g *= h
        f0 = f0^2 * g
        h0 = h0^2 * z * z_frob
        x0, y0 = xt, yt
    end
    g = x + x0 * z
    f0 = f0^2 * g
    h0 = h0^2 * z
    f0 = Frob(f0) * h0 / (f0 * Frob(h0))
    f0 = f0^f
    return f0
end
