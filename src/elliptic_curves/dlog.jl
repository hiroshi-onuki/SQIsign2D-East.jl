# return n1, n2, n3, n4 such that P = [n1]P0 + [n1]Q0, Q = [n3]P0 + [n4]Q0, where (P0, Q0) is a fixed basis of E0[2^ExponentFull]
function ec_bi_dlog_E0(P::Point{FqFieldElem}, Q::Point{FqFieldElem}, E0::E0Data)
    t1 = Tate_pairing_iP0(P, E0.tate_table, Cofactor)
    t2 = Tate_pairing_P0(P, E0.tate_table, Cofactor)
    t3 = Tate_pairing_iP0(Q, E0.tate_table, Cofactor)
    t4 = Tate_pairing_P0(Q, E0.tate_table, Cofactor)

    n1 = -fq_dlog_power_of_2_opt(t1, E0.dlog_data_full)
    n2 = fq_dlog_power_of_2_opt(t2, E0.dlog_data_full)
    n3 = -fq_dlog_power_of_2_opt(t3, E0.dlog_data_full)
    n4 = fq_dlog_power_of_2_opt(t4, E0.dlog_data_full)

    return n1, n2, n3, n4
end

# return n1, n2, n3, n4 such that P = [n1]P0 + [n1]Q0, Q = [n3]P0 + [n4]Q0, where (P0, Q0) is a fixed basis of E0[2^ExponentFull]
function ec_bi_dlog_E0(xP::Proj1{T}, xQ::Proj1{T}, xPQ::Proj1{T}, E0::E0Data) where T <: RingElem
    A0 = E0.A0
    P = Point(A0, xP)
    Q = Point(A0, xQ)
    PQ = add(P, -Q, Proj1(A0))
    if !(xPQ == Proj1(PQ.X, PQ.Z))
        Q = -Q
    end

    return ec_bi_dlog_E0(P, Q, E0)
end

# return n1, n2, n3, n4 such that P = [n1]P0 + [n1]Q0, Q = [n3]P0 + [n4]Q0
# on E_A[2^SQISIGN_challenge_length]
function ec_bi_dlog_commitment(A::T, xP::Proj1{T}, xQ::Proj1{T}, xPQ::Proj1{T}, 
                    xPb::Proj1{T}, xQb::Proj1{T}, xPQb::Proj1{T}, E0::E0Data) where T <: RingElem
    P = Point(A, xP)
    Q = Point(A, xQ)
    PQ = add(P, -Q, Proj1(A))
    if !(xPQ == Proj1(PQ.X, PQ.Z))
        Q = -Q
    end
    Pb = Point(A, xPb)
    Qb = Point(A, xQb)
    PQb = add(Pb, -Qb, Proj1(A))
    if !(xPQb == Proj1(PQb.X, PQb.Z))
        Qb = -Qb
    end

    base = Weil_pairing_2power(A, Pb, Qb, SQISIGN_challenge_length)
    w1 = Weil_pairing_2power(A, P, Qb, SQISIGN_challenge_length)
    w2 = Weil_pairing_2power(A, Pb, P, SQISIGN_challenge_length)
    w3 = Weil_pairing_2power(A, Q, Qb, SQISIGN_challenge_length)
    w4 = Weil_pairing_2power(A, Pb, Q, SQISIGN_challenge_length)

    n0 = fq_dlog_power_of_2_opt(base, E0.dlog_data_chall)
    n1 = fq_dlog_power_of_2_opt(w1, E0.dlog_data_chall)
    n2 = fq_dlog_power_of_2_opt(w2, E0.dlog_data_chall)
    n3 = fq_dlog_power_of_2_opt(w3, E0.dlog_data_chall)
    n4 = fq_dlog_power_of_2_opt(w4, E0.dlog_data_chall)
    n0inv = invmod(n0, BigInt(2)^SQISIGN_challenge_length)
    n1 = (n1 * n0inv) % BigInt(2)^SQISIGN_challenge_length
    n2 = (n2 * n0inv) % BigInt(2)^SQISIGN_challenge_length
    n3 = (n3 * n0inv) % BigInt(2)^SQISIGN_challenge_length
    n4 = (n4 * n0inv) % BigInt(2)^SQISIGN_challenge_length

    return n1, n2, n3, n4
end

# return n1, n2 such that P = [n1]P0 + [n2]Q0
function ec_bi_dlog_challenge(A::T, xP::Proj1{T}, xP0::Proj1{T}, xQ0::Proj1{T}, xPQ0::Proj1{T}, E0::E0Data) where T <: RingElem
    P = Point(A, xP)
    P0 = Point(A, xP0)
    Q0 = Point(A, xQ0)
    PQ0 = add(P0, -Q0, Proj1(A))
    if !(xPQ0 == Proj1(PQ0.X, PQ0.Z))
        Q0 = -Q0
    end
    base = Weil_pairing_2power(A, P0, Q0, SQISIGN_challenge_length)
    w1 = Weil_pairing_2power(A, P, Q0, SQISIGN_challenge_length)
    w2 = Weil_pairing_2power(A, P0, P, SQISIGN_challenge_length)

    n0 = fq_dlog_power_of_2_opt(base, E0.dlog_data_chall)
    n1 = fq_dlog_power_of_2_opt(w1, E0.dlog_data_chall)
    n2 = fq_dlog_power_of_2_opt(w2, E0.dlog_data_chall)
    n0inv = invmod(n0, BigInt(2)^SQISIGN_challenge_length)

    return (n1 * n0inv) % BigInt(2)^SQISIGN_challenge_length, (n2 * n0inv) % BigInt(2)^SQISIGN_challenge_length
end 

# return n1, n2, n3, n4 such that P = [n1]P0 + [n1]Q0, Q = [n3]P0 + [n4]Q0
# on E_A[2^ExponentForTorsion]
function ec_bi_dlog_response(A::T, xP::Proj1{T}, xQ::Proj1{T}, xPQ::Proj1{T}, 
    xPb::Proj1{T}, xQb::Proj1{T}, xPQb::Proj1{T}, E0::E0Data) where T <: RingElem
    P = Point(A, xP)
    Q = Point(A, xQ)
    PQ = add(P, -Q, Proj1(A))
    if !(xPQ == Proj1(PQ.X, PQ.Z))
    Q = -Q
    end
    Pb = Point(A, xPb)
    Qb = Point(A, xQb)
    PQb = add(Pb, -Qb, Proj1(A))
    if !(xPQb == Proj1(PQb.X, PQb.Z))
    Qb = -Qb
    end

    base = Weil_pairing_2power(A, Pb, Qb, ExponentForTorsion)
    w1 = Weil_pairing_2power(A, P, Qb, ExponentForTorsion)
    w2 = Weil_pairing_2power(A, Pb, P, ExponentForTorsion)
    w3 = Weil_pairing_2power(A, Q, Qb, ExponentForTorsion)
    w4 = Weil_pairing_2power(A, Pb, Q, ExponentForTorsion)

    n0 = fq_dlog_power_of_2_opt(base, E0.dlog_data_res)
    n1 = fq_dlog_power_of_2_opt(w1, E0.dlog_data_res)
    n2 = fq_dlog_power_of_2_opt(w2, E0.dlog_data_res)
    n3 = fq_dlog_power_of_2_opt(w3, E0.dlog_data_res)
    n4 = fq_dlog_power_of_2_opt(w4, E0.dlog_data_res)
    n0inv = invmod(n0, BigInt(2)^ExponentForTorsion)
    n1 = (n1 * n0inv) % BigInt(2)^ExponentForTorsion
    n2 = (n2 * n0inv) % BigInt(2)^ExponentForTorsion
    n3 = (n3 * n0inv) % BigInt(2)^ExponentForTorsion
    n4 = (n4 * n0inv) % BigInt(2)^ExponentForTorsion

    return n1, n2, n3, n4
end


# return n s.t. P = [n]Q, where P, Q is a fixed point of order 2^SQISIGN_challenge_length and R is a point s.t. (P, R) is a basis of E_A[2^SQISIGN_challenge_length]
function ec_dlog(A::T, xP::Proj1{T}, xQ::Proj1{T}, xR::Proj1{T}, E0::E0Data) where T <: RingElem
    P = Point(A, xP)
    Q = Point(A, xQ)
    R = Point(A, xR)
    w1 = Weil_pairing_2power(A, P, R, SQISIGN_challenge_length)
    w2 = Weil_pairing_2power(A, Q, R, SQISIGN_challenge_length)
    n1 = fq_dlog_power_of_2_opt(w1, E0.dlog_data_chall)
    n2 = fq_dlog_power_of_2_opt(w2, E0.dlog_data_chall)

    return (n1 * invmod(n2, BigInt(2)^SQISIGN_challenge_length)) % BigInt(2)^SQISIGN_challenge_length
end

# x^(2^e)
function square_e(x::FqFieldElem, e::Int)
    y = x
    for i in 1:e
        y = y^2
    end
    return y
end

# return n such that x = base^e
function fq_dlog_power_of_2(x::FqFieldElem, base::FqFieldElem, e::Integer)
    n = BigInt(0)
    t = x
    for i in 1:e
        if t^(BigInt(2)^(e-i)) == base^(BigInt(2)^(e-1))
            n += BigInt(2)^(i-1)
            t //= square_e(base, i-1)
        end
    end
    return n
end

# make a precomputed table for dlog with base of order 2^e
function make_dlog_table(base::FqFieldElem, e::Int, window_size::Int)
    F = parent(base)
    l = 2^window_size
    f, r = divrem(e, window_size)
    T1 = [[F(1) for _ in 1:l] for _ in 1:(f+1)]
    T2 = [[F(1) for _ in 1:l] for _ in 1:f]

    T1[1][2] = 1/base
    for j in 2:l-1
        T1[1][j+1] = (T1[1][2])^j
    end
    for i in 1:f
        for j in 1:l-1
            T1[i+1][j+1] = square_e(T1[i][j+1], window_size)
        end
    end

    for i in 1:f
        for j in 2:l
            T2[i][j] = square_e(T1[i][j], r)
        end
    end

    return T1, T2
end

# compute x=[x0,x1,...,x{k-1}] s.t. h = (g^(2^r*l^(e-k)))^x
function fq_dlog_subtree(e::Int, h::FqFieldElem, window_size::Int,
        strategy::Vector{Int}, table::Vector{Vector{FqFieldElem}})
    t = length(strategy)
    l = BigInt(2^window_size)
    if t == 0
        h == 1 && return [0]
        for j in 1:l-1
            if h == table[end][j+1]
                return [l - j]
            end
        end 
    end
    n = strategy[1]
    L = strategy[2:t-n+1]
    R = strategy[t-n+2:t]

    hL = square_e(h, n * window_size)
    xL = fq_dlog_subtree(e - n, hL, window_size, L, table)

    hR = h
    for i in 1:e-n
        hR *= table[end-e+i][xL[i]+1]
    end
    xR = fq_dlog_subtree(n, hR, window_size, R, table)

    return vcat(xL, xR)
end

# return n such that h = g^n, where g is a fixed base of order 2^e
function fq_dlog_power_of_2_opt(h::FqFieldElem, dlog_data::DlogData)
    e, window_size, T1, T2, strategy = dlog_data.e, dlog_data.window_size, dlog_data.T1, dlog_data.T2, dlog_data.strategy
    l = BigInt(2^window_size)
    f, r = divrem(e, window_size)
    xw = fq_dlog_subtree(f, square_e(h, r), window_size, strategy, T2)

    hr = h
    for i in 1:f
        hr *= T1[i][xw[i] + 1]
    end
    xr = 0
    for j in 0:2^r-1
        if hr == T1[end][2^r-j+1]
            xr = j
            break
        end
    end

    x = BigInt(0)
    for i in 0:f-1
        x += xw[i+1] * l^i
    end
    x += xr * l^f
    return x
end

function bi_dlog_odd_prime(A::T, P::Point{T}, R::Point{T}, S::Point{T}, l::Int) where T <: RingElem
    Pd = infinity_full_point(parent(A))
    for a in 0:(l-1)
        Pdd = Pd
        for b in 0:(l-1)
            if P == Pdd
                return a, b
            end
            Pdd = add(Pdd, S, Proj1(A))
        end
        Pd = add(Pd, R, Proj1(A))
    end
    @assert false "bi_dlog_odd_prime: no solution"
end

function bi_dlog_odd_prime_power(A::T, P::Point{T}, R::Point{T}, S::Point{T}, l::Int, e::Int) where T <: RingElem
    e == 1 && return bi_dlog_odd_prime(A, P, R, S, l)
    f = div(e, 2)
    Pd = mult(l^(e - f), P, Proj1(A))
    Rd = mult(l^(e - f), R, Proj1(A))
    Sd = mult(l^(e - f), S, Proj1(A))
    a, b = bi_dlog_odd_prime_power(A, Pd, Rd, Sd, l, f)
    aRbS = add(mult(a, R, Proj1(A)), mult(b, S, Proj1(A)), Proj1(A))
    P = add(P, -aRbS, Proj1(A))
    R = mult(l^f, R, Proj1(A))
    S = mult(l^f, S, Proj1(A))
    c, d = bi_dlog_odd_prime_power(A, P, R, S, l, e - f)
    return a + c * l^f, b + d * l^f
end

function bi_dlog_odd_prime_power(A::T, xP::Proj1{T}, xQ::Proj1{T}, xPQ::Proj1{T}, xR::Proj1{T}, xS::Proj1{T}, xRS::Proj1{T}, l::Int, e::Int) where T <: RingElem
    P = Point(A, xP)
    Q = Point(A, xQ)
    PQ = add(P, -Q, Proj1(A))
    if !(xPQ == Proj1(PQ.X, PQ.Z))
        Q = -Q
    end
    R = Point(A, xR)
    S = Point(A, xS)
    RS = add(R, -S, Proj1(A))
    if !(xRS == Proj1(RS.X, RS.Z))
        S = -S
    end
    a, b = bi_dlog_odd_prime_power(A, P, R, S, l, e)
    c, d = bi_dlog_odd_prime_power(A, Q, R, S, l, e)
    return a, b, c, d
end
