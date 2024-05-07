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

# return n1, n2 such that P = [n1]P0 + [n2]Q0
function ec_bi_dlog_challenge(A::T, xP::Proj1{T}, P0::Point{T}, Q0::Point{T}, E0::E0Data) where T <: RingElem
    P = Point(A, xP)
    base = Weil_pairing_2power(A, P0, Q0, SQISIGN_challenge_length)
    w1 = Weil_pairing_2power(A, P, Q0, SQISIGN_challenge_length)
    w2 = Weil_pairing_2power(A, P0, P, SQISIGN_challenge_length)

    n0 = fq_dlog_power_of_2_opt(base, E0.dlog_data_chall)
    n1 = fq_dlog_power_of_2_opt(w1, E0.dlog_data_chall)
    n2 = fq_dlog_power_of_2_opt(w2, E0.dlog_data_chall)
    n0inv = invmod(n0, BigInt(2)^SQISIGN_challenge_length)

    return (n1 * n0inv) % BigInt(2)^SQISIGN_challenge_length, (n2 * n0inv) % BigInt(2)^SQISIGN_challenge_length
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

# return n1, n2, n3, n4 such that P = [n1]P0 + [n1]Q0, Q = [n3]P0 + [n4]Q0, where (P0, Q0) is a fixed basis of E0d[2^ExponentFull]
function ec_bi_dlog_E0d(P::Point{FqFieldElem}, Q::Point{FqFieldElem}, global_data::GlobalData, idx::Int)
    E0 = global_data.E0
    E0d = global_data.orders_data[idx]
    t1 = Tate_pairing_P0(P, E0d.tate_tableQ, Cofactor)
    t2 = Tate_pairing_P0(P, E0d.tate_tableP, Cofactor)
    t3 = Tate_pairing_P0(Q, E0d.tate_tableQ, Cofactor)
    t4 = Tate_pairing_P0(Q, E0d.tate_tableP, Cofactor)

    n1 = -fq_dlog_power_of_2_opt(t1, E0.dlog_data_full) * E0d.dlog_base
    n2 = fq_dlog_power_of_2_opt(t2, E0.dlog_data_full) * E0d.dlog_base
    n3 = -fq_dlog_power_of_2_opt(t3, E0.dlog_data_full) * E0d.dlog_base
    n4 = fq_dlog_power_of_2_opt(t4, E0.dlog_data_full) * E0d.dlog_base

    return n1, n2, n3, n4
end

# return n1, n2, n3, n4 such that P = [n1]P0 + [n1]Q0, Q = [n3]P0 + [n4]Q0, where (P0, Q0) is a fixed basis of E0d[2^ExponentFull]
function ec_bi_dlog_E0d(xP::Proj1{FqFieldElem}, xQ::Proj1{FqFieldElem}, xPQ::Proj1{FqFieldElem}, global_data::GlobalData, idx::Int)
    E0d = global_data.orders_data[idx]
    A0 = E0d.A
    P = Point(A0, xP)
    Q = Point(A0, xQ)
    PQ = add(P, -Q, Proj1(A0))
    if !(xPQ == Proj1(PQ.X, PQ.Z))
        Q = -Q
    end

    return ec_bi_dlog_E0d(P, Q, global_data, idx)
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