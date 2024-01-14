using Test
using Nemo
import SQIsign2D: Proj1, IsInfinity, Affine,
    xDBL, xADD, xDBLADD, Ladder, random_point, Montgomery_coeff,
    two_e_iso

# order counting
function order_count(A::FinFieldElem)
    F = parent(A)
    cnt = 1
    for x in F
        y2 = x^3 + A*x^2 + x
        if y2 == 0
            cnt += 1
        elseif is_square(y2)
            cnt += 2
        end
    end
    return cnt
end

# test xDBL, xADD, and xDBLADD
function test_xDBL_xADD_xDBLADD(F::FinField)
    A = rand(F)
    while A^2 == 4
        A = rand(F)
    end
    a24 = Proj1(A + 2, F(4))
    P = Proj1(rand(F), rand(F))
    P2 = xDBL(P, a24)
    P3 = xADD(P2, P, P)
    P2d, P3d = xDBLADD(P, P2, P, a24)
    @test P2 == P2d
    @test P3 == P3d
end

# test Ladder
function test_Ladder(F::FinField)
    A = rand(F)
    while A^2 == 4
        A = rand(F)
    end
    a24 = Proj1(A + 2, F(4))
    P = random_point(A)
    ord = order_count(A)
    @test IsInfinity(Ladder(ord, P, a24))
end

# test isogenies
function test_isogenies(F::FinField)
    A = rand(F)
    q = order(F)
    while A^2 == 4 || !is_square(A^2 - 4)
        A = rand(F)
    end
    a24 = Proj1(A + 2, F(4))
    P = random_point(A)
    ord = order_count(A)
    N = ord
    while N % 2 == 0
        N = div(N, 2)
    end
    K = Ladder(N, P, a24)
    e = 0
    P = K
    while !IsInfinity(P)
        P = xDBL(P, a24)
        e += 1
    end
    @assert IsInfinity(Ladder(2^e, K, a24))
    Q = random_point(A)
    Q2 = xDBL(Q, a24)
    a24d, (Qd, Q2d) = two_e_iso(a24, K, e, [Q, Q2])
    t = ord - (q + 1)
    td = order_count(Montgomery_coeff(a24d)) - (q + 1)
    @test abs(t) == abs(td)
    @test xDBL(Qd, a24d) == Q2d
end

for _ in 1:10
    p = 1
    while !is_prime(p)
        p = rand(3:1000)
    end
    F, _ = finite_field(p, 2)
    test_xDBL_xADD_xDBLADD(F)
    test_Ladder(F)
    test_isogenies(F)
end
