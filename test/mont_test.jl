using Test
using Nemo
import SQIsign2D: Proj1, xDBL, xADD, xDBLADD, Ladder, IsInfinity, Affine

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
    P = Proj1(rand(F), rand(F))
    while !is_square(Affine(P)^3 + A*Affine(P)^2 + Affine(P))
        P = Proj1(rand(F), rand(F))
    end
    ord = order_count(A)
    @test IsInfinity(Ladder(ord, P, a24))
end

for _ in 1:100
    p = 1
    while !is_prime(p)
        p = rand(3:1000)
    end
    F, _ = finite_field(p, 2)
    test_xDBL_xADD_xDBLADD(F)
    test_Ladder(F)
end
