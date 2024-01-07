using Test
using Nemo
import SQIsign2D: Montgomery, OddIsogeny, Ladder, PointOrder,
    Proj1, IsInfinity, Affine, DefintionField

# test Ladder
function order_check(Mont::Montgomery, n::Integer)
    A = Mont.A
    F = parent(A)
    X = rand(F)
    while !is_square(X^3 + A*X^2 + X)
        X = rand(F)
    end
    P = Proj1(X)
    res = Ladder(Mont, P, n)
    if !IsInfinity(res)
        println("F = ", F)
        println("A = ", A)
        println("n = ", n)
        println("P = ", P)
        println("nP = ", res)
    end
    return IsInfinity(res)
end

# test PointOrder
function point_order_check(Mont::Montgomery, l::Integer, ord::Integer)
    P = PointOrder(Mont, l, ord)
    lP = Ladder(Mont, P, l)
    return !IsInfinity(P) && IsInfinity(lP)
end

# order counting
function order_count(Mont::Montgomery)
    F = DefintionField(Mont)
    A = Mont.A
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

# test IsogCoeff
function isog_check(Mont::Montgomery, l::Integer)
    ord = order_count(Mont)
    P = PointOrder(Mont, l, ord)
    Md = OddIsogeny(Mont, P, l)
    return order_count(Md) == ord
end

for _ in 1:1000
    p = rand(5:10000)
    while !is_prime(p)
        p = rand(1:10000)
    end
    Fp = GF(p)
    A = rand(Fp)
    if A^2 == 4
        continue
    end
    Mont = Montgomery(A)
    ord = order_count(Mont)
    l = rand([f[1] for f in factor(ord)])

    @test order_check(Mont, ord) == true 
    @test point_order_check(Mont, l, ord) == true
    if l % 2 == 1
        @test isog_check(Mont, l) == true
    end
end