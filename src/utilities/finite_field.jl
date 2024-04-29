export Frob, mult_by_i, square_root

# p-th power Frobenius map
function Frob(x::FinFieldElem)
    return coeff(x, 0) - coeff(x, 1)*gen(parent(x))
end

# multiplication by i in Fp2
function mult_by_i(x::FinFieldElem)
    return -coeff(x, 1) + coeff(x, 0)*gen(parent(x))
end

# square root of x in Fp2
function square_root(x::FinFieldElem)
    Fp2 = parent(x)
    Fp = base_field(Fp2)
    i = gen(Fp2)
    p = characteristic(Fp)
    inv2 = div(p + 1, 2)
    a, b = Fp(coeff(x, 0)), Fp(coeff(x, 1))
    if b == 0
        if a^div(p - 1, 2) == 1
            return [a^div(p + 1, 4), Fp2(1)]
        else
            return [(-a)^div(p + 1, 4)*i, Fp2(1)]
        end
    end
    d = (a^2 + b^2)^div(p + 1, 4)
    t = (a + d) * inv2
    x = t^div(p + 1, 4)
    if x^2 != t
        x = ((a - d) * inv2)^div(p + 1, 4)
    end
    y = b*inv2
    return [x^2 + y*i, x]
end

# x < y in lexicographic order, i.e. x < y if and only if x[0] < y[0] or (x[0] == y[0] and x[1] < y[1])
function lex_order(x::FinFieldElem, y::FinFieldElem)
    coeff(x, 0) == coeff(y, 0) && return lift(ZZ, coeff(x, 1)) < lift(ZZ, coeff(y, 1))
    return lift(ZZ, coeff(x, 0)) < lift(ZZ, coeff(y, 0))
end