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
# Algorithm 9 in Adj & Rodríguez-Henríquez, "Square Root Computation over Even Extension Fields".
function square_root(x::FinFieldElem)
    Fp2 = parent(x)
    Fp = base_field(Fp2)
    p = characteristic(Fp)

    a1 = x^div(p - 3, 4)
    a1x = a1*x
    a = a1*a1x
    x0 = a1x
    if a == -1
        x0 = mult_by_i(x0)
    else
        b = (1 + a)^div(p - 1, 2)
        x0 *= b
    end

    return x0
end

# x < y in lexicographic order, i.e. x < y if and only if x[0] < y[0] or (x[0] == y[0] and x[1] < y[1])
function lex_order(x::FinFieldElem, y::FinFieldElem)
    coeff(x, 0) == coeff(y, 0) && return lift(ZZ, coeff(x, 1)) < lift(ZZ, coeff(y, 1))
    return lift(ZZ, coeff(x, 0)) < lift(ZZ, coeff(y, 0))
end