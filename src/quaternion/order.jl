export QOrderElem, involution, norm, quadratic_form, order_elem_from_standard_basis

# element in the maximal order <1, i, (i + j)/2, (1 + ij)/2> in B_{p, infinity}
# represented by the coefficients of the above basis
struct QOrderElem
    a::BigInt
    b::BigInt
    c::BigInt
    d::BigInt
end

const Div_p1_4 = div(p + 1, 4)
const Quaternion_0 = QOrderElem(0, 0, 0, 0)
const Quaternion_1 = QOrderElem(1, 0, 0, 0)
const Quaternion_i = QOrderElem(0, 1, 0, 0)
const Quaternion_j = QOrderElem(0, -1, 2, 0)
const Quaternion_ij = QOrderElem(-1, 0, 0, 2)

function QOrderElem(v::Vector{T}) where T <: Integer
    return QOrderElem(v[1], v[2], v[3], v[4])
end

function QOrderElem(a::Integer)
    return QOrderElem(a, 0, 0, 0)
end

function Base.getindex(x::QOrderElem, i::Integer)
    if i == 1
        return x.a
    elseif i == 2
        return x.b
    elseif i == 3
        return x.c
    elseif i == 4
        return x.d
    else
        throw(BoundsError(x, i))
    end
end

function Base.:(==)(x::QOrderElem, y::QOrderElem)
    return x.a == y.a && x.b == y.b && x.c == y.c && x.d == y.d
end

function Base.:(==)(x::QOrderElem, a::Integer)
    return x.a == a && x.b == 0 && x.c == 0 && x.d == 0
end

function Base.:+(x::QOrderElem, y::QOrderElem)
    return QOrderElem(x.a + y.a, x.b + y.b, x.c + y.c, x.d + y.d)
end

function Base.:-(x::QOrderElem, y::QOrderElem)
    return QOrderElem(x.a - y.a, x.b - y.b, x.c - y.c, x.d - y.d)
end

function Base.:-(x::QOrderElem)
    return QOrderElem(-x.a, -x.b, -x.c, -x.d)
end

function Base.div(x::QOrderElem, a::Integer)
    return QOrderElem(div(x.a, a), div(x.b, a), div(x.c, a), div(x.d, a))
end

function Base.gcd(x::QOrderElem)
    return gcd(x.a, x.b, x.c, x.d)
end

function to_vector(x::QOrderElem)
    return [x.a, x.b, x.c, x.d]
end

function order_elem_from_standard_basis(a::Integer, b::Integer, c::Integer, d::Integer)
    return QOrderElem(a - d, b - c, 2*c, 2*d)
end

function left_mult_matrix(x::QOrderElem)
    return [x.a -x.b -x.b-Div_p1_4*x.c -Div_p1_4*x.d
            x.b x.a -Div_p1_4*x.d x.b+Div_p1_4*x.c
            x.c x.d x.a+x.d -x.b
            x.d -x.c x.b x.a+x.d]
end

function Base.:*(x::QOrderElem, y::QOrderElem)
    a, b, c, d = left_mult_matrix(x) * [y.a, y.b, y.c, y.d]
    return QOrderElem(a, b, c, d)
end

function Base.:*(a::Integer, x::QOrderElem)
    return QOrderElem(a*x.a, a*x.b, a*x.c, a*x.d)
end

function involution(x::QOrderElem)
    return QOrderElem(x.a+x.d, -x.b, -x.c, -x.d)
end

function norm(x::QOrderElem)
    return div((2*x.a + x.d)^2 + (2*x.b + x.c)^2 + p*(x.c^2 + x.d^2), 4)
end

function trace(x::QOrderElem)
    return 2*x.a + x.d
end

function quadratic_form(x::QOrderElem, y::QOrderElem)
    return norm(x + y) - norm(x) - norm(y)
end