export CouplePoint, double, double_iter

struct CouplePoint{T <: RingElem}
    P1::Proj1{T}
    P2::Proj1{T}
end

function Base.getindex(P::CouplePoint{T}, i::Integer) where T <: RingElem
    if i == 1
        return P.P1
    elseif i == 2
        return P.P2
    else
        throw(BoundError(P, i))
    end
end

function double(P::CouplePoint{T}, a24_1::Proj1{T}, a24_2::Proj1{T}) where T <: RingElem
    P1, P2 = P.P1, P.P2
    P1 = xDBL(P1, a24_1)
    P2 = xDBL(P2, a24_2)
    return CouplePoint(P1, P2)
end

function double_iter(P::CouplePoint{T}, a24_1::Proj1{T}, a24_2::Proj1{T}, e::Integer) where T <: RingElem
    for _ in 1:e
        P = double(P, a24_1, a24_2)
    end
    return P
end