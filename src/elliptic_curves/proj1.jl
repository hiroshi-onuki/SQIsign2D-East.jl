export Proj1, infinity_point, is_infinity, affine

# projective point of dimension one
struct Proj1{T <: RingElem}
    X::T
    Z::T
end

Proj1(x::RingElem) = Proj1(x, parent(x)(1))
infinity_point(F::T) where T <: Ring = Proj1(F(1), F(0))
is_infinity(P::Proj1) = P.Z == 0
affine(P::Proj1{T}) where T <: FieldElem = P.Z == 0 ? Nothing : P.X//P.Z
Base.:-(P::Proj1) = Proj1(-P.X, P.Z)
Base.:(==)(P::Proj1, Q::Proj1) = P.X*Q.Z == P.Z*Q.X
Base.:(copy)(P::Proj1) = Proj1(P.X, P.Z)

