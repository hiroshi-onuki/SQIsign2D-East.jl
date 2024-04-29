export det_2x2, invmod_2x2

function det_2x2(M::Matrix{T}) where T <: Integer
    return M[1,1] * M[2,2] - M[1,2] * M[2,1]
end

function invmod_2x2(M::Matrix{T}, q::T) where T <: Integer
    return ([M[2, 2] -M[1, 2]; -M[2, 1] M[1, 1]] * invmod(det_2x2(M), q)) .% q
end