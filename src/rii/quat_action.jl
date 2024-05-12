# the matrix representation of alpha by Ms
function quaternion_to_matrix(alpha::QOrderElem, Ms::Vector{Matrix{BigInt}})
    return alpha[1] * [1 0; 0 1] + alpha[2] * Ms[1] + alpha[3] * Ms[2] + alpha[4] * Ms[3]
end

# compute M[1, 1]*P + M[2, 1]*Q, M[1, 2]*P + M[2, 2]*Q for a basis (P, Q) of E[2^ExponentFull]
function action_of_matrix(M::Matrix{BigInt}, a24::Proj1{T}, xP::Proj1{T}, xQ::Proj1{T}, xPQ::Proj1{T}, Ponly::Bool=false) where T <: RingElem
    xP_new = linear_comb_2_e(M[1, 1], M[2, 1], xP, xQ, xPQ, a24, ExponentFull)
    Ponly && return xP_new
    xQ_new = linear_comb_2_e(M[1, 2], M[2, 2], xP, xQ, xPQ, a24, ExponentFull)
    xPQ_new = linear_comb_2_e(M[1, 1] - M[1, 2], M[2, 1] - M[2, 2], xP, xQ, xPQ, a24, ExponentFull)
    return xP_new, xQ_new, xPQ_new
end

# alpha(P), alpha(Q) for a fixed basis (P, Q) of E[2^ExponentFull]
function action_on_torsion_basis(alpha::QOrderElem, a24::Proj1{T}, xP::Proj1{T}, xQ::Proj1{T}, xPQ::Proj1{T}, E0_data::E0Data) where T <: RingElem
    M = quaternion_to_matrix(alpha, E0_data.Matrices_2e)
    return action_of_matrix(M, a24, xP, xQ, xPQ)
end
