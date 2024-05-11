
function quaternion_to_matrix(alpha::QOrderElem, E0_data::E0Data)
    return alpha[1] * [1 0; 0 1] + alpha[2] * E0_data.Matrices_2e[1] + alpha[3] * E0_data.Matrices_2e[2] + alpha[4] * E0_data.Matrices_2e[3]
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
    M = quaternion_to_matrix(alpha, E0_data)
    return action_of_matrix(M, a24, xP, xQ, xPQ)
end

# return the codomain of a random d-isogeny from E0 and the images of the basis points
function RandIsogImages(d::BigInt, global_data::GlobalData, compute_odd_points::Bool=false)
    deg_dim2 = BigInt(1) << ExponentFull
    E0_data = global_data.E0_data
    a24_0 = E0_data.a24_0
    xP0, xQ0, xPQ0 = E0_data.xP2e, E0_data.xQ2e, E0_data.xPQ2e

    alpha, _ = FullRepresentInteger(d*(deg_dim2 - d))

    a24, xP, xQ, xPQ, odd_images = d2isogeny_form_Esquare(a24_0, d, alpha, xP0, xQ0, xPQ0, global_data, compute_odd_points)
    if compute_odd_points
        return a24, xP, xQ, xPQ, odd_images, LeftIdeal(alpha, d)
    else
        return a24, xP, xQ, xPQ, LeftIdeal(alpha, d)
    end
end

# return the codomain of a random d-isogeny from E and the images of (P, Q),
# where P, Q is the image of the fixed basis of E0[2^ExponentFull] under an isogeny corresponding to I
function GeneralizedRandomIsogImages(d::BigInt, a24::Proj1{T}, xP::Proj1{T}, xQ::Proj1{T}, xPQ::Proj1{T},
            I::LeftIdeal, nI::BigInt, global_data::GlobalData) where T <: RingElem
    N = d*((BigInt(1) << ExponentFull) - d)
    alpha = Quaternion_0
    
    found = false

    # make alpha in I + Z s.t. n(alpha) = N
    C, D = EichlerModConstraint(I, nI, Quaternion_1, Quaternion_1, false)
    N_CD = p * (C^2 + D^2)
    N_N_CD = (N * invmod(N_CD, nI)) % nI
    lambda = sqrt_mod(4*N_N_CD, nI)
    alpha, found = FullStrongApproximation(nI, C, D, lambda, 4*N, KLPT_signing_number_strong_approx)
    @assert found
    @assert norm(alpha) == N

    a24, xP, xQ, xPQ, _ = d2isogeny_form_Esquare(a24, d, alpha, xP, xQ, xPQ, global_data)

    return a24, xP, xQ, xPQ
end