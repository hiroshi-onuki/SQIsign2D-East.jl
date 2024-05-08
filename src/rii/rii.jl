
function quaternion_to_matrix(alpha::QOrderElem, E0_data::E0Data)
    return alpha[1] * [1 0; 0 1] + alpha[2] * E0_data.Matrices_2e[1] + alpha[3] * E0_data.Matrices_2e[2] + alpha[4] * E0_data.Matrices_2e[3]
end

function action_of_matrix(M::Matrix{BigInt}, E0_data::E0Data, Ponly::Bool=false)
    xP, xQ, xPQ = E0_data.xP2e, E0_data.xQ2e, E0_data.xPQ2e
    xP_new = linear_comb_2_e(M[1, 1], M[2, 1], xP, xQ, xPQ, E0_data.a24_0, ExponentFull)
    Ponly && return xP_new
    xQ_new = linear_comb_2_e(M[1, 2], M[2, 2], xP, xQ, xPQ, E0_data.a24_0, ExponentFull)
    xPQ_new = linear_comb_2_e(M[1, 1] - M[1, 2], M[2, 1] - M[2, 2], xP, xQ, xPQ, E0_data.a24_0, ExponentFull)
    return xP_new, xQ_new, xPQ_new
end

# alpha(P), alpha(Q) for a fixed basis (P, Q) of E0[2^ExponentFull]
function action_on_torsion_basis(alpha::QOrderElem, E0_data::E0Data)
    M = quaternion_to_matrix(alpha, E0_data)
    return action_of_matrix(M, E0_data)
end

# return the codomain of a random d-isogeny from E0 and the images of the basis points
function RandIsogImages(d::BigInt, E0_data::E0Data, output_ideal::Bool=false)
    deg_dim2 = BigInt(1) << ExponentFull
    a24_0 = E0_data.a24_0
    xP0, xQ0, xPQ0 = E0_data.xP2e, E0_data.xQ2e, E0_data.xPQ2e

    alpha, _ = FullRepresentInteger(d*(deg_dim2 - d))

    if (alpha + deg_dim2) % 2 == 1
        # the first (2,2)-isogeny is E0^2 -> E0^2 represented by the matrix [1 1; 1 -1]
        xP1, xQ1, xPQ1 = action_on_torsion_basis(alpha - d, E0_data)
        xP2, xQ2, xPQ2 = action_on_torsion_basis(-alpha - d, E0_data)
        xR1, xS1, xRS1 = xP0, xQ0, xPQ0
        xR2, xS2, xRS2 = xP0, xQ0, xPQ0
        n = norm(alpha + d)
        beta1 = alpha - d
        beta2 = -alpha - d
        gamma = Quaternion_1
    elseif (alpha + deg_dim2*Quaternion_i) % 2 == Quaternion_i
        # the first (2,2)-isogeny is E0^2 -> E0^2 represented by the matrix [1 i; i 1]
        xP1, xQ1, xPQ1 = action_on_torsion_basis(Quaternion_i*alpha - d, E0_data)
        xP2, xQ2, xPQ2 = action_on_torsion_basis(alpha - d*Quaternion_i, E0_data)
        xR1, xS1, xRS1 = xP0, xQ0, xPQ0
        xR2, xS2, xRS2 = -xP0, -xQ0, -xPQ0
        n = norm(Quaternion_i * alpha - d)
        beta1 = Quaternion_i * alpha - d
        beta2 = alpha - d * Quaternion_i
        gamma = Quaternion_i
    else
        xP1 = ladder(deg_dim2 - d, xP0, a24_0)
        xQ1 = ladder(deg_dim2 - d, xQ0, a24_0)
        xPQ1 = ladder(deg_dim2 - d, xPQ0, a24_0)
        xP2, xQ2, xPQ2 = action_on_torsion_basis(alpha, E0_data)
        O0 = infinity_point(parent(a24_0.X))
        xR1, xS1, xRS1 = xP0, xQ0, xPQ0
        xR2, xS2, xRS2 = O0, O0, O0
        n = 1
        beta1 = -d * Quaternion_1
        beta2 = alpha
        gamma = Quaternion_0
    end

    a24_1= a24_0
    a24_2 = a24_0
    e = -1
    while n % 2 == 0
        n >>= 1
        e += 1
    end
    e = max(e, 0)

    # compute R + T etc. for the gluing isogeny
    M = quaternion_to_matrix((BigInt(1) << (ExponentFull - e - 2)) * beta1, E0_data)
    xR1_T = action_of_matrix([1 0; 0 0] + M, E0_data, true)
    xS1_T = action_of_matrix([0 0; 1 0] + M, E0_data, true)
    xRS1_T = action_of_matrix([1 0; -1 0] + M, E0_data, true)
    M = quaternion_to_matrix((BigInt(1) << (ExponentFull - e - 2)) * beta2, E0_data)
    Mg = quaternion_to_matrix(gamma, E0_data)
    xR2_T = action_of_matrix(Mg*[1 0; 0 0] + M, E0_data, true)
    xS2_T = action_of_matrix(Mg*[0 0; 1 0] + M, E0_data, true)
    xRS2_T = action_of_matrix(Mg*[1 0; -1 0] + M, E0_data, true)

    # if e > 1 then we compute (2, 2)-isogenies by Velu's formulas
    if e > 1
        if is_infinity(xDBLe(xP1, a24_0, ExponentFull-2))
            K = xDBLe(xQ1, a24_0, ExponentFull-e)
        else
            K = xDBLe(xP1, a24_0, ExponentFull-e)
        end
        a24_1, (xP1, xQ1, xPQ1, xR1, xS1, xRS1, xR1_T, xS1_T, xRS1_T) = two_e_iso(a24_0, K, e-1, [xP1, xQ1, xPQ1, xR1, xS1, xRS1, xR1_T, xS1_T, xRS1_T])
        if is_infinity(xDBLe(xP2, a24_0, ExponentFull-2))
            K = xDBLe(xQ2, a24_0, ExponentFull-e)
        else
            K = xDBLe(xP2, a24_0, ExponentFull-e)
        end
        a24_2, (xP2, xQ2, xPQ2, xR2, xS2, xRS2, xR2_T, xS2_T, xRS2_T) = two_e_iso(a24_0, K, e-1, [xP2, xQ2, xPQ2, xR2, xS2, xRS2, xR2_T, xS2_T, xRS2_T])
    end

    P1P2 = CouplePoint(xP1, xP2)
    Q1Q2 = CouplePoint(xQ1, xQ2)
    PQ1PQ2 = CouplePoint(xPQ1, xPQ2)
    R1R2 = CouplePoint(xR1, xR2)
    S1S2 = CouplePoint(xS1, xS2)
    RS1RS2 = CouplePoint(xRS1, xRS2)
    R1R2_T = CouplePoint(xR1_T, xR2_T)
    S1S2_T = CouplePoint(xS1_T, xS2_T)
    RS1RS2_T = CouplePoint(xRS1_T, xRS2_T)

    if haskey(StrategiesDim2, ExponentFull - e)
        strategy = StrategiesDim2[ExponentFull - e]
    else
        strategy = compute_strategy(ExponentFull - e - 2, 2, 1)
    end
    Es, images = product_isogeny_sqrt(a24_1, a24_2, P1P2, Q1Q2, PQ1PQ2, [R1R2, S1S2, RS1RS2], [R1R2_T, S1S2_T, RS1RS2_T], ExponentFull - e, strategy)

    xP, xQ, xPQ = images[1][1], images[2][1], images[3][1]
    A = Es[1]
    w0 = E0_data.Weil_P2eQ2e
    w1 = Weil_pairing_2power(affine(A), xP, xQ, xPQ, ExponentFull)
    if w1 != w0^d
        xP, xQ, xPQ = images[1][2], images[2][2], images[3][2]
        A = Es[2]
    end
    if output_ideal
        return A, xP, xQ, xPQ, LeftIdeal(alpha, d)
    else
        return A, xP, xQ, xPQ
    end
end

function GeneralizedRandomIsogImages(d::BigInt, a24::Proj1{T}, I::LeftIdeal, nI::BigInt, E0_data::E0Data) where T <: RingElem
    N = d*((BigInt(1) << ExponentFull) - d)
    alpha = Quaternion_0
    found = false
    while !found
        C, D, N_N_CD = 0, 0, 0
        while true
            C, D = EichlerModConstraint(I, nI, Quaternion_1, Quaternion_1, false)
            N_CD = p * (C^2 + D^2)
            return quadratic_residue_symbol(N_CD, nI)
            N_N_CD = (N * invmod(N_CD, nI)) % nI
            quadratic_residue_symbol(N_N_CD, nI) == 1 && break
        end
        lambda = sqrt_mod(4*N_N_CD, nI)

        alpha, found = FullStrongApproximation(nI, C, D, lambda, 4*N, KLPT_signing_number_strong_approx)
    end
    println("alpha: ", alpha, " found: ", found)
end