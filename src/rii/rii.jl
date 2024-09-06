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
    l = FactorForAuxiliaryDegree    # we assume l = 3
    extra_path = quadratic_residue_symbol(-N, nI) != 1
    if extra_path
        N = div(N, l)
    end
    alpha = Quaternion_0

    # make alpha in I + Z s.t. n(alpha) = N
    C, D = EichlerModConstraint(I, nI, Quaternion_1, Quaternion_1, false)
    N_CD = p * (C^2 + D^2)
    N_N_CD = (N * invmod(N_CD, nI)) % nI
    lambda = sqrt_mod(4*N_N_CD, nI)
    tries = KLPT_keygen_number_strong_approx
    found = false
    for _ in 1:10
        alpha, found = FullStrongApproximation(nI, C, D, lambda, 4*N, tries)
        found && break
        tries *= 2
    end
    @assert found

    if extra_path
        println("extra path")
        d_inv = invmod(d, BigInt(1) << ExponentFull)

        M = quaternion_to_matrix(alpha, global_data.E0_data.Matrices_2e)
        s = 0
        if (M - [1 0; 0 1]) .% 4 == [0 0; 0 0]
            s = 1
        elseif (M + [1 0; 0 1]) .% 4 == [0 0; 0 0]
            s = -1
        end

        if s != 0
            # compute an extra l-isogeny
            K = random_point_order_l(a24, p + 1, l)
            a24d, (xPd, xQd, xPQd) = odd_isogeny(a24, K, l, [xP, xQ, xPQ])

            xPtmp, xQtmp, xPQtmp = action_on_torsion_basis(d_inv * alpha, a24d, xPd, xQd, xPQd, global_data.E0_data)
            w1 = Weil_pairing_2power(Montgomery_coeff(a24), xP, xQ, xPQ, ExponentFull)
            w2 = Weil_pairing_2power(Montgomery_coeff(a24d), xPtmp, xQtmp, xPQtmp, ExponentFull)
            @assert w1 * w2 == 1
            
            # the first (4,4)-isogeny is E * E' -> E * E' represented by the matrix [-s hat(phi); s*phi 1]
            beta1 = -s + d_inv * l * alpha
            beta2 = s + d_inv * alpha
            M1 = quaternion_to_matrix(beta1, global_data.E0_data.Matrices_2e)
            M2 = quaternion_to_matrix(beta2, global_data.E0_data.Matrices_2e)
            xP1, xQ1, xPQ1 = action_of_matrix(M1, a24, xP, xQ, xPQ)
            xP2, xQ2, xPQ2 = action_of_matrix(M2, a24d, xPd, xQd, xPQd)

            w1 = Weil_pairing_2power(Montgomery_coeff(a24), xP1, xQ1, xPQ1, ExponentFull-2)
            w2 = Weil_pairing_2power(Montgomery_coeff(a24d), xP2, xQ2, xPQ2, ExponentFull-2)
            @assert w1 * w2 == 1

            a24_1 = a24
            a24_2 = a24d
            e = -2
            n = norm(beta1)
            while n % 2 == 0
                n >>= 1
                e += 1
            end 
            println("e: ", e)
            e = max(e, 0)

            eval1 = [xP1, xQ1, xPQ1, xP, xQ, xPQ]
            eval2 = [xP2, xQ2, xPQ2, xPd, xQd, xPQd]
            # D + 2^(ExponentFull-2)*(P1, P2) for D in [(P, Pd), (Q, Qd)]
            for Mt = [[1 0; 0 0], [0 0; 1 0], [1 0; -1 0]]
                M1act = -s*Mt + (BigInt(1) << (ExponentFull - 2)) * M1
                M2act = -s*Mt + (BigInt(1) << (ExponentFull - 2)) * M2
                xP1T = action_of_matrix(M1act, a24, xP, xQ, xPQ, true)
                xP2T = action_of_matrix(M2act, a24d, xPd, xQd, xPQd, true)
                push!(eval1, xP1T)
                push!(eval2, xP2T)
            end

            # if e > 2 then we compute (2, 2)-isogenies by Velu's formulas
            if e > 2
                if is_infinity(xDBLe(xP1, a24, ExponentFull-3))
                    K = xDBLe(xQ1, a24, ExponentFull-e)
                else
                    K = xDBLe(xP1, a24, ExponentFull-e)
                end
                @assert is_infinity(xDBLe(K, a24, e-2))
                @assert !is_infinity(xDBLe(K, a24, e-3))
                a24_1, eval1 = two_e_iso(a24, K, e-2, eval1)

                if is_infinity(xDBLe(xP2, a24d, ExponentFull-3))
                    K = xDBLe(xQ2, a24d, ExponentFull-e)
                else
                    K = xDBLe(xP2, a24d, ExponentFull-e)
                end
                @assert is_infinity(xDBLe(K, a24d, e-2))
                @assert !is_infinity(xDBLe(K, a24d, e-3))
                a24_2, eval2 = two_e_iso(a24d, K, e-2, eval2)
            end

            # (2, 2)-isogenies by theta model
            xP1P2, xQ1Q2, xPQ1PQ2 = [CouplePoint(eval1[i], eval2[i]) for i in 1:3]
            eval_points = [CouplePoint(eval1[i], eval2[i]) for i in 4:6]
            eval_points_T = [CouplePoint(eval1[i], eval2[i]) for i in 7:9]

            exp = ExponentFull - e
            w1 = Weil_pairing_2power(Montgomery_coeff(a24_1), xP1P2[1], xP1P2[1], xPQ1PQ2[1], exp)
            w2 = Weil_pairing_2power(Montgomery_coeff(a24_2), xP1P2[2], xP1P2[2], xPQ1PQ2[2], exp)
            @assert w1 * w2 == 1
            if haskey(StrategiesDim2, exp)
                strategy = StrategiesDim2[exp]
            else
                strategy = compute_strategy(exp - 2, 2, 1)
            end
            Es, images = product_isogeny_sqrt(a24, a24d, xP1P2, xQ1Q2, xPQ1PQ2, eval_points, eval_points_T, exp, strategy)

            # determine the d-isogeny
            idx = 1
            xPim, xQim, xPQim = images[1][idx], images[2][idx], images[3][idx]
            A = Es[idx]
            w0 = Weil_pairing_2power(Montgomery_coeff(a24), xP, xQ, xPQ, ExponentFull)
            w1 = Weil_pairing_2power(affine(A), xPim, xQim, xPQim, ExponentFull)
            if w1 != w0^d
                idx = 2
            end
            xPim, xQim, xPQim = images[1][idx], images[2][idx], images[3][idx]
            A = Es[idx]
            w2 = Weil_pairing_2power(affine(A), xPim, xQim, xPQim, ExponentFull)
            @assert w2 == w0^d

            return A_to_a24(A), xPim, xQim, xPQim
        else
            xPd, xQd, xPQd = action_on_torsion_basis(d_inv * alpha, a24, xP, xQ, xPQ, global_data.E0_data)

            # compute an extra l-isogeny
            K = random_point_order_l(a24, p + 1, l)
            a24d, (xPd, xQd, xPQd) = odd_isogeny(a24, K, l, [xPd, xQd, xPQd])

            a24, xP, xQ, xPQ, _ = d2isogeny(a24, a24d, xP, xQ, xPQ, xPd, xQd, xPQd, ExponentFull, d, Proj1{FqFieldElem}[], global_data)
        end
    else
        a24, xP, xQ, xPQ, _ = d2isogeny_form_Esquare(a24, d, alpha, xP, xQ, xPQ, global_data)
    end

    return a24, xP, xQ, xPQ
end