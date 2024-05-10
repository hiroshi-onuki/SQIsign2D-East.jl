
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
function RandIsogImages(d::BigInt, global_data::GlobalData, compute_odd_points::Bool=false)
    deg_dim2 = BigInt(1) << ExponentFull
    E0_data = global_data.E0_data
    a24_0 = E0_data.a24_0
    xP0, xQ0, xPQ0 = E0_data.xP2e, E0_data.xQ2e, E0_data.xPQ2e

    alpha, _ = FullRepresentInteger(d*(deg_dim2 - d))
    println("n(alpha) % 3: ", norm(alpha) % 3)

    if (alpha + deg_dim2) % 2 == 1
        # the first (2,2)-isogeny is E0^2 -> E0^2 represented by the matrix [1 -1; 1 1]
        println("1 case")
        xP1, xQ1, xPQ1 = action_on_torsion_basis(-alpha - d, E0_data)
        xP2, xQ2, xPQ2 = action_on_torsion_basis(alpha - d, E0_data)
        xR1, xS1, xRS1 = xP0, xQ0, xPQ0
        xR2, xS2, xRS2 = xP0, xQ0, xPQ0
        n = norm(alpha + d)
        beta1 = -alpha - d
        beta2 = alpha - d
        gamma = Quaternion_1
    elseif (alpha + deg_dim2*Quaternion_i) % 2 == Quaternion_i
        # the first (2,2)-isogeny is E0^2 -> E0^2 represented by the matrix [1 i; i 1]
        println("i case")
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
    M1 = quaternion_to_matrix((BigInt(1) << (ExponentFull - e - 2)) * beta1, E0_data)
    xR1_T = action_of_matrix([1 0; 0 0] + M1, E0_data, true)
    xS1_T = action_of_matrix([0 0; 1 0] + M1, E0_data, true)
    xRS1_T = action_of_matrix([1 0; -1 0] + M1, E0_data, true)
    M2 = quaternion_to_matrix((BigInt(1) << (ExponentFull - e - 2)) * beta2, E0_data)
    Mg = quaternion_to_matrix(gamma, E0_data)
    xR2_T = action_of_matrix(Mg*[1 0; 0 0] + M2, E0_data, true)
    xS2_T = action_of_matrix(Mg*[0 0; 1 0] + M2, E0_data, true)
    xRS2_T = action_of_matrix(Mg*[1 0; -1 0] + M2, E0_data, true)

    # check by full point
    P0 = E0_data.P2e
    Q0 = E0_data.Q2e
    PQ0 = add(P0, -Q0, Proj1(E0_data.A0))
    T1 = add(mult(M1[1, 1], P0, Proj1(E0_data.A0)), mult(M1[2, 1], Q0, Proj1(E0_data.A0)), Proj1(E0_data.A0))
    T2 = add(mult(M2[1, 1], P0, Proj1(E0_data.A0)), mult(M2[2, 1], Q0, Proj1(E0_data.A0)), Proj1(E0_data.A0))
    R1_T = add(P0, T1, Proj1(E0_data.A0))
    S1_T = add(Q0, T1, Proj1(E0_data.A0))
    RS1_T = add(PQ0, T1, Proj1(E0_data.A0))
    if gamma == Quaternion_1
        R2_T = add(P0, T2, Proj1(E0_data.A0))
        S2_T = add(Q0, T2, Proj1(E0_data.A0))
        RS2_T = add(PQ0, T2, Proj1(E0_data.A0))
    elseif gamma == Quaternion_i
        iP0 = Point(-P0.X, global_data.Fp2_i * P0.Y, P0.Z)
        iQ0 = Point(-Q0.X, global_data.Fp2_i * Q0.Y, Q0.Z)
        iPQ0 = Point(-PQ0.X, global_data.Fp2_i * PQ0.Y, PQ0.Z)
        R2_T = add(iP0, T2, Proj1(E0_data.A0))
        S2_T = add(iQ0, T2, Proj1(E0_data.A0))
        RS2_T = add(iPQ0, T2, Proj1(E0_data.A0))
    else
        R2_T = T2
        S2_T = T2
        RS2_T = T2
    end
    @assert xR1_T == Proj1(R1_T.X, R1_T.Z)
    @assert xS1_T == Proj1(S1_T.X, S1_T.Z)
    @assert xRS1_T == Proj1(RS1_T.X, RS1_T.Z)
    @assert xR2_T == Proj1(R2_T.X, R2_T.Z)
    @assert xS2_T == Proj1(S2_T.X, S2_T.Z)
    @assert xRS2_T == Proj1(RS2_T.X, RS2_T.Z)
    # end

    odd_x_points_1 = Proj1{FqFieldElem}[]
    odd_x_points_2 = Proj1{FqFieldElem}[]
    if compute_odd_points
        T1 = Point(E0_data.A0, xDBLe(xP1, a24_0, ExponentFull - e - 2))
        T2 = Point(E0_data.A0, xDBLe(xP2, a24_0, ExponentFull - e - 2))
        for basis in E0_data.OddTorsionBases
            Pb1, Pb2 = basis
            if gamma == Quaternion_1
                Podd_1, Podd_2 = Pb1, Pb1
                Qodd_1, Qodd_2 = Pb2, Pb2
            elseif gamma == Quaternion_i
                iPb1 = Point(-Pb1.X, global_data.Fp2_i * Pb1.Y, Pb1.Z)
                iPb2 = Point(-Pb2.X, global_data.Fp2_i * Pb2.Y, Pb2.Z)
                Podd_1, Podd_2 = Pb1, iPb1
                Qodd_1, Qodd_2 = Pb2, iPb2  
            else
                Podd_1, Podd_2 = Pb1, infinity_full_point(global_data.Fp2)
                Qodd_1, Qodd_2 = Pb2, infinity_full_point(global_data.Fp2)
            end
            PQodd_1 = add(Podd_1, -Qodd_1, Proj1(E0_data.A0))
            PQodd_2 = add(Podd_2, -Qodd_2, Proj1(E0_data.A0))
            PTodd_1 = add(Podd_1, T1, Proj1(E0_data.A0))
            PTodd_2 = add(Podd_2, T2, Proj1(E0_data.A0))
            QTodd_1 = add(Qodd_1, T1, Proj1(E0_data.A0))
            QTodd_2 = add(Qodd_2, T2, Proj1(E0_data.A0))
            PQTodd_1 = add(PQodd_1, T1, Proj1(E0_data.A0))
            PQTodd_2 = add(PQodd_2, T2, Proj1(E0_data.A0))

            for P in [Podd_1, PTodd_1, Qodd_1, QTodd_1, PQodd_1, PQTodd_1]
                push!(odd_x_points_1, Proj1(P.X, P.Z))
            end
            for P in [Podd_2, PTodd_2, Qodd_2, QTodd_2, PQodd_2, PQTodd_2]
                if is_infinity(P)
                    push!(odd_x_points_2, infinity_point(global_data.Fp2))
                else
                    push!(odd_x_points_2, Proj1(P.X, P.Z))
                end
            end
        end
    end

    # if e > 1 then we compute (2, 2)-isogenies by Velu's formulas
    if e > 1
        if is_infinity(xDBLe(xP1, a24_0, ExponentFull-2))
            K = xDBLe(xQ1, a24_0, ExponentFull-e)
        else
            K = xDBLe(xP1, a24_0, ExponentFull-e)
        end
        eval_points = vcat([xP1, xQ1, xPQ1, xR1, xS1, xRS1, xR1_T, xS1_T, xRS1_T], odd_x_points_1)
        a24_1, images = two_e_iso(a24_0, K, e-1, eval_points)
        xP1, xQ1, xPQ1, xR1, xS1, xRS1, xR1_T, xS1_T, xRS1_T = images[1:9]
        odd_x_points_1 = images[10:end]
        if is_infinity(xDBLe(xP2, a24_0, ExponentFull-2))
            K = xDBLe(xQ2, a24_0, ExponentFull-e)
        else
            K = xDBLe(xP2, a24_0, ExponentFull-e)
        end
        eval_points = vcat([xP2, xQ2, xPQ2, xR2, xS2, xRS2, xR2_T, xS2_T, xRS2_T], odd_x_points_2)
        a24_2, images = two_e_iso(a24_0, K, e-1, eval_points)
        xP2, xQ2, xPQ2, xR2, xS2, xRS2, xR2_T, xS2_T, xRS2_T = images[1:9]
        odd_x_points_2 = images[10:end]
    end

    if compute_odd_points
        xT1 = xDBLe(xP1, a24_1, ExponentFull - e - 2)
        xT2 = xDBLe(xP2, a24_2, ExponentFull - e - 2)
        for i in 1:6
            if i % 2 == 1
                xP = odd_x_points_1[i]
                @assert is_infinity(ladder(27, xP, a24_1))
                @assert !is_infinity(ladder(9, xP, a24_1))
                xP = odd_x_points_2[i]
                if !is_infinity(xP)
                    @assert is_infinity(ladder(27, xP, a24_2))
                    @assert !is_infinity(ladder(9, xP, a24_2))
                end
            else
                tmp = x_add_sub(odd_x_points_1[i-1], xT1, a24_1)
                if tmp != odd_x_points_1[i]
                    @assert odd_x_points_1[i] == xADD(odd_x_points_1[i-1], xT1, tmp)
                end
                tmp = x_add_sub(odd_x_points_2[i-1], xT2, a24_2)
                if tmp != odd_x_points_2[i]
                    @assert odd_x_points_2[i] == xADD(odd_x_points_2[i-1], xT2, tmp)
                end
            end
        end
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
    odd_couple_points = CouplePoint{FqFieldElem}[]
    odd_couple_points_T = CouplePoint{FqFieldElem}[]
    for i in 1:div(length(odd_x_points_1), 2)
        n = 2*i - 1
        push!(odd_couple_points, CouplePoint(odd_x_points_1[n], odd_x_points_2[n]))
        push!(odd_couple_points_T, CouplePoint(odd_x_points_1[n+1], odd_x_points_2[n+1]))
    end

    if haskey(StrategiesDim2, ExponentFull - e)
        strategy = StrategiesDim2[ExponentFull - e]
    else
        strategy = compute_strategy(ExponentFull - e - 2, 2, 1)
    end
    eval_points = vcat([R1R2, S1S2, RS1RS2], odd_couple_points)
    eval_points_T = vcat([R1R2_T, S1S2_T, RS1RS2_T], odd_couple_points_T)
    Es, images = product_isogeny_sqrt(a24_1, a24_2, P1P2, Q1Q2, PQ1PQ2, eval_points, eval_points_T, ExponentFull - e, strategy)

    # check the pairing
    for idx in 1:2
        xP, xQ, xPQ = images[1][idx], images[2][idx], images[3][idx]
        A = Es[idx]
        w0 = E0_data.Weil_P2eQ2e
        w1 = Weil_pairing_2power(affine(A), xP, xQ, xPQ, ExponentFull)
        println(idx, ": pairing: ", w1 == w0^d)

        odd_images = [image[idx] for image in images[4:end]]
        for xPdd in odd_images
            println("order check: ", is_infinity(ladder(27, xPdd, A_to_a24(A))))
            println("             ", !is_infinity(ladder(9, xPdd, A_to_a24(A))))
        end
    end
    # end

    idx = 1
    xP, xQ, xPQ = images[1][idx], images[2][idx], images[3][idx]
    A = Es[idx]
    w0 = E0_data.Weil_P2eQ2e
    w1 = Weil_pairing_2power(affine(A), xP, xQ, xPQ, ExponentFull)
    if w1 != w0^d
        idx = 2
    end
    xP, xQ, xPQ = images[1][idx], images[2][idx], images[3][idx]
    A = Es[idx]
    odd_images = [image[idx] for image in images[4:end]]

    println("idx: ", idx, " e:", e)
    return A_to_a24(A), xP, xQ, xPQ, odd_images, LeftIdeal(alpha, d)
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