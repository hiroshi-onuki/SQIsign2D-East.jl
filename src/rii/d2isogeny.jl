# compute a (2^e, 2^e)-isogeny with kernel <(-dP0, alpha(P0)), -dQ0, alpha(Q0)> from E^2,
# where e = ExponentFull, d(2^e - d) = norm(alpha),
# (P0, Q0) is the images of the precomputed torion basis of E0[2^e] under an odd isogeny \phi,
# and the action of quaternion is induced by \phi.
function d2isogeny_form_Esquare(a24::Proj1{T}, d::BigInt, alpha::QOrderElem, xP0::Proj1{T}, xQ0::Proj1{T}, xPQ0::Proj1{T}, global_data::GlobalData, compute_odd_points::Bool=false) where T <: RingElem
    deg_dim2 = BigInt(1) << ExponentFull
    E0_data = global_data.E0_data
    if (alpha + deg_dim2) % 2 == 1
        # the first (2,2)-isogeny is E^2 -> E^2 represented by the matrix [1 -1; 1 1]
        beta1 = -d - alpha
        beta2 = -d + alpha
        gamma = Quaternion_1
        xR1, xS1, xRS1 = xP0, xQ0, xPQ0
        xR2, xS2, xRS2 = xP0, xQ0, xPQ0
    elseif a24 == E0_data.a24_0 && (alpha + deg_dim2*Quaternion_i) % 2 == Quaternion_i
        # this occurs only when a24 = a24_0
        # the first (2,2)-isogeny is E0^2 -> E0^2 represented by the matrix [1 i; i 1]
        beta1 = -d + Quaternion_i * alpha
        beta2 = -d * Quaternion_i + alpha
        gamma = Quaternion_i
        xR1, xS1, xRS1 = xP0, xQ0, xPQ0
        xR2, xS2, xRS2 = -xP0, -xQ0, -xPQ0
    else
        beta1 = -d * Quaternion_1
        beta2 = alpha
        gamma = Quaternion_0
        O0 = infinity_point(global_data.Fp2)
        xR1, xS1, xRS1 = xP0, xQ0, xPQ0
        xR2, xS2, xRS2 = O0, O0, O0
    end
    n = norm(beta1)
    xP1, xQ1, xPQ1 = action_on_torsion_basis(beta1, a24, xP0, xQ0, xPQ0, E0_data)
    xP2, xQ2, xPQ2 = action_on_torsion_basis(beta2, a24, xP0, xQ0, xPQ0, E0_data)

    a24_1 = a24
    a24_2 = a24
    e = -1
    while n % 2 == 0
        n >>= 1
        e += 1
    end
    e = max(e, 0)

    # compute R + T etc. for the gluing isogeny
    M1 = quaternion_to_matrix((BigInt(1) << (ExponentFull - e - 2)) * beta1, E0_data)
    xR1_T = action_of_matrix([1 0; 0 0] + M1, a24, xP0, xQ0, xPQ0, true)
    xS1_T = action_of_matrix([0 0; 1 0] + M1, a24, xP0, xQ0, xPQ0, true)
    xRS1_T = action_of_matrix([1 0; -1 0] + M1, a24, xP0, xQ0, xPQ0, true)
    M2 = quaternion_to_matrix((BigInt(1) << (ExponentFull - e - 2)) * beta2, E0_data)
    Mg = quaternion_to_matrix(gamma, E0_data)
    xR2_T = action_of_matrix(Mg*[1 0; 0 0] + M2, a24, xP0, xQ0, xPQ0, true)
    xS2_T = action_of_matrix(Mg*[0 0; 1 0] + M2, a24, xP0, xQ0, xPQ0, true)
    xRS2_T = action_of_matrix(Mg*[1 0; -1 0] + M2, a24, xP0, xQ0, xPQ0, true)

    # compute R + T for odd torsion points
    odd_x_points_1 = Proj1{FqFieldElem}[]
    odd_x_points_2 = Proj1{FqFieldElem}[]
    if compute_odd_points
        # this is only for E0^2
        P0, Q0 = E0_data.P2e, E0_data.Q2e
        T1 = add(mult(M1[1, 1], P0, Proj1(E0_data.A0)), mult(M1[2, 1], Q0, Proj1(E0_data.A0)), Proj1(E0_data.A0))
        T2 = add(mult(M2[1, 1], P0, Proj1(E0_data.A0)), mult(M2[2, 1], Q0, Proj1(E0_data.A0)), Proj1(E0_data.A0))
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
        if is_infinity(xDBLe(xP1, a24, ExponentFull-2))
            K = xDBLe(xQ1, a24, ExponentFull-e)
        else
            K = xDBLe(xP1, a24, ExponentFull-e)
        end
        eval_points = vcat([xP1, xQ1, xPQ1, xR1, xS1, xRS1, xR1_T, xS1_T, xRS1_T], odd_x_points_1)
        a24_1, images = two_e_iso(a24, K, e-1, eval_points)
        xP1, xQ1, xPQ1, xR1, xS1, xRS1, xR1_T, xS1_T, xRS1_T = images[1:9]
        odd_x_points_1 = images[10:end]
        if is_infinity(xDBLe(xP2, a24, ExponentFull-2))
            K = xDBLe(xQ2, a24, ExponentFull-e)
        else
            K = xDBLe(xP2, a24, ExponentFull-e)
        end
        eval_points = vcat([xP2, xQ2, xPQ2, xR2, xS2, xRS2, xR2_T, xS2_T, xRS2_T], odd_x_points_2)
        a24_2, images = two_e_iso(a24, K, e-1, eval_points)
        xP2, xQ2, xPQ2, xR2, xS2, xRS2, xR2_T, xS2_T, xRS2_T = images[1:9]
        odd_x_points_2 = images[10:end]
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

    idx = 1
    xP, xQ, xPQ = images[1][idx], images[2][idx], images[3][idx]
    A = Es[idx]
    if a24 == E0_data.a24_0
        w0 = E0_data.Weil_P2eQ2e
    else
        w0 = Weil_pairing_2power(Montgomery_coeff(a24), xP0, xQ0, xPQ0, ExponentFull)
    end
    w1 = Weil_pairing_2power(affine(A), xP, xQ, xPQ, ExponentFull)
    if w1 != w0^d
        idx = 2
    end
    xP, xQ, xPQ = images[1][idx], images[2][idx], images[3][idx]
    A = Es[idx]
    odd_images = [image[idx] for image in images[4:end]]

    # check the pairing
    w1 = Weil_pairing_2power(affine(A), xP, xQ, xPQ, ExponentFull)
    @assert w1 == w0^d
    # end

    return A_to_a24(A), xP, xQ, xPQ, odd_images
end

function d2isogeny()
end