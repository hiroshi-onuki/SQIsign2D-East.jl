# alpha(P), alpha(Q) for a fixed basis (P, Q) of E0[2^ExponentFull]
function action_on_torsion_basis(alpha::QOrderElem, E0_data::E0Data)
    xP, xQ, xPQ = E0_data.xP2e, E0_data.xQ2e, E0_data.xPQ2e
    Malpha = alpha[1] * [1 0; 0 1] + alpha[2] * E0_data.Matrices_2e[1] + alpha[3] * E0_data.Matrices_2e[2] + alpha[4] * E0_data.Matrices_2e[3]
    x_alphaP = linear_comb_2_e(Malpha[1, 1], Malpha[2, 1], xP, xQ, xPQ, E0_data.a24_0, ExponentFull)
    x_alphaQ = linear_comb_2_e(Malpha[1, 2], Malpha[2, 2], xP, xQ, xPQ, E0_data.a24_0, ExponentFull)
    x_alphaPQ = linear_comb_2_e(Malpha[1, 1] - Malpha[1, 2], Malpha[2, 1] - Malpha[2, 2], xP, xQ, xPQ, E0_data.a24_0, ExponentFull)
    return x_alphaP, x_alphaQ, x_alphaPQ
end

function RandIsogImages(d::BigInt, E0_data::E0Data)
    deg_dim2 = BigInt(1) << ExponentFull
    a24_0 = E0_data.a24_0

    alpha, _ = FullRepresentInteger(d*(deg_dim2 - d))

    n = 0
    if (alpha + deg_dim2) % 2 == 1
        # the first (2,2)-isogeny is E0^2 -> E0^2 represented by the matrix [1 1; 1 -1]
        xP1, xQ1, xPQ1 = action_on_torsion_basis(alpha + d, E0_data)
        xP2, xQ2, xPQ2 = action_on_torsion_basis(alpha - d, E0_data)
        n = norm(alpha + d)
    elseif (alpha + deg_dim2*Quaternion_i) % 2 == Quaternion_i
        # the first (2,2)-isogeny is E0^2 -> E0^2 represented by the matrix [i 1; 1 i]
        xP1, xQ1, xPQ1 = action_on_torsion_basis(Quaternion_i * alpha - d, E0_data)
        xP2, xQ2, xPQ2 = action_on_torsion_basis(alpha - d * Quaternion_i, E0_data)
        n = norm(Quaternion_i * alpha - d)
    else
        xP1, xQ1, xPQ1 = action_on_torsion_basis(alpha, E0_data)
        xP2 = ladder(deg_dim2 - d, E0_data.xP2e, a24_0)
        xQ2 = ladder(deg_dim2 - d, E0_data.xQ2e, a24_0)
        xPQ2 = ladder(deg_dim2 - d, E0_data.xPQ2e, a24_0)
    end

    e = 0 # e is the number of (2, 2)-isogenies computed without theta functions
    a24_1= a24_0
    a24_2 = a24_0
    if n != 0
        e = 1
        n >>= 2 # alpha - d or i*alpha - d is divisible by 2, so its norm is divisible by 4
        while n % 2 == 0
            n >>= 1
            e += 1
        end

        # if e > 1 then we compute (2, 2)-isogenies by Velu's formulas
        if e > 1
            if is_infinity(xDBLe(xP1, a24_0, ExponentFull-2))
                K = xDBLe(xQ1, a24_0, ExponentFull-e)
                Kd = xDBLe(xPQ1, a24_0, ExponentFull-e)
            else
                K = xDBLe(xP1, a24_0, ExponentFull-e)
                if is_infinity(xDBLe(xQ1, a24_0, ExponentFull-2))
                    Kd = xDBLe(xPQ1, a24_0, ExponentFull-e)
                else
                    Kd = xDBLe(xQ1, a24_0, ExponentFull-e)
                end
            end
            a24_d, _ = two_e_iso(a24_0, Kd, e-1, [xP1, xQ1, xPQ1])
            if is_infinity(xDBLe(xP2, a24_0, ExponentFull-2))
                K = xDBLe(xQ2, a24_0, ExponentFull-e)
            else
                K = xDBLe(xP2, a24_0, ExponentFull-e)
            end
            a24_2, (xP2, xQ2, xPQ2)= two_e_iso(a24_0, K, e-1, [xP2, xQ2, xPQ2])
        end
    end

    P1P2 = CouplePoint(xP1, xP2)
    Q1Q2 = CouplePoint(xQ1, xQ2)
    PQ1PQ2 = CouplePoint(xPQ1, xPQ2)

    strategy = compute_strategy(ExponentFull - e - 2, 2, 1)
    product_isogeny_sqrt(a24_1, a24_2, P1P2, Q1Q2, PQ1PQ2, CouplePoint{FqFieldElem}[], ExponentFull - e, strategy)
end