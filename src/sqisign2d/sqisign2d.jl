using SHA

# Sample a random ideal of prime norm 2^e for test
function sample_random_ideal_2e(e::Int)
    gamma = Quaternion_1
    while norm(gamma) % BigInt(2)^e != 0
        gamma, found = FullRepresentInteger(BigInt(2)^(Log2p + e))
        !found && continue
        gamma = div(gamma, gcd(gamma))
        if gcd(gamma * (Quaternion_1 - Quaternion_i)) % 2 == 0
            gamma = div(gamma * (Quaternion_1 - Quaternion_i), 2)
        end
    end
    I = LeftIdeal(gamma, BigInt(2)^e)
    a = rand(1:BigInt(2)^(e))
    return pushforward((1 + a) * Quaternion_1 + a * Quaternion_j, I)
end

# return a random prime <= 2^KLPT_secret_key_prime_size and = 3 mod 4
function random_secret_prime()
    B = BigInt(floor(p^(1/4)))
    n = rand(1:B)
    while !is_probable_prime(n)
        n = rand(1:B)
    end
    return n
end

function auxiliary_path(a24::Proj1{T}, xP::Proj1{T}, xQ::Proj1{T}, xPQ::Proj1{T}, odd_images::Vector{Proj1{T}},
                        I::LeftIdeal, nI::BigInt, q::BigInt, global_data::GlobalData) where T <: RingElem
    c = ExponentForTorsion
    r = (BigInt(1) << c) - q
    d = q * r
    a24d, xPd, xQd, xPQd = GeneralizedRandomIsogImages(d, a24, xP, xQ, xPQ, I, nI, global_data)

    q_inv = invmod(q, BigInt(1) << c)
    xP = xDBLe(xP, a24, ExponentFull - c)
    xQ = xDBLe(xQ, a24, ExponentFull - c)
    xPQ = xDBLe(xPQ, a24, ExponentFull - c)
    xPd = xDBLe(xPd, a24d, ExponentFull - c)
    xQd = xDBLe(xQd, a24d, ExponentFull - c)
    xPQd = xDBLe(xPQd, a24d, ExponentFull - c)
    xPd = ladder(q_inv, xPd, a24d)
    xQd = ladder(q_inv, xQd, a24d)
    xPQd = ladder(q_inv, xPQd, a24d)

    a24aux, xPaux, xQaux, xPQaux, images = d2isogeny(a24, a24d, xP, xQ, xPQ, xPd, xQd, xPQd, c, r, odd_images, global_data)

    return a24aux, xPaux, xQaux, xPQaux, images
end

function key_gen(global_data::GlobalData)
    D_sec = random_secret_prime()
    a24, xP, xQ, xPQ, odd_images, I_sec = RandIsogImages(D_sec, global_data, true)
    a24, images = Montgomery_normalize(a24, vcat([xP, xQ, xPQ], odd_images))
    xP, xQ, xPQ = images[1:3]
    odd_images = images[4:end]
    A = Montgomery_coeff(a24)

    nl = length(global_data.E0_data.DegreesOddTorsionBases)
    Ms = Vector{Matrix{Int}}(undef, nl)   
    for i in 1:nl
        l = global_data.E0_data.DegreesOddTorsionBases[i]
        e = global_data.E0_data.ExponentsOddTorsionBases[i]
        xPodd, xQodd, xPQodd = torsion_basis(a24, l, e)
        xPim, xQim, xPQim = odd_images[3*(i-1)+1:3*i]
        a, b, c, d = bi_dlog_odd_prime_power(A, xPim, xQim, xPQim, xPodd, xQodd, xPQodd, l, e)
        Ms[i] = [a c; b d]
    end

    return Montgomery_coeff(a24), (I_sec, D_sec, xP, xQ, xPQ, odd_images, Ms)
end

function commitment(global_data::GlobalData)
    D_sec = random_secret_prime()
    a24, xP, xQ, xPQ, I_sec = RandIsogImages(D_sec, global_data, false)
    a24, (xP, xQ, xPQ) = Montgomery_normalize(a24, [xP, xQ, xPQ])
    A = Montgomery_coeff(a24)
    xPc, xQc, xPQc = torsion_basis(a24, SQISIGN_challenge_length)
    xPd = xDBLe(xP, a24, ExponentFull - SQISIGN_challenge_length)
    xQd = xDBLe(xQ, a24, ExponentFull - SQISIGN_challenge_length)
    xPQd = xDBLe(xPQ, a24, ExponentFull - SQISIGN_challenge_length)
    n1, n2, n3, n4 = ec_bi_dlog_commitment(A, xPc, xQc, xPQc, xPd, xQd, xPQd, global_data.E0_data)
    M = [n1 n3; n2 n4]
    return A, (I_sec, D_sec, xP, xQ, xPQ, xPc, xQc, xPQc), M
end


function challenge(A::FqFieldElem, m::String, global_data::GlobalData)
    h = sha3_256(string(A) * m)

    c = BigInt(0)
    len = SQISIGN_challenge_length
    n, r = divrem(len, 8)
    for i in 1:(n+1)
        c += BigInt(h[i]) << (8*(i-1))
    end
    c >>= 8 - r

    return c
end

function signing(pk::FqFieldElem, sk, m::String, global_data::GlobalData)
    A = pk
    a24pub = A_to_a24(A)
    Isec, Dsec, xPsec, xQsec, xPQsec, odd_images, M_odd_images = sk

    # compute commitment
    Acom, (Icom, Dcom, xPcom, xQcom, xPQcom, xPcom_fix, xQcom_fix, xPQcom_fix), Mcom = commitment(global_data)
    a24com = A_to_a24(Acom)

    # compute challenge and the pull-back of the corresponding ideal
    cha = challenge(Acom, m, global_data)
    a, b = Mcom * [1, cha]
    a, b, c, d = global_data.E0_data.Matrix_2ed_inv * [b, 0, -a, 0]
    alpha = SQIsign2D.Level1.QOrderElem(a, b, c, d)
    Icha = SQIsign2D.Level1.LeftIdeal(alpha, BigInt(1) << SQISIGN_challenge_length)
    Kcha = ladder3pt(cha, xPcom_fix, xQcom_fix, xPQcom_fix, a24com)
    a24cha, (xPcha, xQcha, xPQcha) = two_e_iso(a24com, Kcha, SQISIGN_challenge_length, [xPcom, xQcom, xPQcom], StrategyChallenge)
    a24cha, (xPcha, xQcha, xPQcha) = Montgomery_normalize(a24cha, [xPcha, xQcha, xPQcha])
    Acha = Montgomery_coeff(a24cha)

    # find alpha in bar(Isec)IcomIcha suitable for the response
    Icomcha = SQIsign2D.Level1.intersection(Icom, Icha)
    I = SQIsign2D.Level1.involution_product(Isec, Icomcha)
    nI = Dsec * Dcom << SQISIGN_challenge_length
    alpha, d, found = element_for_response(I, nI, ExponentForTorsion, [(3, 3)], Dsec)

    if found

        # compute the image under the response sigma
        Malpha = quaternion_to_matrix(involution(alpha), global_data.E0_data.Matrices_2e)
        xPres, xQres, xPQres = action_of_matrix(Malpha, a24cha, xPcha, xQcha, xPQcha)
        Dcom_inv = invmod(Dcom, BigInt(1) << ExponentForTorsion)
        xPres, xQres, xPQres = ladder(Dcom_inv, xPres, a24cha), ladder(Dcom_inv, xQres, a24cha), ladder(Dcom_inv, xPQres, a24cha)
        xPfix, xQfix, xPQfix = torsion_basis(a24cha, ExponentForTorsion)
        n1, n2, n3, n4 = ec_bi_dlog_response(Acha, xPfix, xQfix, xPQfix, xPres, xQres, xPQres, global_data.E0_data)
        Mres = [n1 n3; n2 n4]

        q = div(norm(alpha), d*nI)
        n_odd_l = length(global_data.E0_data.DegreesOddTorsionBases)
        odd_kernels = Proj1{FqFieldElem}[]
        odd_kernel_coeffs = Tuple{Int, Int}[]
        for i in 1:n_odd_l
            l = global_data.E0_data.DegreesOddTorsionBases[i]
            e = global_data.E0_data.ExponentsOddTorsionBases[i]
            if d % l == 0
                f = Int(log(l, gcd(d, l^e)))
                if f > 0
                    xPodd, xQodd, xPQodd = odd_images[3*(i-1)+1:3*i]
                    xPodd = ladder(l^(e-f), xPodd, a24pub)
                    xQodd = ladder(l^(e-f), xQodd, a24pub)
                    xPQodd = ladder(l^(e-f), xPQodd, a24pub)
                    a, b = kernel_coefficients(involution(alpha), l, f, global_data.E0_data.Matrices_odd[i])
                    a, b = l^(e-f) * M_odd_images[i] * [a, b]
                    Kodd = kernel_generator(xPodd, xQodd, xPQodd, a24pub, involution(alpha), l, f, global_data.E0_data.Matrices_odd[i])
                else
                    a, b = 0, 0
                    Kodd = infinity_point(global_data.Fp2)
                end
                push!(odd_kernels, Kodd)
                push!(odd_kernel_coeffs, (a, b))
            end
        end

        # compute the auxiliary ellitic curve
        a24aux, xPaux, xQaux, xPQaux, images = auxiliary_path(a24pub, xPsec, xQsec, xPQsec, odd_kernels, Isec, Dsec, q, global_data)

        dd = d
        for i in 1:n_odd_l
            l = global_data.E0_data.DegreesOddTorsionBases[i]
            e = 0
            while dd % l == 0
                dd = div(dd, l)
                e += 1
            end
            # determine the order of image[i] 
            K = images[i]
            f = 0
            while !is_infinity(K)
                K = ladder(l, K, a24aux)
                f += 1
            end
            for k in 1:f
                K = ladder(l^(f - k), images[i], a24aux)
                a24aux, tmp = odd_isogeny(a24aux, K, l, vcat([xPaux, xQaux, xPQaux], images))
                xPaux, xQaux, xPQaux = tmp[1:3]
                images = tmp[4:end]
            end
            
            # compute the rest l^(e - f)-isogeny
            n = div(e - f, 2)
            if n > 0
                xPaux, xQaux, xPQaux = ladder(l^n, xPaux, a24aux), ladder(l^n, xQaux, a24aux), ladder(l^n, xPQaux, a24aux)
            end
            if (e - f) % 2 == 1
                K = random_point_order_l(a24aux, p + 1, l)
                a24aux, tmp = odd_isogeny(a24aux, K, l, vcat([xPaux, xQaux, xPQaux], images))
                xPaux, xQaux, xPQaux = tmp[1:3]
                images = tmp[4:end]
            end
        end

        Aaux = Montgomery_coeff(a24aux)
        xPaux, xQaux, xPQaux = action_of_matrix(Mres, a24aux, xPaux, xQaux, xPQaux)

        return (Acom, Aaux, xPaux, xQaux, xPQaux, odd_kernel_coeffs), true
    end
    return nothing, false
end

function verify(pk::FqFieldElem, sign, m::String, global_data::GlobalData)
    Acom, Aaux, xPaux, xQaux, xPQaux, odd_kernel_coeffs = sign
    a24pub = A_to_a24(pk)
    a24com = A_to_a24(Acom)
    a24aux = A_to_a24(Aaux)
    
    c = challenge(Acom, m, global_data)
    xPcom, xQcom, xPQcom = torsion_basis(a24com, SQISIGN_challenge_length)
    Kcha = ladder3pt(c, xPcom, xQcom, xPQcom, a24com)
    a24cha, _ = two_e_iso(a24com, Kcha, SQISIGN_challenge_length, Proj1{FqFieldElem}[], StrategyChallenge)
    a24cha, _ = Montgomery_normalize(a24cha, Proj1{FqFieldElem}[])
    xPres, xQres, xPQres = torsion_basis(a24cha, ExponentForTorsion)

    # isogeny of dimension 3
    P1P2 = CouplePoint(xPres, xPaux)
    Q1Q2 = CouplePoint(xQres, xQaux)
    PQ1PQ2 = CouplePoint(xPQres, xPQaux)
    Es, _ = product_isogeny_sqrt(a24cha, a24aux, P1P2, Q1Q2, PQ1PQ2, CouplePoint{FqFieldElem}[], CouplePoint{FqFieldElem}[], ExponentForTorsion, StrategiesDim2[ExponentForTorsion])

    a24mid = a24pub
    n_odd_l = length(global_data.E0_data.DegreesOddTorsionBases)
    for i in 1:n_odd_l
        l = global_data.E0_data.DegreesOddTorsionBases[i]
        e = global_data.E0_data.ExponentsOddTorsionBases[i]
        a, b = odd_kernel_coeffs[i]
        g = gcd(a, b, l^e)
        d = div(l^e, g)
        println(g)
        if d > 0
            xPodd, xQodd, xPQodd = torsion_basis(a24mid, l, e)
            xPodd = ladder(g, xPodd, a24mid)
            xQodd = ladder(g, xQodd, a24mid)
            xPQodd = ladder(g, xPQodd, a24mid)
            a, b = div(a, g), div(b, g)
            if a % l == 0
                a = (a * invmod(b, d)) % d
                Kfull = ladder3pt(a, xQodd, xPodd, xPQodd, a24mid)
            else
                b = (b * invmod(a, d)) % d
                Kfull = ladder3pt(b, xPodd, xQodd, xPQodd, a24mid)
            end
            e = Int(log(l, d))
            for k in 1:e
                K = ladder(l^(e - k), Kfull, a24mid)
                a24mid, tmp = odd_isogeny(a24mid, K, l, [Kfull])
                Kfull = tmp[1]
            end
        end
    end

    j0 = jInvariant_a24(a24mid)
    j1 = jInvariant_A(Es[1])
    j2 = jInvariant_A(Es[2])

    return j1 == j0 || j2 == j0
end