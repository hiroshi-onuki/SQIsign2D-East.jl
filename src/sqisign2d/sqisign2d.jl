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
        l, e = global_data.E0_data.DegreesOddTorsionBases[i]
        xPodd, xQodd, xPQodd = torsion_basis(A, l, e)
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
    xPc, xQc, xPQc = torsion_basis(A, SQISIGN_challenge_length)
    xPd = xDBLe(xP, a24, ExponentFull - SQISIGN_challenge_length)
    xQd = xDBLe(xQ, a24, ExponentFull - SQISIGN_challenge_length)
    xPQd = xDBLe(xPQ, a24, ExponentFull - SQISIGN_challenge_length)
    n1, n2, n3, n4 = ec_bi_dlog_commitment(A, xPc, xQc, xPQc, xPd, xQd, xPQd, global_data.E0_data)
    M = [n1 n3; n2 n4]
    return A, (I_sec, D_sec, xP, xQ, xPQ, xPc, xQc, xPQc), M
end


function challenge(A::FqFieldElem, m::String)
    if SQISIGN_challenge_length <= 256
        h = sha3_256(string(A) * m)
    else
        h = sha3_512(string(A) * m)
    end

    c = BigInt(0)
    len = SQISIGN_challenge_length
    n, r = divrem(len, 8)
    for i in 1:(n+1)
        c += BigInt(h[i]) << (8*(i-1))
    end
    c >>= 8 - r

    return c
end

function signing(pk::FqFieldElem, sk, m::String, global_data::GlobalData, is_compact::Bool=false)
    A = pk
    a24pub = A_to_a24(A)
    Isec, Dsec, xPsec, xQsec, xPQsec, odd_images, M_odd_images = sk

    while true
        # compute commitment
        Acom, (Icom, Dcom, xPcom, xQcom, xPQcom, xPcom_fix, xQcom_fix, xPQcom_fix), Mcom = commitment(global_data)
        a24com = A_to_a24(Acom)

        # compute challenge and the pull-back of the corresponding ideal
        cha = challenge(Acom, m)
        a, b = Mcom * [1, cha]
        a, b, c, d = global_data.E0_data.Matrix_2ed_inv * [b, 0, -a, 0]
        alpha = QOrderElem(a, b, c, d)
        Icha = LeftIdeal(alpha, BigInt(1) << SQISIGN_challenge_length)
        Kcha = ladder3pt(cha, xPcom_fix, xQcom_fix, xPQcom_fix, a24com)
        if !is_compact
            a24cha, (xPcha, xQcha, xPQcha) = two_e_iso(a24com, Kcha, SQISIGN_challenge_length, [xPcom, xQcom, xPQcom], StrategyChallenge)
            a24cha, (xPcha, xQcha, xPQcha) = Montgomery_normalize(a24cha, [xPcha, xQcha, xPQcha])
        else
            a24cha, (xPcha, xQcha, xPQcha, Kcha_dual) = two_e_iso(a24com, Kcha, SQISIGN_challenge_length, [xPcom, xQcom, xPQcom, xQcom_fix], StrategyChallenge)
            a24cha, (xPcha, xQcha, xPQcha, Kcha_dual) = Montgomery_normalize(a24cha, [xPcha, xQcha, xPQcha, Kcha_dual])
        end
        Acha = Montgomery_coeff(a24cha)

        # find alpha in bar(Isec)IcomIcha suitable for the response
        Icomcha = intersection(Icom, Icha)
        I = involution_product(Isec, Icomcha)
        nI = Dsec * Dcom << SQISIGN_challenge_length
        alpha, d, found = element_for_response(I, nI, ExponentForTorsion, global_data.E0_data.DegreesOddTorsionBases, Dsec)
        !found && continue

        # compute the image under the response sigma
        Malpha = quaternion_to_matrix(involution(alpha), global_data.E0_data.Matrices_2e)
        xPres, xQres, xPQres = action_of_matrix(Malpha, a24cha, xPcha, xQcha, xPQcha)
        Dcom_inv = invmod(Dcom, BigInt(1) << ExponentForTorsion)
        xPres, xQres, xPQres = ladder(Dcom_inv, xPres, a24cha), ladder(Dcom_inv, xQres, a24cha), ladder(Dcom_inv, xPQres, a24cha)

        q = div(norm(alpha), d*nI)
        n_odd_l = length(global_data.E0_data.DegreesOddTorsionBases)
        odd_kernels = Proj1{FqFieldElem}[]
        odd_kernel_coeffs = Tuple{Int, Int}[]
        for i in 1:n_odd_l
            l, e = global_data.E0_data.DegreesOddTorsionBases[i]
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

        # compute the auxiliary ellitic curve
        a24aux, xPaux, xQaux, xPQaux, images = auxiliary_path(a24pub, xPsec, xQsec, xPQsec, odd_kernels, Isec, Dsec, q, global_data)

        dd = d
        for i in 1:n_odd_l
            l, _ = global_data.E0_data.DegreesOddTorsionBases[i]
            e = 0
            while dd % l == 0
                dd = div(dd, l)
                e += 1
            end
            e == 0 && continue
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

        if !is_compact
            # modify xPaux, xQaux, xPQaux to be the images of the fixed torions
            xPfix, xQfix, xPQfix = torsion_basis(Acha, ExponentForTorsion)
            n1, n2, n3, n4 = ec_bi_dlog_response(Acha, xPfix, xQfix, xPQfix, xPres, xQres, xPQres, global_data.E0_data)
            Mres = [n1 n3; n2 n4]
            xPaux, xQaux, xPQaux = action_of_matrix(Mres, a24aux, xPaux, xQaux, xPQaux)
    
            # compress the signature
            sign = Vector{UInt8}(undef, SQISIGN2D_signature_length)
            idx = 1
            Acom_byte = Fq_to_bytes(Acom)
            sign[idx:idx+SQISIGN2D_Fp2_length-1] = Acom_byte
            idx += SQISIGN2D_Fp2_length
            Aaux_byte = Fq_to_bytes(Aaux)
            sign[idx:idx+SQISIGN2D_Fp2_length-1] = Aaux_byte
            idx += SQISIGN2D_Fp2_length

            xPfix, xQfix, xPQfix = torsion_basis(Aaux, ExponentForTorsion)
            n1, n2, n3, n4 = ec_bi_dlog_response(Aaux, xPaux, xQaux, xPQaux, xPfix, xQfix, xPQfix, global_data.E0_data)
            for n in [n1, n2, n3, n4]
                n_byte = integer_to_bytes(n, SQISIGN2D_2a_length)
                sign[idx:idx+SQISIGN2D_2a_length-1] = n_byte
                idx += SQISIGN2D_2a_length
            end

        else
            xPsec, xQsec, xPQsec = xDBLe(xPsec, a24pub, ExponentFull - ExponentForTorsion), xDBLe(xQsec, a24pub, ExponentFull - ExponentForTorsion), xDBLe(xPQsec, a24pub, ExponentFull - ExponentForTorsion)
            eval_points = vcat([xPsec, xQsec, xPQsec], odd_kernels)
            a24mid = a24pub
            for i in 1:n_odd_l
                l, e = global_data.E0_data.DegreesOddTorsionBases[i]
                a, b = odd_kernel_coeffs[i]
                g = gcd(a, b, l^e)
                f = e - Int(log(l, g))
                for k in 1:f
                    K = ladder(l^(f - k), eval_points[i + 3], a24mid)
                    @assert is_infinity(ladder(l, K, a24mid))
                    @assert !is_infinity(K)
                    a24mid, eval_points = odd_isogeny(a24mid, K, l, eval_points)
                end
            end
            xPmid, xQmid, xPQmid = eval_points[1:3]
            @assert is_infinity(xDBLe(xPmid, a24mid, ExponentForTorsion))
            @assert is_infinity(xDBLe(xQmid, a24mid, ExponentForTorsion))
            @assert is_infinity(xDBLe(xPQmid, a24mid, ExponentForTorsion))
            @assert !is_infinity(xDBLe(xPmid, a24mid, ExponentForTorsion-1))
            @assert !is_infinity(xDBLe(xQmid, a24mid, ExponentForTorsion-1))
            @assert !is_infinity(xDBLe(xPQmid, a24mid, ExponentForTorsion-1))
            a24mid, (xPmid, xQmid, xPQmid) = Montgomery_normalize(a24mid, [xPmid, xQmid, xPQmid])
            Amid = Montgomery_coeff(a24mid)
            xPfix, xQfix, xPQfix = torsion_basis(Amid, ExponentForTorsion)
            n1, n2, n3, n4 = ec_bi_dlog_response(Amid, xPfix, xQfix, xPQfix, xPmid, xQmid, xPQmid, global_data.E0_data)

            a24aux, xPaux, xQaux, xPQaux, _ = d2isogeny(a24aux, a24cha, xPaux, xQaux, xPQaux, xPres, xQres, xPQres, ExponentForTorsion, q, Proj1{FqFieldElem}[], global_data)
            Aaux = Montgomery_coeff(a24aux)
            xPfix, xQfix, xPQfix = torsion_basis(Aaux, ExponentForTorsion)
            m1, m2, m3, m4 = ec_bi_dlog_response(Aaux, xPaux, xQaux, xPQaux, xPfix, xQfix, xPQfix, global_data.E0_data)
            c = invmod(BigInt(1) << ExponentForTorsion - q, BigInt(1) << ExponentForTorsion)
            M = c * [m1 m3; m2 m4] * [n1 n3; n2 n4]

            # test
            xP1, xQ1, xPQ1 = torsion_basis(Amid, ExponentForTorsion)
            xP2, xQ2, xPQ2 = action_of_matrix(M, a24aux, xPfix, xQfix, xPQfix)
            a24cha_d, _, _, _, _ = d2isogeny(a24mid, a24aux, xP1, xQ1, xPQ1, xP2, xQ2, xPQ2, ExponentForTorsion, q, Proj1{FqFieldElem}[], global_data)
            @assert jInvariant_a24(a24cha_d) == jInvariant_a24(a24cha)

            # compress the signature
            sign = Vector{UInt8}(undef, CompactSQISIGN2D_signature_length)
            idx = 1
            Aaux_byte = Fq_to_bytes(Aaux)
            sign[idx:idx+SQISIGN2D_Fp2_length-1] = Aaux_byte
            idx += SQISIGN2D_Fp2_length

            for n in [M[1, 1], M[2, 1], M[1, 2], M[2, 2]]
                n_byte = integer_to_bytes(n, SQISIGN2D_2a_length)
                sign[idx:idx+SQISIGN2D_2a_length-1] = n_byte
                idx += SQISIGN2D_2a_length
            end

            xPcha, xQcha, xPQcha = torsion_basis(Acha, SQISIGN_challenge_length)
            a, b = ec_bi_dlog_challenge(Acha, Kcha_dual, xPcha, xQcha, xPQcha, global_data.E0_data)
            if a % 2 == 1
                b = (b * invmod(a, BigInt(1) << SQISIGN_challenge_length)) % (BigInt(1) << SQISIGN_challenge_length)
                sign[idx] = 1
                sign[idx+1:idx+SQISIGN2D_2a_length] = integer_to_bytes(b, SQISIGN2D_2a_length)
                P = xQcha
            else
                a = (a * invmod(b, BigInt(1) << SQISIGN_challenge_length)) % (BigInt(1) << SQISIGN_challenge_length)
                sign[idx] = 0
                sign[idx+1:idx+SQISIGN2D_2a_length] = integer_to_bytes(a, SQISIGN2D_2a_length)
                P = xPcha
            end
            idx += SQISIGN2D_2a_length + 1
            a24com_d, tmp = two_e_iso(a24cha, Kcha_dual, SQISIGN_challenge_length, [P], StrategyChallenge)
            a24com_d, tmp = Montgomery_normalize(a24com_d, [tmp[1]])
            Kcha_d = tmp[1]
            @assert a24com_d == a24com
            r = ec_dlog(Acom, Kcha, Kcha_d, xQcom_fix, global_data.E0_data)
            @assert Kcha == ladder(r, Kcha_d, a24com)
            sign[idx:idx+SQISIGN2D_2a_length-1] = integer_to_bytes(r, SQISIGN2D_2a_length)
            idx += SQISIGN2D_2a_length

        end

        # coefficient (a:b) is of the form (l^f:b), where 0 < f <= e
        for i in 1:n_odd_l
            l, e = global_data.E0_data.DegreesOddTorsionBases[i]
            a, b = odd_kernel_coeffs[i]
            ad = gcd(a, l^e)
            if ad == l^e
                ea = e
                b = b % l^e
            else
                inv = invmod(div(a, ad), l^e)
                a = (a * inv) % l^e
                b = (b * inv) % l^e
                ea = Int(log(l, a))
            end
            ab_byte = integer_to_bytes(ea * l^e + b, 1)
            sign[idx] = ab_byte[1]
            idx += 1
        end

        return sign
    end
end

function verify(pk::FqFieldElem, sign::Vector{UInt8}, m::String, global_data::GlobalData)
    # decompress the signature
    idx = 1
    Acom = bytes_to_Fq(sign[idx:idx+SQISIGN2D_Fp2_length-1], global_data.Fp2)
    a24com = A_to_a24(Acom)
    idx += SQISIGN2D_Fp2_length
    Aaux = bytes_to_Fq(sign[idx:idx+SQISIGN2D_Fp2_length-1], global_data.Fp2)
    a24aux = A_to_a24(Aaux)
    idx += SQISIGN2D_Fp2_length
    n = Vector{BigInt}(undef, 4)
    for i in 1:4
        n[i] = bytes_to_integer(sign[idx:idx+SQISIGN2D_2a_length-1])
        idx += SQISIGN2D_2a_length
    end
    xPfix, xQfix, xPQfix = torsion_basis(Aaux, ExponentForTorsion)
    xPaux = linear_comb_2_e(n[1], n[2], xPfix, xQfix, xPQfix, a24aux, ExponentForTorsion)
    xQaux = linear_comb_2_e(n[3], n[4], xPfix, xQfix, xPQfix, a24aux, ExponentForTorsion)
    xPQaux = linear_comb_2_e(n[1]- n[3], n[2] - n[4], xPfix, xQfix, xPQfix, a24aux, ExponentForTorsion)
    n_odd_l = length(global_data.E0_data.DegreesOddTorsionBases)
    odd_kernel_coeffs = Vector{Tuple{Int, Int}}(undef, n_odd_l)
    for i in 1:n_odd_l
        l, e = global_data.E0_data.DegreesOddTorsionBases[i]
        ab = sign[idx]
        idx += 1
        ea = div(ab, l^e)
        a = l^ea % l^e
        b = ab % l^e
        odd_kernel_coeffs[i] = (a, b)
    end

    a24com = A_to_a24(Acom)
    a24aux = A_to_a24(Aaux)
    a24pub = A_to_a24(pk)

    c = challenge(Acom, m)
    xPcom, xQcom, xPQcom = torsion_basis(Acom, SQISIGN_challenge_length)
    Kcha = ladder3pt(c, xPcom, xQcom, xPQcom, a24com)
    a24cha, _ = two_e_iso(a24com, Kcha, SQISIGN_challenge_length, Proj1{FqFieldElem}[], StrategyChallenge)
    a24cha, _ = Montgomery_normalize(a24cha, Proj1{FqFieldElem}[])
    Acha = Montgomery_coeff(a24cha)
    xPres, xQres, xPQres = torsion_basis(Acha, ExponentForTorsion)

    # isogeny of dimension 2
    P1P2 = CouplePoint(xPres, xPaux)
    Q1Q2 = CouplePoint(xQres, xQaux)
    PQ1PQ2 = CouplePoint(xPQres, xPQaux)
    Es, _ = product_isogeny_sqrt(a24cha, a24aux, P1P2, Q1Q2, PQ1PQ2, CouplePoint{FqFieldElem}[], CouplePoint{FqFieldElem}[], ExponentForTorsion, StrategiesDim2[ExponentForTorsion])

    n_odd_l = length(global_data.E0_data.DegreesOddTorsionBases)
    odd_isog_kers = Proj1{FqFieldElem}[]
    odd_isog_degrees = Tuple{Int, Int}[]
    for i in 1:n_odd_l
        l, e = global_data.E0_data.DegreesOddTorsionBases[i]
        a, b = odd_kernel_coeffs[i]
        g = gcd(a, b, l^e)
        d = div(l^e, g)
        if d > 0
            xPodd, xQodd, xPQodd = torsion_basis(pk, l, e)
            xPodd = ladder(g, xPodd, a24pub)
            xQodd = ladder(g, xQodd, a24pub)
            xPQodd = ladder(g, xPQodd, a24pub)
            a, b = div(a, g), div(b, g)
            if a % l == 0
                a = (a * invmod(b, d)) % d
                Kfull = ladder3pt(a, xQodd, xPodd, xPQodd, a24pub)
            else
                b = (b * invmod(a, d)) % d
                Kfull = ladder3pt(b, xPodd, xQodd, xPQodd, a24pub)
            end
            e = Int(log(l, d))
            push!(odd_isog_kers, Kfull)
            push!(odd_isog_degrees, (l, e))
        end
    end

    a24mid = a24pub
    n_isog = length(odd_isog_kers)
    for i in 1:n_isog
        Kfull = odd_isog_kers[i]
        l, e = odd_isog_degrees[i]
        for k in 1:e
            K = ladder(l^(e - k), Kfull, a24mid)
            a24mid, odd_isog_kers = odd_isogeny(a24mid, K, l, odd_isog_kers)
            Kfull = odd_isog_kers[1]
        end
    end

    j0 = jInvariant_a24(a24mid)
    j1 = jInvariant_A(Es[1])
    j2 = jInvariant_A(Es[2])

    return j1 == j0 || j2 == j0
end

function verify_compact(pk::FqFieldElem, sign::Vector{UInt8}, m::String, global_data::GlobalData)
    # decompress the signature
    idx = 1
    Aaux = bytes_to_Fq(sign[idx:idx+SQISIGN2D_Fp2_length-1], global_data.Fp2)
    idx += SQISIGN2D_Fp2_length
    n = Vector{BigInt}(undef, 4)
    for i in 1:4
        n[i] = bytes_to_integer(sign[idx:idx+SQISIGN2D_2a_length-1])
        idx += SQISIGN2D_2a_length
    end
    xPfix, xQfix, xPQfix = torsion_basis(Aaux, ExponentForTorsion)
    a24aux = A_to_a24(Aaux)
    xPaux = linear_comb_2_e(n[1], n[2], xPfix, xQfix, xPQfix, a24aux, ExponentForTorsion)
    xQaux = linear_comb_2_e(n[3], n[4], xPfix, xQfix, xPQfix, a24aux, ExponentForTorsion)
    xPQaux = linear_comb_2_e(n[1]- n[3], n[2] - n[4], xPfix, xQfix, xPQfix, a24aux, ExponentForTorsion)

    bit_s = sign[idx]
    idx += 1
    s = bytes_to_integer(sign[idx:idx+SQISIGN2D_2a_length-1])
    idx += SQISIGN2D_2a_length
    r = bytes_to_integer(sign[idx:idx+SQISIGN2D_2a_length-1])
    idx += SQISIGN2D_2a_length

    n_odd_l = length(global_data.E0_data.DegreesOddTorsionBases)
    odd_kernel_coeffs = Vector{Tuple{Int, Int}}(undef, n_odd_l)
    for i in 1:n_odd_l
        l, e = global_data.E0_data.DegreesOddTorsionBases[i]
        ab = sign[idx]
        idx += 1
        ea = div(ab, l^e)
        a = l^ea % l^e
        b = ab % l^e
        odd_kernel_coeffs[i] = (a, b)
    end

    a24pub = A_to_a24(pk)

    # isogeny of dimension 1
    n_odd_l = length(global_data.E0_data.DegreesOddTorsionBases)
    odd_isog_kers = Proj1{FqFieldElem}[]
    odd_isog_degrees = Tuple{Int, Int}[]
    for i in 1:n_odd_l
        l, e = global_data.E0_data.DegreesOddTorsionBases[i]
        a, b = odd_kernel_coeffs[i]
        g = gcd(a, b, l^e)
        d = div(l^e, g)
        if d > 0
            xPodd, xQodd, xPQodd = torsion_basis(pk, l, e)
            xPodd = ladder(g, xPodd, a24pub)
            xQodd = ladder(g, xQodd, a24pub)
            xPQodd = ladder(g, xPQodd, a24pub)
            a, b = div(a, g), div(b, g)
            if a % l == 0
                a = (a * invmod(b, d)) % d
                Kfull = ladder3pt(a, xQodd, xPodd, xPQodd, a24pub)
            else
                b = (b * invmod(a, d)) % d
                Kfull = ladder3pt(b, xPodd, xQodd, xPQodd, a24pub)
            end
            e = Int(log(l, d))
            push!(odd_isog_kers, Kfull)
            push!(odd_isog_degrees, (l, e))
        end
    end
    a24mid = a24pub
    n_isog = length(odd_isog_kers)
    for i in 1:n_isog
        Kfull = odd_isog_kers[i]
        l, e = odd_isog_degrees[i]
        for k in 1:e
            K = ladder(l^(e - k), Kfull, a24mid)
            a24mid, odd_isog_kers = odd_isogeny(a24mid, K, l, odd_isog_kers)
            Kfull = odd_isog_kers[1]
        end
    end

    a24mid, _ = Montgomery_normalize(a24mid, Proj1{FqFieldElem}[])
    Amid = Montgomery_coeff(a24mid)
    xPmid, xQmid, xPQmid = torsion_basis(Amid, ExponentForTorsion)

    # isogeny of dimension 2
    P1P2 = CouplePoint(xPmid, xPaux)
    Q1Q2 = CouplePoint(xQmid, xQaux)
    PQ1PQ2 = CouplePoint(xPQmid, xPQaux)
    Es, _ = product_isogeny_sqrt(a24mid, a24aux, P1P2, Q1Q2, PQ1PQ2, CouplePoint{FqFieldElem}[], CouplePoint{FqFieldElem}[], ExponentForTorsion, StrategiesDim2[ExponentForTorsion])

    for E in Es
        a24cha = A_to_a24(E)
        a24cha, _ = Montgomery_normalize(a24cha, Proj1{FqFieldElem}[])
        Acha = Montgomery_coeff(a24cha)
        xPcha, xQcha, xPQcha = torsion_basis(Acha, SQISIGN_challenge_length)
        if bit_s == 1
            Kcha_dual = ladder3pt(s, xPcha, xQcha, xPQcha, a24cha)
            P = xQcha
        else
            Kcha_dual = ladder3pt(s, xQcha, xPcha, xPQcha, a24cha)
            P = xPcha
        end
        a24com, tmp = two_e_iso(a24cha, Kcha_dual, SQISIGN_challenge_length, [P], StrategyChallenge)
        a24com, tmp = Montgomery_normalize(a24com, [tmp[1]])
        Kcha_d = tmp[1]
        Acom = Montgomery_coeff(a24com)
        c = challenge(Acom, m)
        xPcom, xQcom, xPQcom = torsion_basis(Acom, SQISIGN_challenge_length)
        Kcha = ladder3pt(c, xPcom, xQcom, xPQcom, a24com)
        Kcha == ladder(r, Kcha_d, a24com) && return true
    end

    return false
end