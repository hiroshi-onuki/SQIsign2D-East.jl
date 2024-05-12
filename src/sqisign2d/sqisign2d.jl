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
                        I::LeftIdeal, nI::BigInt, q::BigInt, c::Int, global_data::GlobalData) where T <: RingElem
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

    # check the orders
    for xR in odd_images
        @assert is_infinity(ladder(27, xR, a24))
        @assert !is_infinity(ladder(9, xR, a24))
    end

    return Montgomery_coeff(a24), (I_sec, D_sec, xP, xQ, xPQ, odd_images)
end

function commitment(global_data::GlobalData)
    D_sec = random_secret_prime()
    a24, xP, xQ, xPQ, I_sec = RandIsogImages(D_sec, global_data, false)
    a24, (xP, xQ, xPQ) = Montgomery_normalize(a24, [xP, xQ, xPQ])
    A = Montgomery_coeff(a24)
    xPc, xQc, xPQc = torsion_basis(a24, SQISIGN_challenge_length)
    xP = xDBLe(xP, a24, ExponentFull - SQISIGN_challenge_length)
    xQ = xDBLe(xQ, a24, ExponentFull - SQISIGN_challenge_length)
    xPQ = xDBLe(xPQ, a24, ExponentFull - SQISIGN_challenge_length)
    n1, n2, n3, n4 = ec_bi_dlog_commitment(A, xP, xQ, xPQ, xPc, xQc, xPQc, global_data.E0_data)
    M = [n1 n3; n2 n4]
    return A, (I_sec, D_sec, xP, xQ, xPQ), M
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
    Isec, Dsec, xPsec, xQsec, xPQsec, odd_images = sk

    # compute commitment
    Acom, (Icom, Dcom, xPcom, xQcom, xPQcom), Mcom = commitment(global_data)

    # compute challenge and the pull-back of the corresponding ideal
    cha = challenge(Acom, m, global_data)
    a, b = Mcom * [1, cha]
    a, b, c, d = global_data.E0_data.Matrix_2ed_inv * [b, 0, -a, 0]
    alpha = SQIsign2D.Level1.QOrderElem(a, b, c, d)
    Icha = SQIsign2D.Level1.LeftIdeal(alpha, BigInt(1) << SQISIGN_challenge_length)

    # find alpha in bar(Isec)IcomIcha suitable for the response
    Icomcha = SQIsign2D.Level1.intersection(Icom, Icha)
    I = SQIsign2D.Level1.involution_product(Isec, Icomcha)
    nI = Dsec * Dcom << SQISIGN_challenge_length
    alpha, c, d, found = element_for_response(I, nI, ExponentForTorsion, [3, 3 ,3], Dsec)

    if found
        q = div(norm(alpha), d*nI)
        n_odd_l = length(global_data.E0_data.DegreesOddTorsionBases)
        odd_kernels = Proj1{FqFieldElem}[]
        ns = 1
        for i in 1:n_odd_l
            l = global_data.E0_data.DegreesOddTorsionBases[i]
            e = global_data.E0_data.ExponentsOddTorsionBases[i]
            if d % l == 0
                f = Int(log(l, gcd(d, l^e)))
                while alpha % l == Quaternion_0
                    alpha = div(alpha, l)
                    ns *= l
                    f -= 2
                end
                if f > 0
                    xPodd, xQodd, xPQodd = odd_images[3*(i-1)+1:3*i]
                    xPodd = ladder(l^(e-f), xPodd, a24pub)
                    xQodd = ladder(l^(e-f), xQodd, a24pub)
                    xPQodd = ladder(l^(e-f), xPQodd, a24pub)
                    Kodd = kernel_generator(xPodd, xQodd, xPQodd, a24pub, alpha, l, f, global_data.E0_data.Matrices_odd[i])
                else
                    Kodd = infinity_point(global_data.Fp2)
                end
                push!(odd_kernels, Kodd)
            end
        end

        # compute the auxiliary ellitic curve
        a24aux, xPaux, xQaux, xPQaux, images = auxiliary_path(a24pub, xPsec, xQsec, xPQsec, odd_kernels, Isec, Dsec, q, c, global_data)
        dd = div(d, ns^2)
        xPaux, xQaux, xPQaux = ladder(ns, xPaux, a24aux), ladder(ns, xQaux, a24aux), ladder(ns, xPQaux, a24aux)
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
                @assert is_infinity(ladder(l, K, a24aux))
                a24aux, tmp = odd_isogeny(a24aux, K, l, vcat([xPaux, xQaux, xPQaux], images))
                xPaux, xQaux, xPQaux = tmp[1:3]
                images = tmp[4:end]
            end
        end
        @assert is_infinity(xDBLe(xPaux, a24aux, c))
        @assert is_infinity(xDBLe(xQaux, a24aux, c))
        @assert is_infinity(xDBLe(xPQaux, a24aux, c))
        @assert !is_infinity(xDBLe(xPaux, a24aux, c-1))
        @assert !is_infinity(xDBLe(xQaux, a24aux, c-1))
        @assert !is_infinity(xDBLe(xPQaux, a24aux, c-1))

        return Acom, Montgomery_coeff(a24aux), xPaux, xQaux, xPQaux, odd_kernels, q, c, d, true
    end
    return nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, false
end

