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

function auxiliary_path(a24::Proj1{T}, xP::Proj1{T}, xQ::Proj1{T}, xPQ::Proj1{T}, I::LeftIdeal, nI::BigInt,
                        q::BigInt, c::Int, global_data::GlobalData) where T <: RingElem
    d = q * ((BigInt(1) << c) - q)
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

    a24aux, xPaux, xQaux, xPQaux, images = d2isogeny(a24, a24d, xP, xQ, xPQ, xPd, xQd, xPQd, c, q, Proj1{FqFieldElem}[], global_data)

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

    return Montgomery_coeff(a24), (I_sec, D_sec, xP, xQ, xPQ)
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
    Isec, Dsec, xPsec, xQsec, xPQsec = sk

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
        
    
        a24aux, xPaux, xQaux, xPQaux, images = auxiliary_path(a24pub, xPsec, xQsec, xPQsec, Isec, Dsec, q, c, global_data)
    end
    return c, d, found
end

