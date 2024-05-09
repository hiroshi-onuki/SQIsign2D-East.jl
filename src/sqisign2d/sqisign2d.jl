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

function key_gen(global_data::E0Data)
    D_sec = random_secret_prime()
    a24, xP, xQ, xPQ, I_sec = RandIsogImages(D_sec, global_data, true)
    a24, (xP, xQ, xPQ) = Montgomery_normalize(a24, [xP, xQ, xPQ])
    return Montgomery_coeff(a24), (I_sec, D_sec, xP, xQ, xPQ)
end

function commitment(global_data::E0Data)
    A, sk = key_gen(global_data)
    _, _, xP, xQ, xPQ = sk
    a24 = A_to_a24(A)
    xPc, xQc, xPQc = torsion_basis(a24, SQISIGN_challenge_length)
    xP = xDBLe(xP, a24, ExponentFull - SQISIGN_challenge_length)
    xQ = xDBLe(xQ, a24, ExponentFull - SQISIGN_challenge_length)
    xPQ = xDBLe(xPQ, a24, ExponentFull - SQISIGN_challenge_length)
    n1, n2, n3, n4 = ec_bi_dlog_commitment(A, xP, xQ, xPQ, xPc, xQc, xPQc, global_data)
    M = [n1 n3; n2 n4]
    return A, sk, M
end


function challenge(A::FqFieldElem, m::String, global_data::E0Data)
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

function signing(pk::FqFieldElem, sk, m::String, global_data::E0Data)
    A = pk
    a24pub = A_to_a24(A)
    Isec, Dsec, xPpub, xQpub, xPQpub = sk
    Acom, (Icom, Dcom, xPcom, xQcom, xPQcom), Mcom = commitment(global_data)

    c = challenge(Acom, m, global_data)

    a, b = Mcom * [1, c]
    a, b, c, d = global_data.Matrix_2ed_inv * [b, 0, -a, 0]
    alpha = SQIsign2D.Level1.QOrderElem(a, b, c, d)
    Icha = SQIsign2D.Level1.LeftIdeal(alpha, BigInt(1) << SQISIGN_challenge_length)

    Icomcha = SQIsign2D.Level1.intersection(Icom, Icha)
    I = SQIsign2D.Level1.involution_product(Isec, Icomcha)
    nI = Dsec * Dcom << SQISIGN_challenge_length
    alpha, c, d, found = element_for_response(I, nI, ExponentForTorsion, [3, 5], Dsec)
    return c, d, found
end

