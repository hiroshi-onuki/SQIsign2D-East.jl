using SHA

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
    E_pub, xP, xQ, xPQ, I_sec = RandIsogImages(D_sec, global_data, true)
    return E_pub, I_sec, xP, xQ, xPQ, D_sec
end

function commitment(global_data::E0Data)
    D_com = random_secret_prime()
    E_com, xP, xQ, xPQ, I_sec = RandIsogImages(D_com, global_data, true)
    return E_com, I_sec, xP, xQ, xPQ, D_sec
end

function challenge(A::Proj1{T}, global_data::E0Data)
    h = sha3_256(string(A) * m)

    c = BigInt(0)
    len = ExponentFull - ExpornentForTorsion
    n, r = divrem(len, 8)
    for i in 1:n
        c += BigInt(h[i]) << (8*(i-1))
    end
end
