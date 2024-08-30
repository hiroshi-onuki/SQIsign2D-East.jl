export integer_square_root, integer_to_bytes, bytes_to_integer, sqrt_mod_2power

# floor(sqrt(n))
function integer_square_root(n::T) where T <: Integer
    n == 0 && return T(0)
    x = n + 1
    y = n
    while y < x
        x = y
        y = div(x + div(n, x), 2)
    end
    return x
end

function integer_to_bytes(n::Integer, length::Integer)
    bytes = Vector{UInt8}(undef, length)
    for i in 1:length
        bytes[i] = n >> (8*(length - i)) & 0xff
    end
    return bytes
end

# the inverse of the above function
function bytes_to_integer(bytes::Vector{UInt8})
    n = BigInt(0)
    for byte in bytes
        n = (n << 8) | byte
    end
    return n
end

# square root of n modulo 2^e
function sqrt_mod_2power(n::BigInt, e::Int)
    x = 1
    N = BigInt(1) << 3
    for i in 4:e
        N <<= 1
        (x^2 - n) % N != 0 && (x += N >> 2)
    end
    return x
end