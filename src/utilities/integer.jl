export integer_square_root, integer_to_bytes, bytes_to_integer

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