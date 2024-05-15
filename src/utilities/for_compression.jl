function Fq_to_bytes(x::FqFieldElem)
    bytes_len = div(Log2p, 8) + 1
    a, b = BigInt(lift(ZZ, coeff(x, 0))), BigInt(lift(ZZ, coeff(x, 1)))
    return vcat(integer_to_bytes(a, bytes_len), integer_to_bytes(b, bytes_len))
end

function bytes_to_Fq(bytes::Vector{UInt8}, Fq::Field)
    bytes_len = div(Log2p, 8) + 1
    a = bytes_to_integer(bytes[1:bytes_len])
    b = bytes_to_integer(bytes[bytes_len+1:2*bytes_len])
    return Fq(a) + Fq(b)*gen(Fq)
end