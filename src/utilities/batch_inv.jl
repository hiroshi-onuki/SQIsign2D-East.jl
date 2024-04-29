"""
Given k values, compute the inverse using
3(k-1) multiplications and one inversion
"""
function batched_inversion(input::Vector{T}) where T <: RingElem
    # for values a0, a1, ..., an compute
    # a0, a0a1, a0a1a2, ... a0...an using
    # (k-1) multiplications
    multiples = T[input[1]]
    for ai in input[2:end]
        push!(multiples, ai * multiples[end])
    end

    # Compute 1 / (a0a1a2...an)
    last_multiple = multiples[end]
    inverses_multiples = T[1 / last_multiple]

    # Compute (a0a1a2...an)^-1, (a0a1a2...a(n-1)) ... a0^-1
    # using k-1 multiplications
    for ai in reverse(input[2:end])
        push!(inverses_multiples, ai * inverses_multiples[end])
    end

    # Reverse for easy ordering below
    # inverses_multiples = a0^-1, (a0a1)^-1 ...(a0a1...an)^-1
    inverses_multiples = reverse(inverses_multiples)
    inverses = [inverses_multiples[1]]

    # Compute the inverse of each element from multiples and
    # their inverses using k-1 multiplications
    k = length(input)
    for i in 2:k
        push!(inverses, inverses_multiples[i] * multiples[i - 1])
    end

    return inverses
end