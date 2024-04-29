export compute_strategy

# Algorithm 60 in SIKE documentation
function compute_strategy(n::Int, p::T, q::T) where T<:Real
    S = Vector{Int}[[]]
    C = T[0]
    for i in 2:n+1
        b = argmin([C[i-b] + C[b] + b*p + (i-b)*q for b in 1:i-1])
        push!(S, vcat([b],S[i-b],S[b]))
        push!(C, C[i-b] + C[b] + b*p + (i-b)*q)
    end
    return S[n+1]
end

