export get_basis, integral_LLL, HNFmod, Gauss_elimination_mod,
    closest_vector, enumerate_close_vector, short_basis

# Return a Z-module basis from input generators gens
function get_basis(gens::Vector{Vector{T}}) where T <: Integer
    n = length(gens)
    n == 0 && return gens
    m = length(gens[1])
    prod([length(v) == m for v in gens]) || error("lengths of generatars are different")
    
    basis = deepcopy(gens)
    i = m
    j = n
    while true
        while j < 1 || prod([b[i] == 0 for b in basis[1:j]])
            i -= 1
            i == max(m - n, 0) && return filter(x->x!=zeros(m), basis)
        end
        while j > 1 && !prod([b[i] == 0 for b in basis[1:j-1]])
            k = findfirst(x->x==minimum(abs.([b[i] for b in filter(x->x[i]!=0, basis[1:j])])), abs.([b[i] for b in basis[1:j]]))
            basis[j], basis[k] = basis[k], basis[j]
            basis[j][i] < 0 && (basis[j] *= -1)
            for k in 1:j-1
                basis[k] -= div(basis[k][i], basis[j][i]) * basis[j]
            end
        end
        j -= 1
    end
end

# Algorithm 2.6.7 in H. Cohen, A Course in Computational Algebraic Number Theory
function integral_LLL(basis::Vector{Vector{T}}, quadratic_form::Function) where T <: Integer
    b = deepcopy(basis)

    # input check
    n = length(b)
    n == 0 && return b
    m = length(b[1])
    prod([length(v) == m for v in b]) || error("lengths of generatars are different")
    m < n && error("number of vectors greater than these dimension")

    q = quadratic_form
    k = 2
    kmax = 1
    d = zeros(T, n+1)
    d[1] = 1
    d[2] = q(b[1], b[1])
    H = zeros(T, n, n)
    for i in 1:n H[i,i] = 1 end
    lam = zeros(T, n, n)

    while k <= n
        if k > kmax
            kmax = k
            for j in 1:k
                u = q(b[k], b[j])
                for i in 1:j-1
                    u = div(d[i+1]*u - lam[k, i]*lam[j, i], d[i])
                end
                if j < k
                    lam[k, j] = u
                elseif j == k
                    d[k+1] = u
                end
            end
            d[k+1] == 0 && error("vectors are not linearly independent")
        end
        b, H, lam = REDI(k, k-1, b, d, H, lam)
        while d[k+1]*d[k-1] < 3//4*d[k]^2 - lam[k, k-1]^2
            b, d, H, lam = SWAPI(k, kmax, b, d, H, lam)
            k = max(2, k-1)
            b, H, lam = REDI(k, k-1, b, d, H, lam)
        end
        for l in 1:k-2
            b, H, lam = REDI(k, k-1-l, b, d, H, lam)
        end
        k += 1
    end
    return H
end

function REDI(k::Int, l::Int, b::Vector{Vector{T}}, d::Vector{T}, H::Matrix{T}, lam::Matrix{T}) where T <: Integer
    if abs(2*lam[k, l]) > d[l+1]
        q = T(round(lam[k,l]//d[l+1]))
        H[:, k] -= q*H[:, l]
        b[k] -= q*b[l]
        lam[k,l] -= q*d[l+1]
        for i in 1:l-1
            lam[k,i] -= q*lam[l,i]
        end
    end
    return b, H, lam
end

function SWAPI(k::Int, kmax::Int, b::Vector{Vector{T}}, d::Vector{T}, H::Matrix{T}, lam::Matrix{T}) where T <: Integer
    H[:, k], H[:, k-1] = H[:, k-1], H[:, k]
    b[k], b[k-1] = b[k-1], b[k]
    if k > 2
        for j in 1:k-2
            lam[k, j], lam[k-1, j] = lam[k-1, j], lam[k, j]
        end
    end
    lamd = lam[k,k-1]
    B = T((d[k-1]*d[k+1] + lamd^2) // d[k])
    for i in k+1:kmax
        t = lam[i,k]
        lam[i,k] = T((d[k+1]*lam[i,k-1] - lamd*t) // d[k])
        lam[i,k-1] = T((B*t + lamd*lam[i,k]) // d[k+1])
    end
    d[k] = B
    return b, d, H, lam
end

#=
    return a list of coefficients of short vectors in alattice.
    basis: a basis of the target lattice, M: bilinear matrix,
    C; upper bound of the quadratic forms of output vectors
    Algorithm 2.7.5 in H. Cohen, A Course in Computational Algebraic Number Theory.
=#
function short_vectors(basis::Vector{Vector{T}}, quadratic_form::Function, C::T) where T <: Integer
    # input check
    n = length(basis)
    n == 0 && return basis
    m = length(basis[1])
    prod([length(v) == m for v in basis]) || error("lengths of generatars are different")
    m < n && error("number of vectors greater than these dimension")

    # LLL reduction
    H = integral_LLL(basis, quadratic_form)
    LLLmat = hcat([b for b in basis]...) * H
    red_basis = [LLLmat[:, i] for i in 1:n]

    q = make_quadratic_form_coeffs(red_basis, M)
    S = zeros(Rational{T}, n)
    U = zeros(Rational{T}, n)
    L = zeros(T, n)
    x = zeros(T, n)
    S[n] = C
    out_vecs = []

    i = n
    tmp = T(floor(S[i] // q[i,i] * denominator(U[i])^2))
    Z = integer_square_root(tmp) // denominator(U[i])
    L[i] = T(floor(Z - U[i]))
    x[i] = T(ceil(-Z-U[i]) - 1)

    while true
        x[i] += 1
        while x[i] > L[i]
            i += 1
            x[i] += 1
        end
        if i > 1
            S[i-1] = S[i] - q[i,i]*(x[i] + U[i])^2
            i -= 1
            U[i] = sum([q[i,j]*x[j] for j in i+1:n])

            tmp = T(floor(S[i] // q[i,i] * denominator(U[i])^2))
            Z = IntegerSquareRoot(tmp) // denominator(U[i])
            L[i] = T(floor(Z - U[i]))
            x[i] = T(ceil(-Z-U[i]) - 1)
        else
            if x != zeros(T, n)
                v = sum([x[i]*red_basis[i] for i in 1:n])
                push!(out_vecs, [v, T(C - S[1] + q[1,1]*(x[1] + U[1])^2)])
            else
                return out_vecs 
            end
        end
    end
end

# Is LLL reduced?
function LLLcheck(b::Vector{Vector{T}}, M::Matrix{T}) where T <: Integer
    # input check
    n = length(b)
    n == 0 && error("no input vecotor")
    m = length(b[1])
    prod([length(v) == m for v in b]) || error("lengths of generatars are different")
    m == size(M,1) == size(M,2) || error("sizes of vector and bilinear matrix are different")
    m < n && error("number of vectors greater than these dimension")

    q(x, y) = transpose(x)*M*y
    mu = zeros(Rational{T}, n, n)
    bs = [zeros(Rational{T}, n, n)[:,i] for i in 1:n]
    for i in 1:n
        bs[i] = copy(b[i])
        for j in 1:i-1
            mu[i,j] = q(b[i],bs[j]) // q(bs[j],bs[j])
            abs(mu[i,j]) > 1//2 && return false
            bs[i] -= mu[i,j]*bs[j]
        end
        i > 1 && q(bs[i], bs[i]) < (3//4 - mu[i,i-1]^2) * q(bs[i-1],bs[i-1]) && return false
    end
    return true
end

# Compute Hermite normal form mod D. Algorithm 2.4.8 in H. Cohen, A Course in Computational Algebraic Number Theory.
function HNFmod(M::Matrix{T}, D::T) where T <: Integer
    m, n = size(M)
    i = 1
    k = 1
    l = min(m,n)
    M = M .% D

    while true
        while M[i,k+1:end] .% D != zeros(T, n-k)
            j = findfirst(x->x==minimum(abs.(filter(x->x!=0, M[i,k:end]))), abs.(M[i,k:end])) + k - 1
            M[:, j], M[:, k] = M[:, k], M[:, j]
            M[i, k] < 0 && (M = -M .% D)
            b = M[i, k]
            for j in k+1:n
                q = div(M[i, j], b)
                M[:, j] = (M[:, j] - q*M[:, k]) .% D
            end
        end
        b = M[i, k]
        if b == 0
            k -= 1
        else
            for j in 1:k-1
                q = div(M[i, j], b)
                M[:, j] = (M[:, j] - q*M[:, k]) .% D
            end
        end
        i == l && return M
        i += 1
        k += 1
    end
end

# Gauss elimination modulo a prime N by row transformations
function Gauss_elimination_mod(M::Matrix{T}, N::Integer) where T <: Integer
    M = mod.(M, N)
    m, n = size(M)
    r = 1
    for i in 1:n
        k = 0
        for j in r:m
            M[j,i] != 0 && (k = j)
        end
        if k != 0
            M[k, :], M[r, :] = M[r, :], M[k, :]
            M[r, :] = mod.(M[r, :] * invmod(M[r, i], N), N)
            for j in 1:m
                j != i && (M[j, :] = mod.(M[j, :] - M[j, i]*M[r, :], N))
            end
            r += 1
        end
    end
    return M
end

# Algorithm 5 in SQIsign documentation
function short_basis(b0::Vector{T}, b1::Vector{T}) where T <: Integer
    Nq(x) = x[1]^2 + x[2]^2
    Bq(x, y) = x[1]*y[1] + x[2]*y[2]
    Nq(b0) < Nq(b1) && ((b0, b1) = (b1, b0))
    gamma = b0
    while true
        r = T(round(Bq(b0, b1) / Nq(b1)))
        gamma = b0 - r * b1
        if Nq(gamma) < Nq(b1)
            (b0, b1) = (b1, gamma)
        else
            break
        end
    end
    Nq(gamma) < Nq(b0) && (b0 = gamma)
    return b1, b0
end

# Algorithm 6 in SQIsign documentation
# the closest vector of a dim-2 vector t in a lattice <b1, b0>, where Nq(b1) < Nq(b0).
function closest_vector(b0::Vector{T}, b1::Vector{T}, t::Vector{T}) where T <: Integer
    Nq(x) = x[1]^2 + x[2]^2
    Bq(x, y) = x[1]*y[1] + x[2]*y[2]
    mu = Nq(b1) * b0 - Bq(b0, b1) * b1
    c = t - T(round(Bq(mu, t) * Nq(b1) // Nq(mu))) * b0
    c -= T(round(Bq(b1, c) // Nq(b1))) * b1
    return t - c
end

# Algorithm 7 in SQIsign documentation
# enumerate vectors close to a dim-2 vector t in a lattice <b1, b0>,
# vc in <b1, b0> close to t, m is maximal number of tries, B is a bound of the norm
function enumerate_close_vector(b1::Vector{T}, b0::Vector{T}, t::Vector{T}, vc::Vector{T}, m::Int, B::T) where T <: Integer
    Nq(x) = x[1]^2 + x[2]^2
    Bq(x, y) = x[1]*y[1] + x[2]*y[2]

    mu = Nq(b1) * b0 - Bq(b0, b1) * b1
    ret = Vector{T}[]

    i = 0
    d = t - vc
    a, b, c = Nq(b0), 2*Bq(b0, b1), Nq(b1)
    Be = B
    B - Nq(d) > 0 && (Be = B - Nq(d))
    By = T(floor((sqrt(4*a^2*Be) + 1)/sqrt(4*a^2*c - b^2))) + 1
    y = -By - 1
    while y < By && i < m
        y += 1
        Bx = T(floor((2*a*(1 + sqrt(4*a^2*Be + 4*c*a^2*y^2 - b^2*y^2)) - b*y*sqrt(4*a^3)) / (2*a*sqrt(4*a^3)))) + 1
        x = -T(floor((2*a*(1 + sqrt(4*a^2*Be + 4*c*a^2*y^2 - b^2*y^2)) + b*y*sqrt(4*a^3)) / (2*a*sqrt(4*a^3)))) - 2
        while x < Bx && i < m
            x += 1
            i += 1
            Nq(d - x*b0 - y*b1) <= B && push!(ret, vc + x*b0 + y*b1)
        end
    end
    return ret
end

