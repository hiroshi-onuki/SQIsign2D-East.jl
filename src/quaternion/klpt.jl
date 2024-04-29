# Algorithm 10 in SQIsign documentation
# return a quaternion in O0 of norm M
function FullRepresentInteger(M::Integer)
    counter = 0
    found = false
    x, y, z, w = 0, 0, 0, 0
    while !found && counter < KLPT_repres_num_gamma_trial
        m = integer_square_root(div(4*M, p))
        z = rand(-m:m)
        md = integer_square_root(div(4*M - z^2, p))
        w = rand(-md:md)
        Md = 4*M - p*(z^2 + w^2)
        x, y, found = sum_of_two_squares(Md)
        if !found || (x - w) % 2 != 0 || (y - z) % 2 != 0
            found = false
            counter += 1
        end
    end
    if found
        return QOrderElem(div(x - w, 2), div(y - z, 2), z, w), found
    else
        return Quaternion_0, found
    end
end

# Algorithm 8 in SQIsign documentation
# return J ~ I s.t. norm(J) ~ sqrt(p) and prime
function RandomEquivalentPrimeIdeal(I::LeftIdeal)
    counter = 0
    found = false
    J = LeftIdeal([Quaternion_0, Quaternion_0, Quaternion_0, Quaternion_0])
    nJ = 0
    N = norm(I)

    # LLL reduction
    Imatrix = ideal_to_matrix(I)
    q(x, y) = quadratic_form(QOrderElem(x), QOrderElem(y))
    H = integral_LLL([Imatrix[:, i] for i in 1:4], q)
    LLLmat = Imatrix * H
    red_basis = [LLLmat[:, i] for i in 1:4]

    while !found && counter < KLPT_equiv_num_iter
        counter += 1
        c1, c2, c3, c4 = [rand(-KLPT_equiv_bound_coeff:KLPT_equiv_bound_coeff) for _ in 1:4]
        beta = c1 * red_basis[1] + c2 * red_basis[2] + c3 * red_basis[3] + c4 * red_basis[4]
        beta = QOrderElem(beta)
        nJ = div(norm(beta), N)
        if is_probable_prime(nJ)
            found = true
            J = ideal_transform(I, beta, N)
        end
    end
    return J, nJ, found
end

# Isec: a left ideal of O0, I: a left ideal of O0 ofr norm N
# return [Isec]^* L, where L ~ (\bar(Isec)I) and norm(L) is prime
function RandomEquivalentPrimeIdeal_for_signing(Isec::LeftIdeal, I::LeftIdeal, Nsec::BigInt, NI::BigInt)
    counter = 0
    found = false
    J = LeftIdeal([Quaternion_0, Quaternion_0, Quaternion_0, Quaternion_0])
    nJ = 0
    N = Nsec * NI

    # compute \bar(Isec)I
    invIsecI = involution_product(Isec, I)

    # LLL reduction
    Imatrix = ideal_to_matrix(invIsecI)
    q(x, y) = quadratic_form(QOrderElem(x), QOrderElem(y))
    H = integral_LLL([Imatrix[:, i] for i in 1:4], q)
    LLLmat = Imatrix * H
    red_basis = [LLLmat[:, i] for i in 1:4]

    beta = Quaternion_0
    while !found && counter < KLPT_equiv_num_iter
        counter += 1
        c1, c2, c3, c4 = [rand(-KLPT_equiv_bound_coeff:KLPT_equiv_bound_coeff) for _ in 1:4]
        beta = c1 * red_basis[1] + c2 * red_basis[2] + c3 * red_basis[3] + c4 * red_basis[4]
        beta = QOrderElem(beta)
        nJ = div(norm(beta), N)
        if is_probable_prime(nJ)
            found = true
            J = ideal_transform(I, beta, NI) # Isec * L
            J = larger_ideal(J, nJ) # [Isec]^* L
        end
    end
    return J, nJ, beta, found
end

# Algorithm 11 in SQIsign documentation
# return C, D s.t. gamma*j*(C + D*i)*delta in Z + I, where N = n(I), if divisible then N | norm(gamma).
function EichlerModConstraint(I::LeftIdeal, N::Integer, gamma::QOrderElem, delta::QOrderElem, divisible::Bool)
    M = HNFmod(ideal_to_matrix(I), N)
    v1 = to_vector(gamma * Quaternion_j)
    v2 = to_vector(-gamma * Quaternion_ij * delta)

    if divisible
        # gamma*j*(C + D*i)*delta in I
        M = hcat(v1, v2, M[:, 1:2])
    else
        # gamma*j*(C + D*i)*delta in Z + I
        M = hcat(v1, v2, [1, 0, 0, 0], M[:, 1:2])
    end
    M = Gauss_elimination_mod(M, N)
    if M[1, 1] == 0
        C, D = 1, 0
    elseif M[2, 2] == 0
        C, D = 0, 1
    elseif M[3, 3] == 0
        i = findfirst(x!=0 for x in M[1, 3:end])
        if i != nothing
            C, D = M[1, i + 2], M[2, i + 2]
        else
            C, D = 0, 1
        end
    elseif M[4, 4] == 0
        i = findfirst(x!=0 for x in M[1, 4:end])
        if i != nothing
            C, D = M[1, i + 3], M[2, i + 3]
        else
            C, D = 0, 1
        end
    else
        C, D = M[1, 5], M[2, 5]
    end
    return C, D
end

# Algorithm 13 in SQIsign documentation
function FullStrongApproximation(N::Integer, C::Integer, D::Integer, lambda::Integer, N_mu::Integer, max_tries::Int)
    c0 = 2*p*lambda*C
    c1 = invmod(2*p*lambda*D, N)
    c2 = div(N_mu - lambda^2*p*(C^2 + D^2), N)
    c3 = (-c0*c1) % N
    b0 = [0, N^2]
    b1 = [N, c3*N]
    b1, b0 = short_basis(b0, b1)
    v = [lambda*C, lambda*D + N*c1*c2]
    vd = closest_vector(b0, b1, v)
    B = BigInt(2) << Int(min(ceil(log(2, N_mu/p)), 3*ceil(log(2, N)) + 10))
    vs = enumerate_close_vector(b1, b0, v, vd, max_tries, B)
    for vdd in vs
        M = div(N_mu - p*((lambda*C - vdd[1])^2 + (lambda*D + N*c1*c2 - vdd[2])^2), N^2)
        a, b, found = sum_of_two_squares(M)
        if found
            mu = order_elem_from_standard_basis(N*a, N*b, lambda*C - vdd[1], -lambda*D - N*c1*c2 + vdd[2])
            if gcd(mu) == 2
                return div(mu, 2), true
            end
        end
    end
    return Quaternion_0, false
end

# Algorithm 14 in SQIsign documentation
# return alpha in I s.t. norm(alpha)/N_I = 2^KLPT_keygen_length * ExtraDegree
function KeyGenKLPT(I::LeftIdeal, N_I::Integer)
    counter = 0
    found = false
    k = Log2p - Int(floor(log(2, N_I))) + KLPT_gamma_exponent_center_shift
    N_gamma = BigInt(2)^k
    N_mu = BigInt(2)^(KLPT_keygen_length - k)

    gamma = Quaternion_0
    mu = Quaternion_0
    while !found && counter < KLPT_keygen_num_gamma_trial
        counter += 1

        gamma, found_gamma = FullRepresentInteger(N_I*N_gamma)
        !found_gamma && continue

        C, D = EichlerModConstraint(I, N_I, gamma, QOrderElem(1), true)
        N_CD = p * (C^2 + D^2)
        N_mu_N_CD = (N_mu * invmod(N_CD, N_I)) % N_I
        quadratic_residue_symbol(N_mu_N_CD, N_I) != 1 && continue
        lambda = sqrt_mod(4*N_mu_N_CD, N_I)

        mu, found = FullStrongApproximation(N_I, C, D, lambda, 4*N_mu, KLPT_signing_number_strong_approx)
    end
    return gamma*mu, found
end

# KLPT algorihtm for an ideal I of prime norm N_I.
# This algorihtm is for computing an isogeny between E0 and another curve whose endomorphism ring p-etremal.
function KLPT(I::LeftIdeal, N_I::Integer)
    counter = 0
    found = false
    n = Int(ceil(log(2, N_I)))
    k = max(Log2p - n, 0)
    (n + Log2p + k) % 2 == 1 && (k += 1)
    N_gamma = BigInt(2)^k * ExtraDegree
    N_mu = BigInt(2)^(Log2p + 3*n)

    gamma = Quaternion_0
    mu = Quaternion_0
    while !found && counter < KLPT_keygen_num_gamma_trial
        counter += 1

        gamma, found_gamma = FullRepresentInteger(N_I*N_gamma)
        !found_gamma && continue

        C, D = EichlerModConstraint(I, N_I, gamma, QOrderElem(1), true)
        N_CD = p * (C^2 + D^2)
        N_mu_N_CD = (N_mu * invmod(N_CD, N_I)) % N_I
        quadratic_residue_symbol(N_mu_N_CD, N_I) != 1 && continue
        lambda = sqrt_mod(4*N_mu_N_CD, N_I)

        mu, found = FullStrongApproximation(N_I, C, D, lambda, 4*N_mu, KLPT_signing_number_strong_approx)
    end
    return gamma*mu, found
end

# Algorithm 17 in SQIsign documentation
function SigningKLPT(Isec::LeftIdeal, I::LeftIdeal, Nsec::BigInt, N_I::BigInt)
    L, NL, beta, found = RandomEquivalentPrimeIdeal_for_signing(Isec, I, Nsec, N_I)
    @assert found
    k = Log2p - Int(floor(log(2, NL))) + KLPT_gamma_exponent_center_shift
    N_gamma = BigInt(2)^k
    N_mu = BigInt(2)^(KLPT_signing_klpt_length - k)

    counter = 0
    found = false
    gamma_mu = Quaternion_0
    while !found && counter < KLPT_signing_num_gamma_trial
        counter += 1

        gamma, gamma_found = FullRepresentInteger(NL * N_gamma)
        !gamma_found && continue

        C0, D0 = EichlerModConstraint(L, NL, gamma, QOrderElem(1), true)
        N_CD = p * (C0^2 + D0^2)
        N_mu_N_CD = (N_mu * invmod(N_CD, NL)) % NL
        quadratic_residue_symbol(N_mu_N_CD, NL) != 1 && continue
        lambda0 = sqrt_mod(N_mu_N_CD, NL)

        C1, D1 = EichlerModConstraint(Isec, Nsec, gamma, QOrderElem(1), false)
        N_CD = p * (C1^2 + D1^2)
        N_mu_N_CD = (N_mu * invmod(N_CD, Nsec)) % Nsec
        quadratic_residue_symbol(N_mu_N_CD, Nsec) != 1 && continue
        lambda1 = sqrt_mod(N_mu_N_CD, Nsec)

        lam = 2 * crt([lambda0, lambda1], [NL, Nsec])
        C = crt([C0, C1], [NL, Nsec])
        D = crt([D0, D1], [NL, Nsec])

        mu, found = FullStrongApproximation(Nsec*NL, C, D, lam, 4*N_mu, KLPT_signing_number_strong_approx)
        gamma_mu = gamma * mu
        (trace(gamma_mu) % 2 == 0 || gcd(gamma_mu*beta) % 2 == 0)&& (found = false)
    end
    return ideal_transform(L, gamma_mu, NL), found
end