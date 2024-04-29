
const chi_eval = Dict{Tuple{Int, Int}, Int}(
    (0, 0) => 1,
    (0, 1) => 1,
    (0, 2) => 1,
    (0, 3) => 1,
    (1, 0) => 1,
    (1, 1) => -1,
    (1, 2) => 1,
    (1, 3) => -1,
    (2, 0) => 1,
    (2, 1) => 1,
    (2, 2) => -1,
    (2, 3) => -1,
    (3, 0) => 1,
    (3, 1) => -1,
    (3, 2) => -1,
    (3, 3) => 1
)

const even_indices = [
    [0, 0],
    [0, 1],
    [0, 2],
    [0, 3],
    [1, 0],
    [1, 2],
    [2, 0],
    [2, 1],
    [3, 0],
    [3, 3]
]

const splitting_map = Dict{Tuple{Int, Int}, Vector{Int}}(
    (0, 2) => [1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, -1, 0],
    (3, 3) => [1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1],
    (0, 3) => [1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, -1],
    (2, 1) => [1, 1, 1, 1, 1, -1, 1, -1, 1, -1, -1, 1, 1, 1, -1, -1],
    (0, 1) => [1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, -1, 0, 0],
    (1, 2) => [1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0],
    (2, 0) => [1, 1, 1, 1, 1, -1, 1, -1, 1, -1, -1, 1, -1, -1, 1, 1],
    (3, 0) => [1, 1, 1, 1, 1, -1, 1, -1, 1, 1, -1, -1, -1, 1, 1, -1],
    (1, 0) => [1, 1, 1, 1, 1, -1, -1, 1, 1, 1, -1, -1, -1, 1, -1, 1]
)

function level_22_constants_sqr(tnull::ThetaNullLv2, chi::Integer, i::Integer)
    Ucontant = 0
    for t in 0:3
        Ucontant += chi_eval[(chi, t)] * tnull[t + 1] * tnull[(i ⊻ t) + 1]
    end
    return Ucontant
end

function identify_even_index(null_point::ThetaNullLv2)
    for (chi, i) in even_indices
        U_sqr = level_22_constants_sqr(null_point, chi, i)
        if U_sqr == 0
            return chi, i
        end
    end
    error("Not a product of elliptic_curves")
end

function compute_splitting_matrix(tnull::ThetaNullLv2)
    chi, i = identify_even_index(tnull)
    if chi == 0 && i == 0
        zeta = gen(parent(tnull[1]))    # zeta = sqrt(-1)
        return [1, zeta, 1, zeta, 1, -zeta, -1, zeta, 1, zeta, -1, -zeta, -1, zeta, -1, zeta]
    else
        return splitting_map[(chi, i)]
    end
end

function splitting_isomorphism(tnull::ThetaNullLv2{T}, image_points::Vector{ThetaPtLv2{T}}) where T <: RingElem
    M = compute_splitting_matrix(tnull)
    tnull = apply_base_chagne(tnull, M)
    
    ret = Vector{ThetaPtLv2{T}}(undef, length(image_points))
    for i in 1:length(image_points)
        ret[i] = apply_base_chagne(image_points[i], M)
    end
    return tnull, ret
end

function split_theta_point(P::ThetaLv2)
    a, b, d = P.a, P.b, P.d
    return [ThetaDim1(a, b), ThetaDim1(b, d)]
end

# theta null corresponding to a diagonal period matrix to Montgomery curves
function split_to_product(tnull::ThetaNullLv2{T}, image_points::Vector{ThetaPtLv2{T}}) where T <: RingElem
    O1, O2 = split_theta_point(tnull)
    E1 = theta_to_Montgomery(O1)
    E2 = theta_to_Montgomery(O2)

    images = Vector{CouplePoint{T}}(undef, length(image_points))
    for i in 1:length(image_points)
        P1, P2 = split_theta_point(image_points[i])
        xP1 = theta_point_to_Montgomery(O1, P1)
        xP2 = theta_point_to_Montgomery(O2, P2)
        images[i] = CouplePoint(xP1, xP2)
    end

    return [E1, E2], images
end