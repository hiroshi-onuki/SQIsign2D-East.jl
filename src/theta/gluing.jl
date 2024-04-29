
function get_base_submatrix(a24::Proj1{T}, P::Proj1{T}) where T <: RingElem
    P2 = xDBL(P, a24)
    x, z = P.X, P.Z
    u, w = P2.X, P2.Z

    wx = w * x
    wz = w * z
    ux = u * x
    uz = u * z
    det = wx - uz

    inverse = batched_inversion([det, z])

    d = uz * inverse[1]
    a = -d
    b = -wz * inverse[1]
    c = ux * inverse[1] - x * inverse[2]

    return [a, b, c, d]
end

function get_base_matrix(a24_1::Proj1{T}, a24_2::Proj1{T}, T1::CouplePoint{T}, T2::CouplePoint{T}) where T <: RingElem
    F = parent(a24_1.X)

    P1, P2 = T1.P1, T1.P2
    Q1, Q2 = T2.P1, T2.P2
    
    g00_1, g01_1, g10_1, g11_1 = get_base_submatrix(a24_1, P1)
    g00_2, g01_2, g10_2, g11_2 = get_base_submatrix(a24_2, P2)
    h00_1, _, h10_1, _ = get_base_submatrix(a24_1, Q1)
    h00_2, h01_2, h10_2, h11_2 = get_base_submatrix(a24_2, Q2)

    gh00_1 = g00_1 * h00_1 + g01_1 * h10_1
    gh10_1 = g10_1 * h00_1 + g11_1 * h10_1
    gh00_2 = g00_2 * h00_2 + g01_2 * h10_2
    gh10_2 = g10_2 * h00_2 + g11_2 * h10_2
    
    a, b, c, d = F(1), F(0), F(0), F(0)

    # T1
    a += g00_1 * g00_2
    b += g00_1 * g10_2
    c += g10_1 * g00_2
    d += g10_1 * g10_2

    # T2
    a += h00_1 * h00_2
    b += h00_1 * h10_2
    c += h10_1 * h00_2
    d += h10_1 * h10_2

    # T1+T2
    a += gh00_1 * gh00_2
    b += gh00_1 * gh10_2
    c += gh10_1 * gh00_2
    d += gh10_1 * gh10_2

    # Now we act by (0, Q2)
    a1 = h00_2 * a + h01_2 * b
    b1 = h10_2 * a + h11_2 * b
    c1 = h00_2 * c + h01_2 * d
    d1 = h10_2 * c + h11_2 * d

    # Now we act by (P1, 0)
    a2 = g00_1 * a + g01_1 * c
    b2 = g00_1 * b + g01_1 * d
    c2 = g10_1 * a + g11_1 * c
    d2 = g10_1 * b + g11_1 * d

    # Now we act by (P1, Q2)
    a3 = g00_1 * a1 + g01_1 * c1
    b3 = g00_1 * b1 + g01_1 * d1
    c3 = g10_1 * a1 + g11_1 * c1
    d3 = g10_1 * b1 + g11_1 * d1

    return [a, b, c, d, a1, b1, c1, d1, a2, b2, c2, d2, a3, b3, c3, d3]
end

function apply_base_chagne(P::ThetaPtLv2{T}, M::Vector) where T <: RingElem
    x, y, z, t = P.a, P.b, P.c, P.d
    a = M[1] * x + M[2] * y + M[3] * z + M[4] * t
    b = M[5] * x + M[6] * y + M[7] * z + M[8] * t
    c = M[9] * x + M[10] * y + M[11] * z + M[12] * t
    d = M[13] * x + M[14] * y + M[15] * z + M[16] * t
    return ThetaPtLv2(a, b, c, d)
end

function apply_base_chagne(P::ThetaNullLv2{T}, M::Vector) where T <: RingElem
    x, y, z, t = P.a, P.b, P.c, P.d
    a = M[1] * x + M[2] * y + M[3] * z + M[4] * t
    b = M[5] * x + M[6] * y + M[7] * z + M[8] * t
    c = M[9] * x + M[10] * y + M[11] * z + M[12] * t
    d = M[13] * x + M[14] * y + M[15] * z + M[16] * t
    return ThetaNullLv2(a, b, c, d)
end

function base_change_couple_point(P::CouplePoint{T}, M::Vector{T}) where T <: RingElem
    P1, P2 = P.P1, P.P2
    X1, Z1 = P1.X, P1.Z
    X2, Z2 = P2.X, P2.Z

    X = X1 * X2;
    Y = X1 * Z2;
    Z = Z1 * X2;
    W = Z1 * Z2;
    return apply_base_chagne(ThetaPtLv2(X, Y, Z, W), M)
end

# codomain of the gluing (2, 2)-isogeny with kernel <4*T1, 4*T2>
function gluing_codomain(T1::ThetaPtLv2{T}, T2::ThetaPtLv2{T}) where T <: RingElem
    F = parent(T1.a)

    xAxByCyD = Hadamard(square(T1))
    zAzBtYtD = Hadamard(square(T2))

    zero_idx = 1
    while xAxByCyD[zero_idx] != 0
        zero_idx += 1
    end
    zero_idx -= 1

    t1 = zAzBtYtD[(1 ⊻ zero_idx) + 1]
    t2 = xAxByCyD[(2 ⊻ zero_idx) + 1]
    t3 = zAzBtYtD[(3 ⊻ zero_idx) + 1]
    t4 = xAxByCyD[(3 ⊻ zero_idx) + 1]
    inverse = batched_inversion([t1, t2, t3, t4])
    
    ABCD = zeros(F, 4)
    ABCD[(0 ⊻ zero_idx) + 1] = F(0)
    ABCD[(1 ⊻ zero_idx) + 1] = t1 * inverse[3]
    ABCD[(2 ⊻ zero_idx) + 1] = t2 * inverse[4]
    ABCD[(3 ⊻ zero_idx) + 1] = F(1)

    a_inverse = t3 * inverse[1]
    b_inverse = t4 * inverse[2]

    A, B, C, D = Hadamard(ABCD)
    return ThetaNullLv2(A, B, C, D), a_inverse, b_inverse, zero_idx
end

# the iamge of P under a gluing (2, 2)-isogeny, where PT is P + T with T one of generators of the kernel
function gluing_image(P::ThetaPtLv2{T}, PT::ThetaPtLv2{T}, a_inv::T, b_inv::T, zero_idx::Integer) where T <: RingElem
    F = parent(P.a)
    AxByCzDt = Hadamard(square(P))
    AyBxCtDz = Hadamard(square(PT))

    y = AxByCzDt[(1 ⊻ zero_idx) + 1] * a_inv
    z = AxByCzDt[(2 ⊻ zero_idx) + 1] * b_inv
    t = AxByCzDt[(3 ⊻ zero_idx) + 1]

    zb = AyBxCtDz[(3 ⊻ zero_idx) + 1]
    tb = AyBxCtDz[(2 ⊻ zero_idx) + 1] * b_inv

    if z != 0
        lam = z / zb
    else
        lam = t / tb
    end

    xb = AyBxCtDz[(1 ⊻ zero_idx) + 1] * a_inv
    x = xb * lam

    xyzt = zeros(F, 4)
    xyzt[(0 ⊻ zero_idx) + 1] = x
    xyzt[(1 ⊻ zero_idx) + 1] = y
    xyzt[(2 ⊻ zero_idx) + 1] = z
    xyzt[(3 ⊻ zero_idx) + 1] = t

    return ThetaPtLv2(Hadamard(xyzt))
end

# the gluing (2, 2)-isogeny with kernel <(P1, P2), (Q1, Q2)>
function gluing_isogeny(a24_1::Proj1{T}, a24_2::Proj1{T},
        T1_8::CouplePoint{T}, T2_8::CouplePoint{T}, P1Q1P2Q2::CouplePoint{T},
        image_points::Vector{CouplePoint{T}}, n::Integer) where T <: RingElem
    T1_4 = double(T1_8, a24_1, a24_2)
    T2_4 = double(T2_8, a24_1, a24_2)
    M = get_base_matrix(a24_1, a24_2, T1_4, T2_4)

    T1 = base_change_couple_point(T1_8, M)
    T2 = base_change_couple_point(T2_8, M)

    codomain, a_inv, b_inv, zero_idx = gluing_codomain(T1, T2)

    images = Vector{ThetaPtLv2{T}}(undef, length(image_points))
    for i in 1:length(image_points)
        P = image_points[i]

        P1, P2 = P.P1, P.P2

        # the last two points are the images of generators of the kernel
        if i == length(image_points) - 1
            # P1 = 2^n*T1_4.P1, P2 = 2^n*T2_4.P2
            PT1 = ladder(BigInt(2)^n + 1, P1, a24_1)
            PT2 = ladder(BigInt(2)^n + 1, P2, a24_2)
        elseif i == length(image_points)
            # P1 + T1_4.P1 = P1 + 2^n * image_points[end-1].P1,
            # P2 + T1_4.P2 = P2 + 2^n * image_points[end-1].P2
            P1Q1, P2Q2 = P1Q1P2Q2.P1, P1Q1P2Q2.P2
            PT1 = ladder3pt(BigInt(2)^n, P1, image_points[end-1].P1, P1Q1, a24_1)
            PT2 = ladder3pt(BigInt(2)^n, P2, image_points[end-1].P2, P2Q2, a24_2)
        else
            # require two squre roots
            PT1 = x_add_sub(P1, T1_4.P1, a24_1)
            PT2 = x_add_sub(P2, T1_4.P2, a24_2)
        end
        PT = CouplePoint(PT1, PT2)
        Ptheta = base_change_couple_point(P, M)
        PTtheta = base_change_couple_point(PT, M)
        
        Pimage = gluing_image(Ptheta, PTtheta, a_inv, b_inv, zero_idx)
        images[i] = Pimage
    end

    return codomain, images
end

