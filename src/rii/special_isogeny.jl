
function isogeny_from_special_product(P1P2::CouplePoint{T}, Q1Q2::CouplePoint{T}, PQ1PQ2::CouplePoint{T}, E0_data::E0Data, image_points::Vector{CouplePoint{T}}, n::Int) where T <: RingElem
    a24 = E0_data.a24_0
    T1T2 = double_iter(P1P2, a24, a24, e-1)
    if (T1T2.P1.X == 0 && T1T2.P2.X != 0) || (T1T2.P1.X != 0 && T1T2.P2.X == 0)
        return nothing
    end
    S1S2 = double_iter(Q1Q2, a24, a24, e-1)
    if (S1S2.P1.X == 0 && S1S2.P2.X != 0) || (S1S2.P1.X != 0 && S1S2.P2.X == 0)
        return nothing
    end
    if T1T2.P1 == T1T2.P2 && S1S2.P1 == S1S2.P2
        # kernel is {(R, R) R in E0[2]}
        ret = CouplePoint{T}[]
        for R1R2 in image_points
            
            if R1R2.P1 == T1T2.P1 && R1R2.P2 == S1S2.P1
                return nothing
            end
        end
    end

end