struct DlogData
    e::Int
    window_size::Int
    T1::Vector{Vector{FqFieldElem}}
    T2::Vector{Vector{FqFieldElem}}
    strategy::Vector{Int}
end

# structure for precomputed values
struct E0Data{FqFieldElem <: RingElem}
    A0::FqFieldElem
    A0d::FqFieldElem
    A0dd::FqFieldElem
    a24_0::Proj1{FqFieldElem}
    j0::FqFieldElem
    P2e::Point{FqFieldElem}
    Q2e::Point{FqFieldElem}
    xP2e::Proj1{FqFieldElem}
    xQ2e::Proj1{FqFieldElem}
    xPQ2e::Proj1{FqFieldElem}
    Matrices_2e::Vector{Matrix{BigInt}}
    Matrix_2ed_inv::Matrix{BigInt}
    isomorphism_to_A0::Function
    dlog_data_full::DlogData
    dlog_data_chall::DlogData
    tate_table::Vector{Vector{FqFieldElem}}
end
