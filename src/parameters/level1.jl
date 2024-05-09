include("level1/prime.jl")
include("level1/klpt_constants.jl")

include("../quaternion/order.jl")
include("../quaternion/cornacchia.jl")
include("../quaternion/ideal.jl")
include("../quaternion/klpt.jl")

include("global_data.jl")

include("../elliptic_curves/dlog.jl")
include("../rii/rii.jl")
include("../sqisign2d/sqisign2d.jl")

const StrategiesDim2 = Dict(
    ExponentFull => compute_strategy(ExponentFull-2, 2, 1),
    ExponentFull-1 => compute_strategy(ExponentFull-3, 2, 1),
    ExponentFull-2 => compute_strategy(ExponentFull-4, 2, 1),
    ExponentFull-3 => compute_strategy(ExponentFull-5, 2, 1),
    ExponentFull-4 => compute_strategy(ExponentFull-6, 2, 1),
    ExponentFull-5 => compute_strategy(ExponentFull-7, 2, 1)
)

function make_E0_data()
    _, T = polynomial_ring(GF(p), "T")
    Fp2, Fp2_i = finite_field(T^2 + 1, "i")

    A0 = Fp2(0)

    # constatns from precompute/level1torsion.sage
    P2e = Point(90915926031238959600021966182790273474590221405448741162072108620129331769439*Fp2_i + 1331831870677393868798541274452788250548570814894740416092387468610971155686, 11425578182889347356787341719085235637662500290587984622496057771050050020109*Fp2_i + 99213047816543008104484745736407684862998470109110190454228972499617487688574)
    Q2e = Point(32113168783409498037522205388940628619509137301794358129851574388278368472992*Fp2_i + 121697262943971063768745630297278113843550787892348358875831295539796729086745, 99213047816543008104484745736407684862998470109110190454228972499617487688574*Fp2_i + 111603516631759110280756829852645666456436858416655114669427625237357650222322)
    M_i_2e = [0 7237005577332262213973186563042994240829374041602535252466099000494570602495; 1 0]
    M_ij_2e = [2765022130526353248995196748405302602234424172011193721451773769877439650028 7086431377706405983283581248380179825110081150604599905442226954277120595312; 7086431377706405983283581248380179825110081150604599905442226954277120595313 4471983446805908964977989814637691638594949869591341531014325230617130952468]
    M_1k_2e = [150574199625856230689605314662814415719292890997935347023872046217450007184 2765022130526353248995196748405302602234424172011193721451773769877439650028; 2765022130526353248995196748405302602234424172011193721451773769877439650028 7086431377706405983283581248380179825110081150604599905442226954277120595313]
    M44inv = [25020076890174559025267592708576204145 8242280446756730142174088273756153620 8242280446756730142174088273756153620 17515218974942748907654233220394822288; 8242280446756730142174088273756153620 17515218974942748907654233220394822288 17515218974942748907654233220394822287 34293015418360577790747737655214872812; 26050734971603847648573649381458719192 7504857915231810117613359488181381857 7504857915231810117613359488181381857 16484560893513460284348176547512307240; 35030437949885497815308466440789644575 26050734971603847648573649381458719192 26050734971603847648573649381458719192 7504857915231810117613359488181381857]
    
    a24_0 = A_to_a24(A0)
    xP2e = Proj1(P2e.X, P2e.Z)
    xQ2e = Proj1(Q2e.X, Q2e.Z)
    PQ2e = add(P2e, -Q2e, Proj1(A0))
    xPQ2e = Proj1(PQ2e.X, PQ2e.Z)

    # precomputed values for discrete logarithm
    tp_table = make_pairing_table(A0, P2e, ExponentFull)
    tp_P2e_Q2e = Tate_pairing_P0(Q2e, tp_table, Cofactor)
    window_size = 3
    fq_dlog_table1, fq_dlog_table2 = make_dlog_table(tp_P2e_Q2e, ExponentFull, window_size)
    strategy_dlog = compute_strategy(div(ExponentFull, window_size) - 1, window_size, 1)
    dlog_data_full = DlogData(ExponentFull, window_size, fq_dlog_table1, fq_dlog_table2, strategy_dlog)
    base = tp_P2e_Q2e^(BigInt(2)^(ExponentFull - SQISIGN_challenge_length))
    fq_dlog_table1_c, fq_dlog_table2_c = make_dlog_table(base, SQISIGN_challenge_length, window_size)
    strategy_dlog_c = compute_strategy(div(SQISIGN_challenge_length, window_size) - 1, window_size, 1)
    dlog_data_chall = DlogData(SQISIGN_challenge_length, window_size, fq_dlog_table1_c, fq_dlog_table2_c, strategy_dlog_c)

    Matrices_2e = [M_i_2e, M_ij_2e, M_1k_2e]

    w = Weil_pairing_2power(A0, P2e, Q2e, ExponentFull)

    # make constants for isomorphism to the curve E_A0
    _, T = polynomial_ring(Fp2, "T")
    As = roots((256 * (T^2 - 3)^3 - 1728 * (T^2 - 4))/T^2)
    A0d = As[1]
    beta = -A0d/3
    gamma = square_root(1 / (1 - 3*beta^2))
    A0dd = As[2]
    beta_d = -A0dd/3
    gamma_d = square_root(1 / (1 - 3*beta_d^2))
    function isomorphism_to_A0(A::Proj1{FqFieldElem}, Ps::Vector{Proj1{FqFieldElem}})
        if A == Proj1(A0)
            return Ps
        elseif A == Proj1(A0d)
            return [Proj1(gamma*(P.X - beta*P.Z), P.Z) for P in Ps]
        elseif A == Proj1(A0dd)
            return [Proj1(gamma_d*(P.X - beta_d*P.Z), P.Z) for P in Ps]
        else
            throw(ArgumentError("A is not A0d or A0dd"))
        end
    end

    return Fp2, E0Data(A0, A0d, A0dd, a24_0, jInvariant_A(A0), P2e, Q2e, xP2e, xQ2e, xPQ2e, Matrices_2e, M44inv, w, isomorphism_to_A0, dlog_data_full, dlog_data_chall, tp_table)
end
