include("level1/prime.jl")
include("level1/constants.jl")

include("../quaternion/order.jl")
include("../quaternion/cornacchia.jl")
include("../quaternion/ideal.jl")
include("../quaternion/klpt.jl")

include("global_data.jl")

include("../elliptic_curves/dlog.jl")
include("../rii/quat_action.jl")
include("../rii/d2isogeny.jl")
include("../rii/rii.jl")
include("../utilities/for_compression.jl")
include("../sqisign2d/sqisign2d.jl")

StrategyChallenge = compute_strategy(div(SQISIGN_challenge_length, 2) - 1, 1, 1)

const StrategiesDim2 = Dict(
    ExponentFull => compute_strategy(ExponentFull-2, 2, 1),
    ExponentFull-1 => compute_strategy(ExponentFull-3, 2, 1),
    ExponentFull-2 => compute_strategy(ExponentFull-4, 2, 1),
    ExponentFull-3 => compute_strategy(ExponentFull-5, 2, 1),
    ExponentFull-4 => compute_strategy(ExponentFull-6, 2, 1),
    ExponentFull-5 => compute_strategy(ExponentFull-7, 2, 1),
    ExponentForTorsion => compute_strategy(ExponentForTorsion-2, 2, 1)
)

function make_precomputed_values()
    _, T = polynomial_ring(GF(p), "T")
    Fp2, Fp2_i = finite_field(T^2 + 1, "i")

    A0 = Fp2(0)

    # constatns from precompute/level1torsion.sage
    P2e = Point(316623125071715518772158340021045282959479262663943283253543156979590879681937*Fp2_i + 126014172852369434262151758375756662898243266424750582023096561213827028220808, 341706843175207212785549465416045120066109766272815901261142766545754177295340*Fp2_i + 157384708942011414488406961093001680377352947193020742482870421848371076661789)
    Q2e = Point(74175176104226640782393734383276406045306935582593620379626189047115932852846*Fp2_i + 264784128323572725292400316028565026106542931821786321610072784812879784313975, 157384708942011414488406961093001680377352947193020742482870421848371076661789*Fp2_i + 49091458000734946769002608988276568938676431973721002372026579480952635239443)
    M_i_2e = [0 14474011154664524427946373126085988481658748083205070504932198000989141204991; 1 0]
    M_ij_2e = [13551627187396668629744157997382890111556584978623772704122722262201241471332 4363428768907237784387094857376449674935230259589858134694636875197695799311; 4363428768907237784387094857376449674935230259589858134694636875197695799312 922383967267855798202215128703098370102163104581297800809475738787899733660]
    M_1k_2e = [10110582385757286643559278268709538806723517823615212370237561125791445405681 13551627187396668629744157997382890111556584978623772704122722262201241471332; 13551627187396668629744157997382890111556584978623772704122722262201241471332 4363428768907237784387094857376449674935230259589858134694636875197695799312]
    M44inv = BigInt[27445936713386080586175713661918785552 75334857984268849943091041782148331164 75334857984268849943091041782148331164 57624655016848535279667938196023267313; 75334857984268849943091041782148331164 57624655016848535279667938196023267313 57624655016848535279667938196023267312 9735733745965765922752610075793721700; 19471467491931531845505220151587443400 54891873426772161172351427323837571103 54891873426772161172351427323837571103 65599124238303084020338431706354609464; 30178718303462454693492224534104481761 19471467491931531845505220151587443400 19471467491931531845505220151587443400 54891873426772161172351427323837571103]
    P3 = Point(71428194509393699580050219996362648677259120466327103178011703226172995427130*Fp2_i + 189870362096014084817533691451031521581100286578627598598640358389670292262548, 29239894617524128590023903517008385749718093073792154822488100145278003406417*Fp2_i + 290110459003652114366579003210248426450355831916748350534412193446526156835924)
    Q3 = Point(74360760996323894980339737470266779337689959973463792034765462360845485885175*Fp2_i + 195669981078806460791590142331658076427261987639973205277662314191173644180146, 343444816337799379055043162838699909716554641950230906479035464922009651512108*Fp2_i + 260134696925494400784580590983388889729399381669209389967125251870137560453135)
    M_i_3 = [3 8; 19 24]
    M_ij_3 = [12 18; 19 15]
    M_1k_3 = [0 12; 9 1]

    a24_0 = A_to_a24(A0)
    xP2e = Proj1(P2e.X, P2e.Z)
    xQ2e = Proj1(Q2e.X, Q2e.Z)
    PQ2e = add(P2e, -Q2e, Proj1(A0))
    xPQ2e = Proj1(PQ2e.X, PQ2e.Z)

    # precomputed values for discrete logarithm
    window_size = 3
    gen = 1
    while gen^(BigInt(2)^(ExponentFull - 1)) == 1
        gen = rand(Fp2)^((p^2 - 1) >> ExponentFull)
    end
    base = gen^(BigInt(2)^(ExponentFull - SQISIGN_challenge_length))
    fq_dlog_table1_c, fq_dlog_table2_c = make_dlog_table(base, SQISIGN_challenge_length, window_size)
    strategy_dlog_c = compute_strategy(div(SQISIGN_challenge_length, window_size) - 1, window_size, 1)
    dlog_data_chall = DlogData(SQISIGN_challenge_length, window_size, fq_dlog_table1_c, fq_dlog_table2_c, strategy_dlog_c)
    base = gen^(BigInt(2)^(ExponentFull - ExponentForTorsion))
    fq_dlog_table1_res, fq_dlog_table2_res = make_dlog_table(base, ExponentForTorsion, window_size)
    strategy_dlog_res = compute_strategy(div(ExponentForTorsion, window_size) - 1, window_size, 1)
    dlog_data_res = DlogData(ExponentForTorsion, window_size, fq_dlog_table1_res, fq_dlog_table2_res, strategy_dlog_res)

    DegreesOddTorsionBases = [(3, 3)]
    OddTorsionBases = [[P3, Q3]]

    Matrices_2e = [M_i_2e, M_ij_2e, M_1k_2e]
    Matrices_odd = [[M_i_3, M_ij_3, M_1k_3]]

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

    return GlobalData(Fp2, Fp2_i, E0Data(A0, A0d, A0dd, a24_0, jInvariant_A(A0), P2e, Q2e, xP2e, xQ2e, xPQ2e, DegreesOddTorsionBases, OddTorsionBases, Matrices_2e, M44inv, Matrices_odd, w, isomorphism_to_A0, dlog_data_chall, dlog_data_res))
end
