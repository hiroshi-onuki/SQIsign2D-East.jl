using Nemo
import Pkg
Pkg.activate(@__DIR__)
using SQIsign2D

function check_exponents(a::Int, b::Int, is_odd::Bool=true)

    # compute an ideal corresponding to a reponse
    Fp2, global_data = SQIsign2D.Level1.make_E0_data()
    _, I_sec, D_sec, _, _, _ = SQIsign2D.Level1.key_gen(global_data)
    I = SQIsign2D.Level1.sample_random_ideal_2e(256)
    J = SQIsign2D.Level1.involution_product(I_sec, I)
    nJ = SQIsign2D.Level1.norm(J)

    # shortest vector
    alpha, e, found = SQIsign2D.Level1.element_for_response(J, nJ, nJ*BigInt(2)^a, 7)

    if found
        q = div(SQIsign2D.Level1.norm(alpha), e)
        return q !=0 && q*(BigInt(2)^a - q) < BigInt(2)^b
    else
        return false
    end
end

a = 127
b = 252
cnt = 0
for _ in 1:1000
    if check_exponents(a, b)
        global cnt += 1
    end
end
println(cnt)
