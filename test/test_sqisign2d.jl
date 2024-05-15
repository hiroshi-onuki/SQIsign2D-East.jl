using SQIsign2D

function check(param::Module, num::Int)
    global_data = param.make_precomputed_values()
    for _ in 1:num
        pk, sk = param.key_gen(global_data)
        m = "hello"
        sign = param.signing(pk, sk, m, global_data)
        println("sign len: ", length(sign))
        @assert param.verify(pk, sign, m, global_data)
    end
end

check(SQIsign2D.Level1, 10)
check(SQIsign2D.Level3, 10)
check(SQIsign2D.Level5, 10)