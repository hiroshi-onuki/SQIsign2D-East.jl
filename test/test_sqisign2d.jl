using SQIsign2D

function check(param::Module, num::Int, is_compact::Bool)
    global_data = param.make_precomputed_values()
    for i in 1:num
        pk, sk = param.key_gen(global_data)
        m = "hello"
        sign = param.signing(pk, sk, m, global_data, is_compact)
        i == 1 && println("sign len: ", length(sign))
        if !is_compact
            @assert param.verify(pk, sign, m, global_data)
        else
            @assert param.verify_compact(pk, sign, m, global_data)
        end
    end
end

#check(SQIsign2D.Level1, 100, false)
#check(SQIsign2D.Level1, 100, true)
#check(SQIsign2D.Level3, 100, false)
#check(SQIsign2D.Level3, 100, true)
check(SQIsign2D.Level5, 100, false)
check(SQIsign2D.Level5, 100, true)
