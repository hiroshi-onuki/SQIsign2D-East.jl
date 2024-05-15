using SQIsign2D

function check(param::Module, num::Int)
    global_data = param.make_precomputed_values()
    cnt = 0
    for _ in 1:num
        pk, sk = param.key_gen(global_data)
        m = "hello"
        sign, found = param.signing(pk, sk, m, global_data)
        if found
            println("sign len: ", length(sign))
            cnt += 1
            println(param.verify(pk, sign, m, global_data))
        end
        println("found: ", found)
    end
    println(cnt)
end

check(SQIsign2D.Level1, 10)
check(SQIsign2D.Level3, 10)