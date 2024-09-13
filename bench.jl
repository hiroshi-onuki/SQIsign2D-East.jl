using Nemo
import Pkg
Pkg.activate(@__DIR__)
using SQIsign2D

function benchmark_test(param::Module, num::Int, is_compact::Bool)
    global_data = param.make_precomputed_values()

    # for compilation
    pk, sk = param.key_gen(global_data)
    m = "message to sign"
    sign = param.signing(pk, sk, m, global_data, is_compact)
    if !is_compact
        verif = param.verify(pk, sign, m, global_data)
    else
        varif = param.verify_compact(pk, sign, m, global_data)
    end

    t_gen = 0
    t_sign = 0
    t_verif = 0
    println("Benchmark test for $(param) $((is_compact ? "compact" : "")) start")
    for i in 1:num
        (pk, sk), t, _, _, _ = @timed param.key_gen(global_data)
        t_gen += t

        m = "message to sign"
        sign, t, _, _, _ = @timed param.signing(pk, sk, m, global_data, is_compact)
        t_sign += t

        if !is_compact
            verif, t, _, _, _ = @timed param.verify(pk, sign, m, global_data)
            t_verif += t
        else
            verif, t, _, _, _ = @timed param.verify_compact(pk, sign, m, global_data)
            t_verif += t
        end

        @assert verif
        print("\r($i/$num) done.")
    end
    print("\n")
    println("Average time for key generation: ", t_gen / num)
    println("Average time for signing: ", t_sign / num)
    println("Average time for verification: ", t_verif / num)
    print("\n")
end

num = 100
benchmark_test(SQIsign2D.Level1, num, false)
benchmark_test(SQIsign2D.Level1, num, true)
benchmark_test(SQIsign2D.Level3, num, false)
benchmark_test(SQIsign2D.Level3, num, true)
benchmark_test(SQIsign2D.Level5, num, false)
benchmark_test(SQIsign2D.Level5, num, true)
