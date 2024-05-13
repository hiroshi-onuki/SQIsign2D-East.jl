using SQIsign2D
global_data = SQIsign2D.Level1.make_precomputed_values()

e_full = SQIsign2D.Level1.ExponentFull
chall_len = SQIsign2D.Level1.SQISIGN_challenge_length

cnt = 0
for _ in 1:100
    pk, sk = SQIsign2D.Level1.key_gen(global_data)
    m = "hello"
    sign, found = SQIsign2D.Level1.signing(pk, sk, m, global_data)
    if found
        global cnt += 1
        println(SQIsign2D.Level1.verify(pk, sign, m, global_data))
    end
    println("found: ", found)
end
println(cnt)