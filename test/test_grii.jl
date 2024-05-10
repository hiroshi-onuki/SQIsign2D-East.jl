using SQIsign2D
global_data = SQIsign2D.Level1.make_precomputed_values()

e_full = SQIsign2D.Level1.ExponentFull
chall_len = SQIsign2D.Level1.SQISIGN_challenge_length

cnt = 0
for _ in 1:10
    pk, sk = SQIsign2D.Level1.key_gen(global_data)
    c, d, found = SQIsign2D.Level1.signing(pk, sk, "hello", global_data)
    if found
        global cnt += 1
    end
    println(c, " ", d, " ", found)
end
println(cnt)