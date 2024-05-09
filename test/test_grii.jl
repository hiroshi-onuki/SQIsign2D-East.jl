using SQIsign2D
_, E0_data = SQIsign2D.Level1.make_E0_data()

e_full = SQIsign2D.Level1.ExponentFull
chall_len = SQIsign2D.Level1.SQISIGN_challenge_length

for _ in 1:100
    pk, sk = SQIsign2D.Level1.key_gen(E0_data)
    ret = SQIsign2D.Level1.signing(pk, sk, "hello", E0_data)
    println(ret)
end