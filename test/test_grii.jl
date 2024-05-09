using SQIsign2D
_, E0_data = SQIsign2D.Level1.make_E0_data()

e_full = SQIsign2D.Level1.ExponentFull
chall_len = SQIsign2D.Level1.SQISIGN_challenge_length

for _ in 1:10
    Epub, I, nI, xPpub, xQpub, xPQpub = SQIsign2D.Level1.key_gen(E0_data)
    Ecom, Icom, nIcom, xPcom, xQcom, xPQcom = SQIsign2D.Level1.key_gen(E0_data)
    c = SQIsign2D.Level1.challenge(Ecom, "", E0_data)

    a24_pub = A_to_a24(Epub)
    xPc, xQc, xPQc = torsion_basis(a24_pub, chall_len)
    xPpub_s = xDBLe(xPpub, a24_pub, e_full - chall_len)
    xQpub_s = xDBLe(xQpub, a24_pub, e_full - chall_len)
    xPQpub_s = xDBLe(xPQpub, a24_pub, e_full - chall_len)
    n1, n2, n3, n4 = SQIsign2D.Level1.ec_bi_dlog_pubkey(Epub, xPpub_s, xQpub_s, xPQpub_s, xPc, xQc, xPQc, E0_data)
    Mpub = [n1 n3; n2 n4]

    @assert xPpub_s == linear_comb_2_e(n1, n2, xPc, xQc, xPQc, a24_pub, chall_len)

    a, b = Mpub * [1, c]
    a, b, c, d = E0_data.Matrix_2ed_inv * [b, 0, -a, 0]
    alpha = SQIsign2D.Level1.QOrderElem(a, b, c, d)
    Icha = SQIsign2D.Level1.LeftIdeal(alpha, BigInt(1) << chall_len)
    

end