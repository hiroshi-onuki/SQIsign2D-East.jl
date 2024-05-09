using SQIsign2D
_, E0_data = SQIsign2D.Level1.make_E0_data()

e_full = SQIsign2D.Level1.ExponentFull
chall_len = SQIsign2D.Level1.SQISIGN_challenge_length

for _ in 1:10
    Epub, I, nI, xPpub, xQpub, xPQpub = SQIsign2D.Level1.key_gen(E0_data)
    Ecom, Icom, nIcom, xPcom, xQcom, xPQcom = SQIsign2D.Level1.key_gen(E0_data)
    c = SQIsign2D.Level1.challenge(Ecom, "", E0_data)

    a24_com = A_to_a24(Ecom)
    xPc, xQc, xPQc = torsion_basis(a24_com, chall_len)
    xPcom_s = xDBLe(xPcom, a24_com, e_full - chall_len)
    xQcom_s = xDBLe(xQcom, a24_com, e_full - chall_len)
    xPQcom_s = xDBLe(xPQcom, a24_com, e_full - chall_len)
    n1, n2, n3, n4 = SQIsign2D.Level1.ec_bi_dlog_commitment(Ecom, xPcom_s, xQcom_s, xPQcom_s, xPc, xQc, xPQc, E0_data)


end