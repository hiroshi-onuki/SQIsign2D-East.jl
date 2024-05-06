# alpha(P), alpha(Q) for a fixed basis (P, Q) of E0[2^ExponentFull]
function action_on_torsion_basis(alpha::QOrderElem, E0_data::E0Data)
    xP, xQ, xPQ = E0_data.xP2e, E0_data.xQ2e, E0_data.xPQ2e
    Malpha = alpha[1] * [1 0; 0 1] + alpha[2] * E0_data.Matrices_2e[1] + alpha[3] * E0_data.Matrices_2e[2] + alpha[4] * E0_data.Matrices_2e[3]
    x_alphaP = linear_comb_2_e(Malpha[1, 1], Malpha[2, 1], xP, xQ, xPQ, E0_data.a24_0, ExponentFull)
    x_alphaQ = linear_comb_2_e(Malpha[1, 2], Malpha[2, 2], xP, xQ, xPQ, E0_data.a24_0, ExponentFull)
    x_alphaPQ = linear_comb_2_e(Malpha[1, 1] - Malpha[1, 2], Malpha[2, 1] - Malpha[2, 2], xP, xQ, xPQ, E0_data.a24_0, ExponentFull)
    return x_alphaP, x_alphaQ, x_alphaPQ
end

function RandIsogImages(d::BigInt, E0_data::E0Data)
    deg_dim2 = BigInt(2) << ExponentFull
    alpha = FullRepresentInteger(d*(deg_dim2 - d))
    xP, xQ, xPQ = action_on_torsion_basis(alpha, E0_data)


end