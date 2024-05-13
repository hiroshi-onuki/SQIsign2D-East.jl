import endomorphism as end
import torsion_points as tp

def matrix_to_str(M):
    return "[%s %s; %s %s]" % (M[0,0], M[0,1], M[1,0], M[1,1])

def matrix44_to_str(M):
    return "[%s %s %s %s; %s %s %s %s; %s %s %s %s; %s %s %s %s]" % (
        M[0,0], M[0,1], M[0,2], M[0,3],
        M[1,0], M[1,1], M[1,2], M[1,3],
        M[2,0], M[2,1], M[2,2], M[2,3],
        M[3,0], M[3,1], M[3,2], M[3,3]
    )

def make_constants(p, e, ed, degs, file_name):
    Fp4, Fp2, zeta4 = tp.calcFields(p)
    E0 = EllipticCurve(Fp4, [1, 0])
    Fp2d.<Fp2_i> = GF(p^2, modulus=x^2+1)

    out_file = open(file_name, "w")

    # 2^e-torsion in E(Fp2)
    P, Q = tp.basis_2e_special(E0, Fp2, zeta4, e)
    Ms = end.action_matrices([P, Q], 2^e, zeta4, Fp4)
    Px, Py = [tp.Fp2ToFp2d(v, zeta4, Fp2_i) for v in P.xy()]
    Qx, Qy = [tp.Fp2ToFp2d(v, zeta4, Fp2_i) for v in Q.xy()]
    out_file.write("P2e = Point(%s, %s)\n" % (Px, Py))
    out_file.write("Q2e = Point(%s, %s)\n" % (Qx, Qy))
    out_file.write("M_i_2e = %s\n" % matrix_to_str(Ms[0]))
    out_file.write("M_ij_2e = %s\n" % matrix_to_str(Ms[1]))
    out_file.write("M_1k_2e = %s\n" % matrix_to_str(Ms[2]))
    Msd = [[[1, 0], [0, 1]]] + Ms
    Z2e = quotient(ZZ, 2^ed)
    M44 = matrix(Z2e, [[Msd[i][j][k] for i in range(4)] for j, k in [(0, 0), (1, 0), (0, 1), (1, 1)]])
    out_file.write("M44inv = %s\n" % matrix44_to_str(M44^(-1)))

   # odd torsion in E(Fp2)
    for l, e in factor(degs):
        P, Q = tp.basis(E0, Fp2, False, l, e)
        Ms = end.action_matrices([P, Q], l^e, zeta4, Fp4)
        Px, Py = [tp.Fp2ToFp2d(v, zeta4, Fp2_i) for v in P.xy()]
        Qx, Qy = [tp.Fp2ToFp2d(v, zeta4, Fp2_i) for v in Q.xy()]
        out_file.write("P%d = Point(%s, %s)\n" % (l, Px, Py))
        out_file.write("Q%d = Point(%s, %s)\n" % (l, Qx, Qy))
        out_file.write("M_i_%d = %s\n" % (l, matrix_to_str(Ms[0])))
        out_file.write("M_ij_%d = %s\n" % (l, matrix_to_str(Ms[1])))
        out_file.write("M_1k_%d = %s\n" % (l, matrix_to_str(Ms[2])))

    out_file.close()

# level1
set_random_seed(0)
p = 2^253 * 3^3 - 1
e = 253
ed = 126
degs = 3^3
make_constants(p, e, ed, degs, "level1torsion.txt")

# level3
set_random_seed(0)
p = 2^380 * 35 - 1
e = 380
ed = 189
degs = 35
make_constants(p, e, ed, degs, "level3torsion.txt")

# level5
set_random_seed(0)
p = 2^520 * 2 - 1
e = 520
ed = 259
degs = 1
make_constants(p, e, ed, degs, "level5torsion.txt")