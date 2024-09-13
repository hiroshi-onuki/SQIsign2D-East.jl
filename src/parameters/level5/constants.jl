# basically the same as SQIsign reference implementation
const KLPT_equiv_bound_coeff = 7
const KLPT_equiv_num_iter = 50625
const KLPT_gamma_exponent_center_shift = 15
const KLPT_repres_num_gamma_trial = 32768
const KLPT_signing_number_strong_approx = 6784
const KLPT_secret_key_prime_size = 64
const KLPT_keygen_num_gamma_trial = 64
const KLPT_keygen_number_strong_approx = 5218
const SQISIGN_challenge_length = ExponentFull - ExponentForTorsion

const SQISIGN2D_Fp2_length = 130
const SQISIGN2D_2a_length = Int(ceil(ExponentForTorsion/8))
const SQISIGN2D_signature_length = 2 * SQISIGN2D_Fp2_length + 3 * SQISIGN2D_2a_length + 4
const CompactSQISIGN2D_signature_length = SQISIGN2D_Fp2_length + 5 * SQISIGN2D_2a_length + 6

const SmallPrimes = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97]
const FactorForAuxiliaryDegree = 3
const FactorInTwist = true