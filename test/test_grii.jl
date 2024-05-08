using SQIsign2D
_, E0_data = SQIsign2D.Level1.make_E0_data()

a24, I, nI, xP, xQ, xPQ = SQIsign2D.Level1.key_gen(E0_data)
d = rand(1:nI^3)
SQIsign2D.Level1.GeneralizedRandomIsogImages(d, a24, I, nI, E0_data)