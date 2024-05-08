using SQIsign2D

_, E0_data = SQIsign2D.Level1.make_E0_data()
for i in 1:100
    SQIsign2D.Level1.RandIsogImages(BigInt(1241231298203467), E0_data)
    println(i)
end