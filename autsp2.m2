-- find automorphisms of E which extend to automorphisms of P^2
-- f3 and f4 for automorphism x -> t1 - x
-- f5 and f4 for x -> t2 - x
-- f6 and f4 for x -> t3 - x

R = QQ[a_(0,0)..a_(2,2)]

A = transpose(genericMatrix(R, 3, 3))
B1 = (A * matrix{ {0}, {1}, {0} })
B2 = (A * matrix{ {0}, {0}, {1} })
B3 = (A * matrix{ {2}, {0}, {1} })
B4 = (A * matrix{ {1}, {0}, {1} })
I1 = minors(2, B1 | matrix{ {0}, {0}, {1}}) -- p -> s1
I2 = minors(2, B2 | matrix{ {0}, {1}, {0}}) -- s1 -> p
I3 = minors(2, B3 | matrix{ {1}, {0}, {1}}) -- s3 -> s2
I4 = minors(2, B4 | matrix{ {2}, {0}, {1}}) -- s2 -> s3

