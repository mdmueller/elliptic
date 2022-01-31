-- find automorphisms x -> s_i - x, s_i are 2-torsion translates of q = [121/16:165/64:1]

loadPackage "MultiprojectiveVarieties"

R = QQ[a_(0,0)..a_(3,3)]

A = transpose(genericMatrix(R, 4, 4))
q = transpose(matrix{ {0,0,0,1} })
p = transpose(matrix{ {0,0,1,0} })
t1 = transpose(matrix{ {-121/16, 0, 0, 165/64} })
t2 = transpose(matrix{ {1,6,0,(165/64)/(6-121/6)} })
s1 = transpose(matrix{ {-44/15, -896/55, 672/121, 1} })
s2 = transpose(matrix{ {-20/33, -72/55, -96/25, 1} })
I1 = minors(2, (A*p) | q)
I2 = minors(2, (A*q) | p)
I3 = minors(2, (A*t1) | s1)
I4 = minors(2, (A*s1) | t1)
I5 = minors(2, (A*t2) | s2)
I6 = minors(2, (A*s2) | t2)
I = I1+I2+I3+I4

