-- find automorphisms of E which extend to automorphisms of P^5
-- f3 and f4 for automorphism x -> t1 - x
-- f5 and f4 for x -> t2 - x
-- f6 and f4 for x -> t3 - x

R = QQ[a_(0,0)..a_(5,5)]

A = transpose(genericMatrix(R, 6, 6))
p = transpose(matrix{ {0,0,0,0,0,1} })
t1 = transpose(matrix{ {1,0,0,0,0,0} })
t2 = transpose(matrix{ {1,6,0,36,0,0} })
t3 = transpose(matrix{ {1,7,0,49,0,0} })
q1 = transpose(matrix{ {1,3,6,9,18,36} })
q2 = transpose(matrix{ {1,3,-6,9,-18,36} })
q3 = transpose(matrix{ {1,14,28,14^2,14*28,28^2} })
q4 = transpose(matrix{ {1,14,-28,14^2,-14*28,28^2} })
r1 = transpose(matrix{ {1,8,-4,8^2,-8*4,4^2} })
r2 = transpose(matrix{ {1,8,4,8^2,8*4,4^2} })
s1 = transpose(matrix{ {1,21/4,21/8,(21/4)^2,21/4*21/8,(21/8)^2} })
s2 = transpose(matrix{ {1,21/4,-21/8,(21/4)^2,-21/4*21/8,(21/8)^2} })

I1 = minors(2, (A*p) | t1)
I2 = minors(2, (A*t1) | p)
I3 = minors(2, (A*t2) | t3)
I4 = minors(2, (A*t3) | t2)
I5 = minors(2, (A*q1) | q3)
I6 = minors(2, (A*q2) | q4)
I7 = minors(2, (A*q3) | q1)
I8 = minors(2, (A*q4) | q2)
I = I1+I2+I3+I4+I5+I6+I7+I8 -- x |-> t1 - x

I1 = minors(2, (A*p) | t2)
I2 = minors(2, (A*t1) | t3)
I3 = minors(2, (A*t2) | p)
I4 = minors(2, (A*t3) | t1)
I5 = minors(2, (A*q1) | r1)
I6 = minors(2, (A*r1) | q1)
I7 = minors(2, (A*q2) | r2)
I8 = minors(2, (A*r2) | q2)
J = I1+I2+I3+I4+I5+I6+I7+I8 -- x |-> t2 - x

I1 = minors(2, (A*p) | t3)
I2 = minors(2, (A*t1) | t2)
I3 = minors(2, (A*t2) | t1)
I4 = minors(2, (A*t3) | p)
I5 = minors(2, (A*q1) | s1)
I6 = minors(2, (A*s1) | q1)
I7 = minors(2, (A*q2) | s2)
I8 = minors(2, (A*s2) | q2)
K = I1+I2+I3+I4+I5+I6+I7+I8 -- x |-> t3 - x

SOL1 = matrix{
    { 0, -1/42, 0, 13/1764, 0, 1/1764},
    { 0, 0, 0, 1/42, 0, 0},
    { 0, 0, 0, 0, 1/42, 0},
    { 0, 1, 0, 0, 0, 0},
    { 0, 0, 1, 0, 0, 0},
    {42, -13, 0, 1, 0, 0} }

SOL2 = matrix{
    {-6, 11/6, 0, -5/36, 0, 1/36},
    {-42, 13, 0, -1, 0, 1/6},
    {0, 0, 1, 0, -1/6, 0},
    {-294, 91, 0, -7, 0, 1},
    {0, 0, 7, 0, -1, 0},
    {0, -7, 0, 1, 0, 0} }

SOL3 = matrix{
    {-7, 15/7, 0, -8/49, 0, 1/49},
    {-42, 13, 0, -1, 0, 1/7},
    {0, 0, -1, 0, 1/7, 0},
    {-252, 78, 0, -6, 0, 1},
    {0, 0, -6, 0, 1, 0},
    {0, -6, 0, 1, 0, 0}}
