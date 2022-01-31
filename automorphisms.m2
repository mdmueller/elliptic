-- find automorphisms of E which extend to automorphisms of P^3
-- f3 and f4 for automorphism x -> t1 - x
-- f5 and f4 for x -> t2 - x
-- f6 and f4 for x -> t3 - x

R = QQ[A11, A12, A13, A14, A21, A22, A23, A24, A31, A32, A33, A34, A41, A42, A43, A44]
S = R[x,y,z,w]
phi = map(R, S)

A = matrix{ {A11, A12, A13, A14}, {A21, A22, A23, A24}, {A31, A32, A33, A34}, {A41, A42, A43, A44} }
p = A * transpose(matrix{ {x,y,z,w} })
X = p_(0,0)
Y = p_(1,0)
Z = p_(2,0)
W = p_(3,0)

f1 = (x^2-x*y+y^2+z^2-w^2) - (X^2-X*Y+Y^2+Z^2-W^2)
f2 = (3*x^2-10*x*y-5*x*w+4*y^2+z^2+6*y*w+2*w^2) - (3*X^2-10*X*Y-5*X*W+4*Y^2+Z^2+6*Y*W+2*W^2)
f3 = (x-y-w) - (X-Y-W)
f4 = z - Z
f5 = (2*x-y-w) - (2*X-Y-W)
f6 = (x+6*y+6*w) - (X+6*Y+6*W)

I1 = phi ideal((coefficients f1)_1)
I2 = phi ideal((coefficients f2)_1)
I3 = phi ideal((coefficients f3)_1)
I4 = phi ideal((coefficients f4)_1)
I5 = phi ideal((coefficients f5)_1)
I6 = phi ideal((coefficients f6)_1)
