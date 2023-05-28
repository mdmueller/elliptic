-- find the 2-torsion points of E by solving H \cap E = 2p + 2t
-- H = {ax + b(y + w) = 0}

loadPackage "NumericalAlgebraicGeometry"
loadPackage "MultiprojectiveVarieties"

R = QQ[q0,q1,q2]
q3=1

-* OLD NUMBERS
a=1
b=2
c=3
d=4
e=5
f=11
g=7
h=8*-

a=1
b=-1
c=0
d=0
e=3
f=-10
g=0
h=-5

I1 = ideal(a*q0^2+b*q0*q1+c*q0*q2+d*q0*q3+q1^2+q2^2-q3^2, e*q0^2+f*q0*q1+g*q0*q2+h*q0*q3+4*q1^2+q2^2+6*q1*q3+2*q3^2)
M = matrix{ {2*a*q0+b*q1+c*q2+d*q3, b*q0+2*q1, c*q0+2*q2, d*q0-2*q3},
     {2*e*q0+f*q1+g*q2+h*q3, f*q0+8*q1+6*q3, g*q0+2*q2, h*q0+6*q1+4*q3},
     {1,0,0,0},
     {0,1,0,1}}

I = I1+ideal(det(M))
--torsions = solveSystem(first entries gens I)
