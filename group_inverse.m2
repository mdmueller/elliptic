loadPackage "NumericalAlgebraicGeometry"
loadPackage "MultiprojectiveVarieties"
loadPackage "NumericalImplicitization"

-- compute t-x for x in E

R = CC[q0,q1,q2,A]
q3 = 1

t0 = 1
t1 = 0
t2 = 0
t3 = 1

a=1
b=-1
c=0
d=0
e=3
f=-10
g=0
h=-5

I = ideal(a*q0^2+b*q0*q1+c*q0*q2+d*q0*q3+q1^2+q2^2-q3^2, e*q0^2+f*q0*q1+g*q0*q2+h*q0*q3+4*q1^2+q2^2+6*q1*q3+2*q3^2)

L = numericalSourceSample(I, 5)
Lout = {}

for i from 0 to 4 do (
    P = coordinates(L#i);
    x0 = P#0;
    x1 = P#1;
    x2 = P#2;
    x3 = 1;
    results = solveSystem({a*q0^2+b*q0*q1+c*q0*q2+d*q0*q3+q1^2+q2^2-q3^2, e*q0^2+f*q0*q1+g*q0*q2+h*q0*q3+4*q1^2+q2^2+6*q1*q3+2*q3^2,
	    A*(q0-q1-q3)+q2, A*(x0-x1-x3)+x2});
    Lout = append(Lout, results)
    )
    

