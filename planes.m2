loadPackage "QuillenSuslin"
loadPackage "MultiprojectiveVarieties"

R = QQ[alpha,beta,gamma,delta,q0,q1,q2,q3,r0,r1,r2,r3,Degrees=>{4:{1,0,0},4:{0,1,0},4:{0,0,1}}]

-- q,r in Q_1
I1 = ideal(q0*q0+2*q0*q1+3*q0*q2+4*q0*q3+q1*q1+q2*q2-q3*q3, r0*r0+2*r0*r1+3*r0*r2+4*r0*r3+r1*r1+r2*r2-r3*r3)

-- q,r in Q_2
I2 = ideal(5*q0*q0+11*q0*q1+7*q0*q2+8*q0*q3+4*q1*q1+q2*q2+6*q1*q3+2*q3*q3,
           5*r0*r0+11*r0*r1+7*r0*r2+8*r0*r3+4*r1*r1+r2*r2+6*r1*r3+2*r3*r3)

-- T_pQ_1 \cap T_pQ_2 \subset H_1 where p=q,r
M1 = matrix{ {alpha, beta, gamma, delta},
     	    {2*q0+2*q1+3*q2+4*q3, 2*q0+2*q1, 3*q0+2*q2, 4*q0-2*q3},
	    {10*q0+11*q1+7*q2+8*q3, 6*q0+8*q1+6*q3, 7*q0+2*q2, 8*q0+6*q1+4*q3} }
M2 = matrix{ {alpha, beta, gamma, delta},
     	    {2*r0+2*r1+3*r2+4*r3, 2*r0+2*r1, 3*r0+2*r2, 4*r0-2*r3},
	    {10*r0+11*r1+7*r2+8*r3, 6*r0+8*r1+6*r3, 7*r0+2*r2, 8*r0+6*r1+4*r3} }
I3 = maxMinors(M1)
I4 = maxMinors(M2)

I = I1+I2+I3+I4 -- ideal which should have 1- and 2-dimensional components

-- construct ideal of 2-dimensional component Z (q=r)
Z = ideal(q0*q0+2*q0*q1+3*q0*q2+4*q0*q3+q1*q1+q2*q2-q3*q3, 5*q0*q0+11*q0*q1+7*q0*q2+8*q0*q3+4*q1*q1+q2*q2+6*q1*q3+2*q3*q3) + maxMinors(matrix{ {q0, q1, q2, q3}, {r0, r1, r2, r3} }) + I3


-- find multidegree of I3 (codim 2 in P^3 x P^3 x P^3):
-- I3 in Chow is _alpha^2 + _alpha*beta + ...

-- coeff of beta^2
X0 = projectiveVariety(I3 + ideal(random({1,0,0},R), random({1,0,0},R), random({1,0,0},R), random({0,1,0},R), random({0,0,1},R), random({0,0,1},R), random({0,0,1},R)))
-- coeff of beta*gamma
X1 = projectiveVariety(I3 + ideal(random({1,0,0},R), random({1,0,0},R), random({1,0,0},R), random({0,1,0},R), random({0,1,0},R), random({0,0,1},R), random({0,0,1},R)))
-- coeff of alpha*beta
X2 = projectiveVariety(I3 + ideal(random({1,0,0},R), random({1,0,0},R), random({0,1,0},R), random({0,1,0},R), random({0,0,1},R), random({0,0,1},R), random({0,0,1},R)))
-- coeff of alpha^2
X3 = projectiveVariety(I3 + ideal(random({1,0,0},R), random({0,1,0},R), random({0,1,0},R), random({0,1,0},R), random({0,0,1},R), random({0,0,1},R), random({0,0,1},R)))
-- coeff of gamma^2
X4 = projectiveVariety(I3 + ideal(random({1,0,0},R), random({1,0,0},R), random({1,0,0},R), random({0,1,0},R), random({0,1,0},R), random({0,1,0},R), random({0,0,1},R)))
-- coeff of alpha*gamma
X5 = projectiveVariety(I3 + ideal(random({1,0,0},R), random({1,0,0},R), random({0,1,0},R), random({0,1,0},R), random({0,1,0},R), random({0,0,1},R), random({0,0,1},R)))

-- turns out: X1, X4, X5 are empty, degree(X0)=3, degree(X2)=2, degree(X3)=1
-- so [I3] = 3beta^2 + 2alpha*beta + alpha^2, [I4] = 3gamma^2 + 2alpha*gamma + alpha^2
