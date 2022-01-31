loadPackage "MultiprojectiveVarieties"

R = QQ[q0,q1,q2,r0,r1,r2,a,Degrees=>{3:{1,0,0},3:{0,1,0},1:{0,0,1}}]
q3 = 1
r3 = 1

I1 = ideal(q0*q0+2*q0*q1+3*q0*q2+4*q0+q1*q1+q2*q2-1, r0*r0+2*r0*r1+3*r0*r2+4*r0+r1*r1+r2*r2-1)
I2 = ideal(5*q0*q0+11*q0*q1+7*q0*q2+8*q0+4*q1*q1+q2*q2+6*q1+2,
           5*r0*r0+11*r0*r1+7*r0*r2+8*r0+4*r1*r1+r2*r2+6*r1+2)
       
I3 = ideal(a*q0+(q1+1), a*r0+(r1+1))
p = (I1+I2+I3):intersect(ideal(r0,r2,r1+r3),ideal(q0,q2,q1+q3),ideal(q0-r0,q1-r1,q2-r2))
p = p:ideal(r0,r2,r1+r3)
p = p:ideal(q0,q2,q1+q3)

I4 = ideal(a*(q0-q1-q3)+q2, a*(r0-r1-r3)+r2)
t1 = (I1+I2+I4):intersect(ideal(r0,r2,r1+r3),ideal(q0,q2,q1+q3),ideal(q0-r0,q1-r1,q2-r2))
t1 = t1:ideal(r0,r2,r1+r3)
t1 = t1:ideal(q0,q2,q1+q3)
t1 = t1:ideal(q0-q3,q1,q2)
t1 = t1:ideal(r0-r3,r1,r2)

I5 = ideal(a*(2*q0-q1-q3)+q2, a*(2*r0-r1-r3)+r2)
t2 = (I1+I2+I5):intersect(ideal(r0,r2,r1+r3),ideal(q0,q2,q1+q3),ideal(q0-r0,q1-r1,q2-r2))
t2 = t2:ideal(r0,r2,r1+r3)
t2 = t2:ideal(q0,q2,q1+q3)
t2 = t2:ideal(q0-q1,q3-q1,q2)
t2 = t2:ideal(r0-r1,r3-r1,r2)

I6 = ideal(a*(q0+6*q1+6*q3)+q2, a*(r0+6*r1+6*r3)+r2)
t3 = (I1+I2+I6):intersect(ideal(r0,r2,r1+r3),ideal(q0,q2,q1+q3),ideal(q0-r0,q1-r1,q2-r2))
t3 = t3:ideal(r0,r2,r1+r3)
t3 = t3:ideal(q0,q2,q1+q3)
t3 = t3:ideal(q0+48/43*q3, q1+35/43*q3, q2)

q3 = symbol q3
r3 = symbol r3

A = QQ[q0,q1,q2,q3,r0,r1,r2,r3,a,b,Degrees=>{8:0,2:1}]
F = map(A, R)
P = homogenize(F p, b)
T1 = homogenize(F t1, b)
T2 = homogenize(F t2, b)
T3 = homogenize(F t3, b)
B = QQ[q0,q1,q2,q3,r0,r1,r2,r3,a,b,Degrees=>{4:1,6:0}]
F = map(B, A)
P = homogenize(F P, q3)
T1 = homogenize(F T1, q3)
T2 = homogenize(F T2, q3)
T3 = homogenize(F T3, q3)
C = QQ[q0,q1,q2,q3,r0,r1,r2,r3,a,b,Degrees=>{4:0,4:1,2:0}]
F = map(C, B)
P = saturate(homogenize(F P, r3), q3*r3*b)
T1 = saturate(homogenize(F T1, r3), q3*r3*b)
T2 = saturate(homogenize(F T2, r3), q3*r3*b)
T3 = saturate(homogenize(F T3, r3), q3*r3*b)

S = QQ[q0,q1,q2,q3,r0,r1,r2,r3,Degrees=>{4:{1,0},4:{0,1}}]
F = map(C, S)
Pfinal = preimage(F, P)
T1final = preimage(F, T1)
T2final = preimage(F, T2)
T3final = preimage(F, T3)

M = matrix{  {2*q0+2*q1+3*q2+4*q3, 2*q0+2*q1, 3*q0+2*q2, 4*q0-2*q3},
    {10*q0+11*q1+7*q2+8*q3, 6*q0+8*q1+6*q3, 7*q0+2*q2, 8*q0+6*q1+4*q3},
    {2*r0+2*r1+3*r2+4*r3, 2*r0+2*r1, 3*r0+2*r2, 4*r0-2*r3},
    {10*r0+11*r1+7*r2+8*r3, 6*r0+8*r1+6*r3, 7*r0+2*r2, 8*r0+6*r1+4*r3} }
-- forget about last column: D_1*M_1-D_2*M_2+D_3*M_3-D_4*M_4=0, so D_1*M_1-D_2*M_2 represents the plane containing these lines
D1 = det(submatrix'(M, {0}, {3}))
D2 = det(submatrix'(M, {1}, {3}))
-- H = {Ax+By+Cz+Dw=0}
A = M_(0,0)*D1 - M_(1,0)*D2
B = M_(0,1)*D1 - M_(1,1)*D2
C = M_(0,2)*D1 - M_(1,2)*D2
D = M_(0,3)*D1 - M_(1,3)*D2
-- map P^3xP^3->(P^3)* is [A:B:C:D], then composition with (P^3)*->P^2 is [B:C:D] <-> {x=By+Cz+Dw=0}

F = rationalMap {A,B,C,D}
