-- P^5 version of four_curves_linear.m2
-- first attempt, working entirely in P^5

loadPackage "MultiprojectiveVarieties"

-- x->-x
A1 = matrix{ {1, 0, 0, 0, 0, 0},
    {0, 1, 0, 0, 0, 0},
    {0, 0, -1, 0, 0, 0},
    {0, 0, 0, 1, 0, 0},
    {0, 0, 0, 0, -1, 0},
    {0, 0, 0, 0, 0, 1} }

-- x->t1-x
A2 = matrix{
    { 0, -1/42, 0, 13/1764, 0, 1/1764},
    { 0, 0, 0, 1/42, 0, 0},
    { 0, 0, 0, 0, 1/42, 0},
    { 0, 1, 0, 0, 0, 0},
    { 0, 0, 1, 0, 0, 0},
    {42, -13, 0, 1, 0, 0} }

-- x->t2-x
A3 = matrix{
    {-6, 11/6, 0, -5/36, 0, 1/36},
    {-42, 13, 0, -1, 0, 1/6},
    {0, 0, 1, 0, -1/6, 0},
    {-294, 91, 0, -7, 0, 1},
    {0, 0, 7, 0, -1, 0},
    {0, -7, 0, 1, 0, 0} }

-- x->t3-x
A4 = matrix{
    {-7, 15/7, 0, -8/49, 0, 1/49},
    {-42, 13, 0, -1, 0, 1/7},
    {0, 0, -1, 0, 1/7, 0},
    {-252, 78, 0, -6, 0, 1},
    {0, 0, -6, 0, 1, 0},
    {0, -6, 0, 1, 0, 0}}

R = QQ[q0,q1,q2,q3,q4,q5,r0,r1,r2,r3,r4,r5,s0,s1,s2,s3,s4,s5,Degrees=>{6:{1,0,0},6:{0,1,0},6:{0,0,1}}]
I1 = ideal(q_4^2-q_3*q_5,q_2*q_4-q_1*q_5,42*q_0*q_4-13*q_1*q_4+q_3*q_4-q_2*q_5,q_2*q_3-q_1*q_4,42*q_0*q_3-13*q_1*q_3+q_3^2-q_1*q_5,
    q_2^2-q_0*q_5,42*q_1*q_2-13*q_1*q_4+q_3*q_4-q_2*q_5,42*q_1^2-13*q_1*q_3+q_3^2-q_1*q_5,1764*q_0*q_1-127*q_1*q_3+13*q_3^2-42*q_0*q_5-13*q_1*q_5)
I2 = ideal(r_4^2-r_3*r_5,r_2*r_4-r_1*r_5,42*r_0*r_4-13*r_1*r_4+r_3*r_4-r_2*r_5,r_2*r_3-r_1*r_4,42*r_0*r_3-13*r_1*r_3+r_3^2-r_1*r_5,
    r_2^2-r_0*r_5,42*r_1*r_2-13*r_1*r_4+r_3*r_4-r_2*r_5,42*r_1^2-13*r_1*r_3+r_3^2-r_1*r_5,1764*r_0*r_1-127*r_1*r_3+13*r_3^2-42*r_0*r_5-13*r_1*r_5)
I3 = ideal(s_4^2-s_3*s_5,s_2*s_4-s_1*s_5,42*s_0*s_4-13*s_1*s_4+s_3*s_4-s_2*s_5,s_2*s_3-s_1*s_4,42*s_0*s_3-13*s_1*s_3+s_3^2-s_1*s_5,
    s_2^2-s_0*s_5,42*s_1*s_2-13*s_1*s_4+s_3*s_4-s_2*s_5,42*s_1^2-13*s_1*s_3+s_3^2-s_1*s_5,1764*s_0*s_1-127*s_1*s_3+13*s_3^2-42*s_0*s_5-13*s_1*s_5)

q = matrix{ {q0}, {q1}, {q2}, {q3}, {q4}, {q5} }
r = matrix{ {r0}, {r1}, {r2}, {r3}, {r4}, {r5} }
s = matrix{ {s0}, {s1}, {s2}, {s3}, {s4}, {s5} }

-- ideals of the four components of B = {(x,y) in E^2: 2x+2y+2z=0}
B1 = I1+minors(2, (A1*q) | r)
B2 = I1+minors(2, (A2*q) | r)
B3 = I1+minors(2, (A3*q) | r)
B4 = I1+minors(2, (A4*q) | r)

M = matrix{  {2*q0-q1, -q0+2*q1, 2*q2, -2*q3},
    {6*q0-10*q1-5*q3, -10*q0+8*q1+6*q3, 2*q2, -5*q0+6*q1+4*q3},
    {2*r0-r1, -r0+2*r1, 2*r2, -2*r3},
    {6*r0-10*r1-5*r3, -10*r0+8*r1+6*r3, 2*r2, -5*r0+6*r1+4*r3} }
-- forget about last column: D_1*M_1-D_2*M_2+D_3*M_3-D_4*M_4=0, so D_1*M_1-D_2*M_2 represents the plane containing these lines
D1 = det(submatrix'(M, {0}, {3}))
D2 = det(submatrix'(M, {1}, {3}))
-- H = {Ax+By+Cz+Dw=0}
R0 = M_(0,0)*D1 - M_(1,0)*D2
R1 = M_(0,1)*D1 - M_(1,1)*D2
R2 = M_(0,2)*D1 - M_(1,2)*D2
R3 = M_(0,3)*D1 - M_(1,3)*D2
-- map P^3xP^3->(P^3)* is [R0:R1:R2:R3], then composition with (P^3)*->P^2 is [R1:R2:R3] <-> {x=R1y+R2z+R3w=0}

F = rationalMap {R0,R1,R2,R3}
C1 = projectiveVariety(F B1)
C2 = projectiveVariety(F B2)
C3 = projectiveVariety(F B3)
C4 = projectiveVariety(F B4)
C = C1+C2+C3+C4

S = target F
t0 = (gens S)_0
t1 = (gens S)_1
t2 = (gens S)_2
t3 = (gens S)_3

G = multirationalMap {rationalMap {t1,t2,t3} }
D1 = G C1
D2 = G C2
D3 = G C3
D4 = G C4
D = D1*D2+D1*D3+D1*D4+D2*D3+D2*D4+D3*D4

-- use the plane H_0 = {-2/3*x+5*y+w=0} (meets E at 2 distinct tangent points) to define another map (P^3)*->P^2, intersecting with H_0
H = multirationalMap {rationalMap {t0+2/3*t3, t1-5*t3, t2} }
M1 = H C1
M2 = H C2
M3 = H C3
M4 = H C4
M = M1*M2+M1*M3+M1*M4+M2*M3+M2*M4+M3*M4

-- use the plane H_0 = {27x-6y+8sqrt(10)*z-38w=0} (meets E in (2,1,1)) to define another map (P^3)*->P^2, intersecting with H_0
H = multirationalMap {rationalMap {t0+27/38*t3, t1-6/38*t3, t2+8*w/38*t3}}
N1 = H C1
N2 = H C2
N3 = H C3
N4 = H C4
N = N1*N2+N1*N3+N1*N4+N2*N3+N2*N4+N3*N4

-- use the plane H_0 = {x+y+z+w=0} (meets E transversely) to define another map (P^3)*->P^2, intersecting with H_0
H = multirationalMap {rationalMap {t0-t3, t1-t3, t2-t3}}
O1 = H C1
O2 = H C2
O3 = H C3
O4 = H C4
O = O1*O2+O1*O3+O1*O4+O2*O3+O2*O4+O3*O4

-- use the plane H_0 = {-6451/7264*x+45/454*y+30*sqrt(10)/227*z+w=0} (meets E at a,b, one with order 3) to define another map (P^3)*->P^2, intersecting with H_0
H = multirationalMap {rationalMap {t0+6451/7264*t3, t1-45/454*t3, t2-30*w/(227)*t3} }
P1 = H C1
P2 = H C2
P3 = H C3
P4 = H C4
P = P1*P2+P1*P3+P1*P4+P2*P3+P2*P4+P3*P4
