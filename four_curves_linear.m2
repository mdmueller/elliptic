-- four_curves_affine.m2, but using linear automorphisms for x->t-x

loadPackage "MultiprojectiveVarieties"

-- x->-x
A1 = matrix{ {1, 0, 0, 0},
    {0, 1, 0, 0},
    {0, 0, -1, 0},
    {0, 0, 0, 1} }

-- x->t1-x
A2 = matrix{ {-2/7, -6/7, 0, 2/7},
    {-9/8, 1/4, 0, 1/4},
    {0, 0, 1, 0},
    {-9/56, -3/28, 0, 29/28} }

-- x->t2-x
A3 = matrix{ {20/13, -6/13, 0, -14/13},
    {-35/104, 67/52, 0, 35/52},
    {0, 0, 1, 0},
    {147/104, -63/52, 0, -95/52} }

-- x->t3-x
A4 = matrix{ {-23/91, 120/91, 0, 72/91},
    {19/13, -7/13, 0, -12/13},
    {0, 0, 1, 0},
    {-114/91, 120/91, 0, 163/91} }

k = toField(QQ[w]/(w^2-10))
--k = QQ
R = k[q0,q1,q2,q3,r0,r1,r2,r3,Degrees=>{4:{1,0},4:{0,1}}]
I1 = ideal(q0^2 - q0*q1 + q1^2 + q2^2 - q3^2, 3*q0^2 - 10*q0*q1 - 5*q0*q3 + 4*q1^2 + q2^2 + 6*q1*q3 + 2*q3^2)
I2 = ideal(r0^2 - r0*r1 + r1^2 + r2^2 - r3^2, 3*r0^2 - 10*r0*r1 - 5*r0*r3 + 4*r1^2 + r2^2 + 6*r1*r3 + 2*r3^2)

q = matrix{ {q0}, {q1}, {q2}, {q3} }
r = matrix{ {r0}, {r1}, {r2}, {r3} }

-- ideals of the four components of B = {(x,y) in E^2: 2x+2y=0}
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
E = ideal(t0^2 - t0*t1 + t1^2 + t2^2 - t3^2, 3*t0^2 - 10*t0*t1 - 5*t0*t3 + 4*t1^2 + t2^2 + 6*t1*t3 + 2*t3^2)

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

-- understanding nonreduced pts
-- A = {(l, x1, x2): H_0, x_1,x_2 in l} where H_0 = {27x-6y+8sqrt(10)*z-38w=0}
-- x1 = [q0:q1:q2:q3] in l
T = k[q0,q1,q2,q3,r0,r1,r2,r3,l0,l1,l2, Degrees=>{4:{1,0,0},4:{0,1,0},3:{0,0,1}}]
M = matrix{{27,-6,8*w,-38},{l0,0,l1,l2}} -- l
A1 = minors(3,M||matrix{{q0,q1,q2,q3}}) -- x1 in l
A2 = minors(3,M||matrix{{r0,r1,r2,r3}}) -- x2 in l
A = A1+A2

h1 = map(T,S,{q0,q1,q2,q3})
h2 = map(T,S,{r0,r1,r2,r3})
-- B = A \cap {x1 in C_1, x2 in C_2} (should be 0-dimensional)
B1 = h1(ideal(C1))
B2 = h2(ideal(C2))
B = A+B1+B2
Y = projectiveVariety(B) -- 0 dim, degree 4 (as expected)

-- D = B \cap {l=Va}
Va = minors(2,matrix{{l0,l1,l2},{-76,0,64}})
D = B + Va
Z = projectiveVariety(D) -- 0 dim, degree 1 (?)

-- E = B \cap {l\subset Ua}, Ua={H: a\in H}
Ua = ideal(16*l0+4*w*l1+19*l2)
E = B + Ua
W = projectiveVariety(E)
