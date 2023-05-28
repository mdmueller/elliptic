-- like six_torsion.m2, but use an elliptic curve y^2 = x(x-6)(x-7) with rational 2-torsion points

loadPackage "MultiprojectiveVarieties"
loadPackage "NumericalAlgebraicGeometry"

-- takes projective 0-dim ideal I, returns lengths of points
getLengths = I -> (
    K := primaryDecomposition(I);
    L := {};
    for j from 1 to #K do (
	J = K_(j-1);
	D1 = degree(projectiveVariety(J));
	D2 = degree(projectiveVariety(radical(J)));
	L = append(L, {D1,D2});
	);
    L)

-- for t a 2-torsion, the map t-x, E->E should correspond to a linear map P^5->P^5; let's find it
-- here E = {y^2*z = x*(x-6*z)*(x-7*z)} \subset P^2 -> P^5 by v(x:y:z) = [z^2:xz:yz:x^2:xy:y^2]
-- take some points in E, compute their images a_i in P^5 by v and images b_i by v\circ (t-x)

p = (0,1,0) -- t0
t1 = (0,0,1)
t2 = (6,0,1)
t3 = (1,0,1/7) -- hacky: it's (7,0,1), but want a rational # somewhere so types will be QQ not ZZ
q1 = (3,6,1) 
r1 = (3,-6,1)
q2 = (8,4,1)
r2 = (8,-4,1)
q3 = (14,28,1)
r3 = (14,-28,1)
q4 = (21/4,21/8,1)
r4 = (21/4,-21/8,1)
ellipt = {p,t1,t2,t3,q1,q2,q3}

-- images_{i,j} = t_i - ellipt_j    
images = {{p,t1,t2,t3,r1,r2,r3}, {t1,p,t3,t2,q3,q4,q1}, {t2,t3,p,t1,r2,r1,r4}, {t3,t2,t1,p,q4,q3,q2}}
matrices = {}
solns = {}
for i from 0 to 3 do (
    a = {};
    b = {};
    for j from 0 to 6 do (
    	x = (ellipt_j)_0;
    	y = (ellipt_j)_1;
    	z = (ellipt_j)_2;
    	X = ((images_i)_j)_0;
    	Y = ((images_i)_j)_1;
    	Z = ((images_i)_j)_2;
    	a = append(a, {z^2,x*z,y*z,x^2,x*y,y^2});
    	b = append(b, {Z^2,X*Z,Y*Z,X^2,X*Y,Y^2});
    	);
    -- now we want an automorphism of P^5 sending a_i -> b_i
    U = QQ[l0,l1,l2,l3,l4,l5];
    f = map(U,QQ);
    A = f(transpose(matrix(a_{0..5})));
    B = f(transpose(matrix(b_{0..5})))*diagonalMatrix{l0,l1,l2,l3,l4,l5};
    N = B*inverse(A);
    q = transpose(matrix{a_6});
    q2 = transpose(matrix{b_6});
    sol = trim(minors(2, q2|(N*q)));
    solns = append(solns,sol);
    -- values from sol
    l = {{1,1,1,1,1,1}, {1/512,9261/64,1323/64,7/512,27/512,1}, {1/8,-1323,-27,1/392,-27/8,1}, {49,-343,-1,7,-64,1}}_i;
    Asol = transpose(matrix(a_{0..5}));
    Bsol = transpose(matrix(b_{0..5}))*diagonalMatrix(l);
    M = Bsol*inverse(Asol); -- matrix defining the automorphism of P^5
    matrices = append(matrices, M);
)

R = QQ[symbol x,symbol y,symbol z]
I = ideal(y^2*z - x*(x-6*z)*(x-7*z))
G = rationalMap {z^2,x*z,y*z,x^2,x*y,y^2}
E = G I -- ideal of E inside P^5
t = gens(ring(E)) -- t_0,...,t_5 are coordinates

w = (transpose(matrices_0)*transpose(matrix{t}))_0
H0 = rationalMap {w_0,w_1,w_2,w_3,w_4,w_5}
w = (transpose(matrices_1)*transpose(matrix{t}))_0
H1 = rationalMap {w_0,w_1,w_2,w_3,w_4,w_5}
w = (transpose(matrices_2)*transpose(matrix{t}))_0
H2 = rationalMap {w_0,w_1,w_2,w_3,w_4,w_5}
w = (transpose(matrices_3)*transpose(matrix{t}))_0
H3 = rationalMap {w_0,w_1,w_2,w_3,w_4,w_5}

F = rationalMap {z^2,2*x*z,2*y*z,x^2,2*x*y,y^2}
S1 = projectiveVariety(image F) -- surface S_1 inside (P^5)*
S2 = projectiveVariety(H1(image F)) -- surface S_2
S3 = projectiveVariety(H2(image F)) -- surface S_3
S4 = projectiveVariety(H3(image F)) -- surface S_4

-- project S1, S2 to the hyperplane H_0 = {t_0-t_1+t_2-t_3+t_4-t_5 = 0} (intersects E transversely)
P1 = multirationalMap {rationalMap {t_0+t_5,t_1-t_5,t_2+t_5,t_3-t_5,t_4+t_5} }
A1 = P1 S1 -- nonsingular, deg 4
A2 = P1 S2 -- nonsingular, deg 4
A = A1*A2 -- 16 distinct points

-- project to H_0, which is tangent to E at a_4 = q1 and transverse elsewhere (2,1,1,1,1):
-- H_0 = {5*t_0 + 3*t_1 + t_2 + 442/99*t_3 - 529/99*t_4 + t_5 = 0}
P2 = multirationalMap {rationalMap {
	t_0 - 5*t_5,
	t_1 - 3*t_5,
	t_2 - t_5,
	t_3 - 442/99*t_5,
	t_4 + 529/99*t_5}}
B1 = P2 S1 -- nonsingular, deg 4
B2 = P2 S2 -- nonsingular, deg 4
B = B1*B2 -- 8 regular points and 4 tangent points (8+4*2=16)

-- project to H_0 which is tangent to E at a_4 = q1 and a_5 = q2 and transverse at 2 other points (2,2,1,1):
-- H_0 = {159984/37655*t_0 + 159984/37655*t_1 - 65208/7531*t_2 - 13924/37655*t_3 + 956/7531*t_4 + t_5 = 0}
P3 = multirationalMap {rationalMap {
	t_0 - 159984/37655*t_5,
	t_1 - 159984/37655*t_5,
	t_2 + 65208/7531*t_5,
	t_3 + 13924/37655*t_5,
	t_4 - 956/7531*t_5}}
C1 = P3 S1 -- nonsingular, deg 4
C2 = P3 S2 -- nonsingular, deg 4
C = C1*C2 -- 4 regular points, 4 points of order 2, 1 point of order 4 (4*1 + 4*2 + 1*4 = 16)

-- project to H_0 which is triply tangent to E at a_4 = q1 and transverse at a_5 = q2, a_6 = q3, and another point (1029/5476, -1105293/405224) (3,1,1,1):
-- H_0 = {-882/37*t_0 + 684/37*t_1 - 327/74*t_2 - 38/37*t_3 - 131/74*t_4 + t_5 = 0}
P4 = multirationalMap {rationalMap {
    t_0 + 882/37*t_5,
    t_1 - 684/37*t_5,
    t_2 + 327/74*t_5,
    t_3 + 38/37*t_5,
    t_4 + 131/74*t_5}}
D1 = P4 S1 -- nonsingular, deg 4
D2 = P4 S2 -- nonsingular, deg 4
D = D1*D2 -- 4 regular points, 4 points of order 3 (4*1 + 4*3 = 16)

-- project to H_0 which is tangent to E at a_4 = q1, a_5 = q2, and the point (54/25,792/125,1) (2,2,2):
-- H_0 = {1296/25*t_0 - 144/25*t_1 - 72/5*t_2 + 4/25*t_3 + 4/5*t_4 + t_5 = 0}
P5 = multirationalMap {rationalMap {
	t_0 - 1296/25*t_5,
	t_1 + 144/25*t_5,
	t_2 + 72/5*t_5,
	t_3 - 4/25*t_5,
	t_4 - 4/5*t_5}}
G1 = P5 S1 -- nonsingular, deg 3 (?)
G2 = P5 S2 -- nonsingular, deg 4
G3 = P5 S3
G4 = P5 S4
G = G1*G2 -- TODO: something weird here, it's all nonreduced?

-- project to H_0 which is tritangent to E at a_4 = q1, tangent at a_1 = t1, and transverse at the point (31827/5329,160062/389017,1) (3,2,1):
-- H_0 = {-3993/73*t_1+655/73*t_3+192/73*t_4+t_5 = 0}
P6 = multirationalMap {rationalMap {
	t_0,
	t_1 + 3993/73*t_5,
	t_2,
	t_3 - 655/73*t_5,
	t_4 - 192/73*t_5}}
W1 = P6 S1 -- nonsingular, deg 4
W2 = P6 S2 -- nonsingular, deg 4
W3 = P6 S3
W4 = P6 S4
W = W1*W2

-- project to H_0 which is 2-tangent to E at one point and 4-tangent at another (4,2):
-- H_0 = {-6*t_1 + t_3 + t_5 = 0}
P7 = multirationalMap {rationalMap {
	t_0,
	t_1 + 6*t_5,
	t_2,
	t_3 - t_5,
	t_4}}
Y1 = P7 S1 -- nonsingular, deg 3 (?)
Y2 = P7 S2 -- nonsingular, deg 4
Y3 = P7 S3
Y4 = P7 S4
Y = Y1*Y2

-- project to H_0 which is 2-tangent to E at one point and 4-tangent at another (4,2):
-- H_0 = {t_0 = 0}
P8 = multirationalMap {rationalMap {t_1,t_2,t_3,t_4,t_5}}
Z1 = P8 S1
Z2 = P8 S2
Z3 = P8 S3
Z4 = P8 S4
Z = Z1*Z2

-- project to H_0 which is 4-tangent to E at a_4=q1 and transverse elsewhere (4,1,1):
-- H_0 = {327129/5042*t_0 - 39292/2521*t_1 - 39292/2521*t_2 + 17197/15126*t_3 + 4084/2521*t_4 + t_5 = 0}
P9 = multirationalMap {rationalMap {
	t_0 - 327129/5042*t_5,
	t_1 + 39292/2521*t_5,
	t_2 + 39292/2521*t_5,
	t_3 - 17197/15126*t_5,
	t_4 - 4084/2521*t_5}}
Q1 = P9 S1
Q2 = P9 S2
Q3 = P9 S3
Q4 = P9 S4
Q = Q1*Q2


-- find tangent lines to E at our points
x = 54/25
y = 792/125
z = 1
P = {z^2,x*z,y*z,x^2,x*y,y^2}
h0 = map(QQ,ring(E),a_0)
h1 = map(QQ,ring(E),a_1)
h2 = map(QQ,ring(E),a_2)
h3 = map(QQ,ring(E),a_3)
h4 = map(QQ,ring(E),a_4)
h5 = map(QQ,ring(E),a_5)
hP = map(QQ,ring(E),P)
J0 = h0(jacobian E)
J1 = h1(jacobian E)
J2 = h2(jacobian E)
J3 = h3(jacobian E)
J4 = h4(jacobian E)
J5 = h5(jacobian E)
JP = hP(jacobian E)

T = ring(ideal(S1)) -- (P^5)*
tr = transpose(matrix{{t_0,t_1,t_2,t_3,t_4,t_5}})
I0 = minors(5, J0 | tr) -- hyperplane contains J0
I1 = minors(5, J1 | tr) -- hyperplane contains J1
I2 = minors(5, J2 | tr) -- hyperplane contains J2
I3 = minors(5, J3 | tr) -- hyperplane contains J3
I4 = minors(5, J4 | tr) -- hyperplane contains J4
I5 = minors(5, J5 | tr) -- hyperplane contains J5
IP = minors(5, JP | tr) -- hyperplane contains JP

-- here a = a_4 = [1:3:6:9:18:36], b = a_5 = [1:8:4:64:32:16]
Ua = ideal(t_0+3*t_1+6*t_2+9*t_3+18*t_4+36*t_5) -- hyperplanes containing a
Va = trim I4 -- hyperplanes tangent to E at a
Wa = Va + ideal(-73/192*t_2+t_3-121/64*t_4-4*t_5) -- hyperplanes with contact order >=3 at a
Ub = ideal(t_0+8*t_1+4*t_2+64*t_3+32*t_4+16*t_5)
U = QQ[x0,x1,x2,x3,x4] -- homogeneous ring for Ua = P^4
phi = map(U, T, {x0,x1,x2,x3,x4,-1/36*(x0+3*x1+6*x2+9*x3+18*x4)}) -- inclusion Ua -> P^5
C1 = phi(ideal(S1)) -- S_1 \cap Ua (nonreduced curve)
C2 = phi(ideal(S2)) -- S_2 \cap Ua (nonreduced curve)

-- H_0 = {5*t_0 + 3*t_1 + t_2 + 442/99*t_3 - 529/99*t_4 + t_5}, or [5:3:1:442/99:-529/99] in Ua (2,1,1,1,1)
h = {1, 3/5, 1/5, (442/99)/5, (-529/99)/5, 1/5} -- easiest if h_0 = 1
Pi = rationalMap{ x1-x0*h_1, x2-x0*h_2, x3-x0*h_3, x4-x0*h_4 } -- quotient map pi : Ua-{H_0} -> P^3
C1tilde = projectiveVariety(Pi(C1))
C2tilde = projectiveVariety(Pi(C2))
-- C1tilde \cap C2tilde = degree 8

-- now with H_0 = {-882/37*t_0 + 684/37*t_1 - 327/74*t_2 - 38/37*t_3 - 131/74*t_4 + t_5 = 0} (intersects E in (3,1,1,1))
g = -882/37
h = {1, (684/37)/g, (-327/74)/g, (-38/37)/g, (-131/74)/g, 1/g}
Pi = rationalMap{ x1-x0*h_1, x2-x0*h_2, x3-x0*h_3, x4-x0*h_4 } -- quotient map pi : Ua-{H_0} -> P^3
C1tilde = projectiveVariety(Pi(C1))
C2tilde = projectiveVariety(Pi(C2))
-- C1tilde \cap C2tilde = degree 8? Shouldn't it be 12? ...

-- now with H_0 = {159984/37655*t_0 + 159984/37655*t_1 - 65208/7531*t_2 - 13924/37655*t_3 + 956/7531*t_4 + t_5 = 0} (2,2,1,1)
g = 159984/37655
h = {1, 1, (-65208/7531)/g, (-13924/37655)/g, (956/7531)/g, 1/g}
Pi = rationalMap{ x1-x0*h_1, x2-x0*h_2, x3-x0*h_3, x4-x0*h_4 } -- quotient map pi : Ua-{H_0} -> P^3
C1tilde = projectiveVariety(Pi(C1))
C2tilde = projectiveVariety(Pi(C2))

-- now with H_0 = {327129/5042*t_0 - 39292/2521*t_1 - 39292/2521*t_2 + 17197/15126*t_3 + 4084/2521*t_4 + t_5 = 0} (4,1,1)
g = 327129/5042
h = {1, (-39292/2521)/g, (-39292/2521)/g, (17197/15126)/g, (4084/2521)/g, 1/g}
Pi = rationalMap{ x1-x0*h_1, x2-x0*h_2, x3-x0*h_3, x4-x0*h_4 } -- quotient map pi : Ua-{H_0} -> P^3
C1tilde = projectiveVariety(Pi(C1))
C2tilde = projectiveVariety(Pi(C2))

-- now with H_0 = {-2429811/135577*t_0 - 13442199/135577*t_1 + 1112256/135577*t_2 + 2180884/135577*t_3 + 643008/135577*t_4 + t_5 = 0} (5,1)
g = -2429811/135577
h = {1, (-13442199/135577)/g, (1112256/135577)/g, (2180884/135577)/g, (643008/135577)/g, 1/g}
Pi = rationalMap{ x1-x0*h_1, x2-x0*h_2, x3-x0*h_3, x4-x0*h_4 } -- quotient map pi : Ua-{H_0} -> P^3
C1tilde = projectiveVariety(Pi(C1))
C2tilde = projectiveVariety(Pi(C2))

-- now let's try working in all of (P^5)* rather than just Ua
Pi = rationalMap{ t_1-t_0*h_1, t_2-t_0*h_2, t_3-t_0*h_3, t_4-t_0*h_4, t_5-t_0*h_5 }
S1tilde = projectiveVariety(Pi(ideal(S1)))
S2tilde = projectiveVariety(Pi(ideal(S2)))

-*
-- H_0 = {62410/1171*t_0 + 12482/1171*t_1 + 6241/1171*t_2 + 12482/1171*t_3 - 8216/1171*t_4 + t_5 = 0}
-- find pairs (H_1,H_2) where H_1,H_2 both contain a and H_0,H_1,H_2 collinear
S = QQ[r0,r1,r2,r3,r4,r5,s0,s1,s2,s3,s4,s5,Degrees=>{6:{1,0},6:{0,1}}]
f=map(S,ring(ideal(S1)),{r0,r1,r2,r3,r4,r5})
T1 = f(ideal(S1))
g=map(S,ring(ideal(S2)),{s0,s1,s2,s3,s4,s5})
T2 = g(ideal(S2))
-- H_1 contains a
T3 = ideal(r0+21*r1+79*r2+441*r3+1659*r4+6241*r5)
-- H_2 contains a
T4 = ideal(s0+21*s1+79*s2+441*s3+1659*s4+6241*s5)
M = matrix{ {r0,r1,r2,r3,r4,r5},
    {s0,s1,s2,s3,s4,s5},
--    {62410/1171, 12482/1171, 6241/1171, 12482/1171, -8216/1171, 1} }
    {46916/115, 23458/115, -2331/46, -153/230, -4, 1} }
I = T1+T2+T3+T4+minors(3,M)
X = projectiveVariety(I)
*-
