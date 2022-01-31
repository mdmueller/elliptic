loadPackage "MultiprojectiveVarieties"
loadPackage "NumericalAlgebraicGeometry"

-- the map t-x, E->E should correspond to a linear map P^5->P^5; let's find it
-- here E = {y^2*z = x^3 + x^2*z - x*z^2} \subset P^2 -> P^5 by v(x:y:z) = [z^2:xz:yz:x^2:xy:y^2]

-- take points p_1,...,p_6 in E, compute their images a_i in P^5 by v and images b_i by v\circ (t-x)
p1 = (21,-101,1)
p2 = (1,-41,1) -- 2p1
p3 = (-15,7,1) -- 3p1
p4 = (1,39,1) -- 4p1
p5 = (21,79,1) -- 5p1
p6 = (0,1,0) -- 6p1
p7 = (-9,49,1) -- square root of p1
p = {p1,p2,p3,p4,p5,p6,p7}
q = {p6,p5,p4,p3,p2,p1,p7} -- q_i = p1 - p_i
a = {}
b = {}
for i from 0 to 6 do (
    x = (p_i)_0;
    y = (p_i)_1;
    z = (p_i)_2;
    X = (q_i)_0;
    Y = (q_i)_1;
    Z = (q_i)_2;
    a = append(a, {z^2,x*z,y*z,x^2,x*y,y^2});
    b = append(b, {Z^2,X*Z,Y*Z,X^2,X*Y,Y^2});
    )
-- now we want an automorphism of P^5 sending a_i -> b_i
U = QQ[l0,l1,l2,l3,l4,l5]
f = map(U,ZZ)
A = f(transpose(matrix(a_{0..5})))
B = f(transpose(matrix(b_{0..5})))*diagonalMatrix{l0,l1,l2,l3,l4,l5}
N = B*inverse(A)
q = transpose(matrix{a_6})
sol = trim(minors(2, q|(N*q)))
l = {466560000, 6400, 46656, 10000, 72900, 1} -- values from sol
g = map(QQ,ZZ)
Asol = g(transpose(matrix(a_{0..5})))
Bsol = g(transpose(matrix(b_{0..5})))*diagonalMatrix(l)
M = Bsol*inverse(Asol) -- matrix defining the automorphism of P^5

R = QQ[symbol x,symbol y,symbol z]
I = ideal(y^2*z + x*y*z + y*z^2 - x^3 + x^2*z + 122*x*z^2 - 1721*z^3)
G = rationalMap {z^2,x*z,y*z,x^2,x*y,y^2}
E = G I -- ideal of E inside P^5
t = gens(ring(E)) -- t_0,...,t_5 are coordinates

q = (M*transpose(matrix{t}))_0
q2 = (transpose(M)*transpose(matrix{t}))_0
H = rationalMap {q_0,q_1,q_2,q_3,q_4,q_5} -- map E->E in P^5 sending a_i -> b_i
H2 = rationalMap {q2_0,q2_1,q2_2,q2_3,q2_4,q2_5} -- map S_1 -> S_2

F = rationalMap {z^2,2*x*z,2*y*z,x^2,2*x*y,y^2}
S1 = projectiveVariety(image F) -- surface S_1 inside (P^5)*
S2 = projectiveVariety(H2(image F)) -- surface S_2

-- project S1, S2 to the hyperplane H_0 = {t_4 = 0} (intersects E transversely)
P1 = multirationalMap {rationalMap {t_0,t_1,t_2,t_3,t_5} }
A1 = P1 S1 -- singular (why?), deg 4
A2 = P1 S2 -- nonsingular, deg 4
A = A1*A2 -- 16 distinct points

-- project to H_0, which is tangent to E at a_4 = 5t and transverse elsewhere (2,1,1,1,1):
-- H_0 = {62410/1171*t_0 + 12482/1171*t_1 + 6241/1171*t_2 + 12482/1171*t_3 - 8216/1171*t_4 + t_5 = 0}
P2 = multirationalMap {rationalMap {
	t_0 - 62410/1171*t_5,
	t_1 - 12482/1171*t_5,
	t_2 - 6241/1171*t_5,
	t_3 - 12482/1171*t_5,
	t_4 + 8216/1171*t_5}}
B1 = P2 S1 -- nonsingular, deg 4
B2 = P2 S2 -- nonsingular, deg 4
B = B1*B2 -- 8 regular points and 4 tangent points (8+4*2=16)

-- project to H_0 which is tangent to E at a_3 = 4t and a_4 = 5t and transverse at 2 other points (2,2,1,1):
-- H_0 = {46916/115*t_0 + 23458/115*t_1 - 2331/46*t_2 - 153/230*t_3 - 4*t_4 + t_5 = 0}
P3 = multirationalMap {rationalMap {
	t_0 - 46916/115*t_5,
	t_1 - 23458/115*t_5,
	t_2 + 2331/46*t_5,
	t_3 + 153/230*t_5,
	t_4 + 4*t_5}}
C1 = P3 S1
C2 = P3 S2
C = C1*C2

-- project to H_0 which is triply tangent to E at a_3 = 4t and transverse at p, t, -t (3,1,1,1):
-- H_0 = {861*t_0-83*t_1-21*t_2+2*t_3+t_4 = 0}
P4 = multirationalMap {rationalMap {
	t_0 - 861*t_4,
	t_1 + 83*t_4,
	t_2 + 21*t_4,
	t_3 - 2*t_4,
    	t_5}}
D1 = P4 S1
D2 = P4 S2
D = D1*D2

-- project to H_0 which is tangent to E at p, t, -t (2,2,2):
-- H_0 = {441*t_0 - 42*t_1 + t_3 = 0}
P5 = multirationalMap {rationalMap {
	t_0 - 441*t_3,
	t_1 + 42*t_3,
	t_2,
	t_4,
	t_5}}
F1 = P5 S1
F2 = P5 S2
F = F1*F2 -- TODO: something weird here, it's all nonreduced?


-- find tangent line to E at a_5=p, a_4=5t, a_3=4t, a_0=t
h0 = map(QQ,ring(E),a_0)
h3 = map(QQ,ring(E),a_3)
h4 = map(QQ,ring(E),a_4)
h5 = map(QQ,ring(E),a_5)
J0 = (h0(jacobian E))
J3 = (h3(jacobian E))
J4 = (h4(jacobian E))
J5 = (h5(jacobian E))

T = QQ[c0,c1,c2,c3,c4,c5] -- (P^5)*
I0 = minors(5, J0 | transpose(matrix{{c0,c1,c2,c3,c4,c5}})) -- hyperplane contains J0
I3 = minors(5, J3 | transpose(matrix{{c0,c1,c2,c3,c4,c5}})) -- hyperplane contains J3
I4 = minors(5, J4 | transpose(matrix{{c0,c1,c2,c3,c4,c5}})) -- hyperplane contains J4
I5 = minors(5, J5 | transpose(matrix{{c0,c1,c2,c3,c4,c5}})) -- hyperplane contains J5
I = I3+I4+I5

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
