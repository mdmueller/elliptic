loadPackage "MultiprojectiveVarieties"

-- the map t-x, E->E should correspond to a linear map P^5->P^5; let's find it
-- here E = {y^2*z = x^3 + x^2*z - x*z^2} \subset P^2 -> P^5 by v(x:y:z) = [z^2:xz:yz:x^2:xy:y^2]

-- take points p_1,...,p_6 in E, compute their images a_i in P^5 by v and images b_i by v\circ (t-x)
p1 = (-1,1,1)
p2 = (1,-1,1) -- 2p1
p3 = (0,0,1) -- 3p1
p4 = (1,1,1) -- 4p1
p5 = (-1,-1,1) -- 5p1
p6 = (0,1,0) -- 6p1
p = {p1,p2,p3,p4,p5,p6}
q = {p6,p5,p4,p3,p2,p1} -- q_i = p1 - p_i
a = {}
b = {}
for i from 0 to 5 do (
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
-- take M = BA^{-1}
f = map(QQ,ZZ)
A = f(transpose(matrix a))
B = f(transpose(matrix b))
M = B*inverse(A)

R = QQ[symbol x,symbol y,symbol z]
I = ideal(y^2*z - x^3 - x^2*z + x*z^2)
G = rationalMap {z^2,x*z,y*z,x^2,x*y,y^2}
E = G I -- ideal of E inside P^5
t = gens(ring(E)) -- t_0,...,t_5 are coordinates
F = rationalMap {z^2,2*x*z,2*y*z,x^2,2*x*y,y^2}
S1 = image F -- ideal of surface S_1 inside (P^5)*

-- find tangent line to E at a_5=p, a_4=5t, a_3=4t
h3 = map(QQ,ring(E),a_3)
h4 = map(QQ,ring(E),a_4)
h5 = map(QQ,ring(E),a_5)
J3 = (h3(jacobian E))
J4 = (h4(jacobian E))
J5 = (h5(jacobian E))

T = QQ[c0,c1,c2,c3,c4,c5] -- (P^5)*
I3 = minors(5, J3 | transpose(matrix{{c0,c1,c2,c3,c4,c5}})) -- hyperplane contains J3
I4 = minors(5, J4 | transpose(matrix{{c0,c1,c2,c3,c4,c5}})) -- hyperplane contains J4
I5 = minors(5, J5 | transpose(matrix{{c0,c1,c2,c3,c4,c5}})) -- hyperplane contains J5
I = I3+I4+I5

