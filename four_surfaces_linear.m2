-- P^5 version of four_curves_linear.m2

loadPackage "MultiprojectiveVarieties"

S = QQ[x,y,z]
I = ideal(y^2*z - x*(x-6*z)*(x-7*z))
F = rationalMap {z^2,x*z,y*z,x^2,x*y,y^2}
f = toMap F
E = F I -- ideal of E inside P^5
J = transpose(jacobian E)

R = QQ[x1,y1,z1,x2,y2,z2,Degrees=>{3:{1,0},3:{0,1}}]
p1 = map(R,S,{x1,y1,z1})
p2 = map(R,S,{x2,y2,z2})

-- c = -(a+b)
w = (x2*z1-x1*z2)
m = (y2*z1-y1*z2)
x3 = w^2*(13*z1*z2 - x2*z1 - x1*z2) + z1*z2*m^2
y3 = m*w*(x3-x1*z2)+y1*z2*w^2
z3 = z1*z2*w^2
p3 = map(R,S,{x3,y3,z3})

N = J^{3..6} -- rows 3,4,5,6 should generally be linearly independent
A1 = p1(f(N))
A2 = p2(f(N))
A3 = p3(f(N))
A = A1 || A2 -- A is an 8x6 matrix with rank 6

B = submatrix'(A,{7},{}) -- remove last row from A; B is 7x6
D1 = det(submatrix'(B,{0},{}))
D2 = det(submatrix'(B,{1},{}))
D3 = det(submatrix'(B,{2},{}))
D4 = det(submatrix'(B,{3},{}))
D5 = det(submatrix'(B,{4},{}))
D6 = det(submatrix'(B,{5},{}))
D7 = det(submatrix'(B,{6},{}))
-- D1A1 - D2A2 + ... + D7A7 = 0

C = submatrix'(A,{0},{}) -- remove first row from A; C is 7x6
E1 = det(submatrix'(C,{0},{}))
E2 = det(submatrix'(C,{1},{}))
E3 = det(submatrix'(C,{2},{}))
E4 = det(submatrix'(C,{3},{}))
E5 = det(submatrix'(C,{4},{}))
E6 = det(submatrix'(C,{5},{}))
E7 = det(submatrix'(C,{6},{}))
-- E1A2 - E2A3 + ... + E7A8 = 0

M = matrix{ {D1*A^{0}-D2*A^{1}+D3*A^{2}-D4*A^{3}}, {E1*A^{1}-E2*A^{2}+E3*A^{3}} } || A3 -- M is 6x6, rank 5
-- remove the 1st column to get 6x5; should hopefully be rank 5 still
M2 = submatrix'(M,{},{0})-*
F1 = det(submatrix'(M2,{0},{}))
F2 = det(submatrix'(M2,{1},{}))
F3 = det(submatrix'(M2,{2},{}))
F4 = det(submatrix'(M2,{3},{}))
F5 = det(submatrix'(M2,{4},{}))
F6 = det(submatrix'(M2,{5},{}))
-- F1M1 - F2M2 + ... - F6M6 = 0
H = F1*M^{0} - F2*M^{1} *-

-- TODO: would it help to work out the group law inside P^5xP^5 (rather than P^2xP^2) directly?

-*
-- forget about last column: D_1*M_1-D_2*M_2+D_3*M_3-D_4*M_4+D_5*M_5-D_6*M_6=0, so D_1*M_1-D_2*M_2 represents the plane containing these lines
D1 = det(submatrix'(M, {0}, {5}))
D2 = det(submatrix'(M, {1}, {5}))
-- H = {Ax+By+Cz+Dw=0}
R0 = M_(0,0)*D1 - M_(1,0)*D2
R1 = M_(0,1)*D1 - M_(1,1)*D2
R2 = M_(0,2)*D1 - M_(1,2)*D2
R3 = M_(0,3)*D1 - M_(1,3)*D2
R4 = M_(0,4)*D1 - M_(1,4)*D2
R5 = M_(0,5)*D1 - M_(1,5)*D2
*-
