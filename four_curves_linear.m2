-- four_curves_affine.m2, but using linear automorphisms for x->t-x

needsPackage "MultiprojectiveVarieties"

-- x->-x
Aut1 = matrix{ {1, 0, 0, 0},
    {0, 1, 0, 0},
    {0, 0, -1, 0},
    {0, 0, 0, 1} }

-- x->t1-x
Aut2 = matrix{ {-2/7, -6/7, 0, 2/7},
    {-9/8, 1/4, 0, 1/4},
    {0, 0, 1, 0},
    {-9/56, -3/28, 0, 29/28} }

-- x->t2-x
Aut3 = matrix{ {20/13, -6/13, 0, -14/13},
    {-35/104, 67/52, 0, 35/52},
    {0, 0, 1, 0},
    {147/104, -63/52, 0, -95/52} }

-- x->t3-x
Aut4 = matrix{ {-23/91, 120/91, 0, 72/91},
    {19/13, -7/13, 0, -12/13},
    {0, 0, 1, 0},
    {-114/91, 120/91, 0, 163/91} }

R := QQ[q0,q1,q2,q3,r0,r1,r2,r3,Degrees=>{4:{1,0},4:{0,1}}]
I1 := ideal(q0^2 - q0*q1 + q1^2 + q2^2 - q3^2, 3*q0^2 - 10*q0*q1 - 5*q0*q3 + 4*q1^2 + q2^2 + 6*q1*q3 + 2*q3^2)
I2 := ideal(r0^2 - r0*r1 + r1^2 + r2^2 - r3^2, 3*r0^2 - 10*r0*r1 - 5*r0*r3 + 4*r1^2 + r2^2 + 6*r1*r3 + 2*r3^2)

q := matrix{ {q0}, {q1}, {q2}, {q3} }
r := matrix{ {r0}, {r1}, {r2}, {r3} }

-- ideals of the four components of B = {(x,y) in E^2: 2x+2y=0}
B1 := I1+minors(2, (Aut1*q) | r)
B2 := I1+minors(2, (Aut2*q) | r)
B3 := I1+minors(2, (Aut3*q) | r)
B4 := I1+minors(2, (Aut4*q) | r)

M := matrix{  {2*q0-q1, -q0+2*q1, 2*q2, -2*q3},
    {6*q0-10*q1-5*q3, -10*q0+8*q1+6*q3, 2*q2, -5*q0+6*q1+4*q3},
    {2*r0-r1, -r0+2*r1, 2*r2, -2*r3},
    {6*r0-10*r1-5*r3, -10*r0+8*r1+6*r3, 2*r2, -5*r0+6*r1+4*r3} }
-- forget about last column: D_1*M_1-D_2*M_2+D_3*M_3-D_4*M_4=0, so D_1*M_1-D_2*M_2 represents the plane containing these lines
D1 := det(submatrix'(M, {0}, {3}))
D2 := det(submatrix'(M, {1}, {3}))
-- H = {Ax+By+Cz+Dw=0}
R0 := M_(0,0)*D1 - M_(1,0)*D2
R1 := M_(0,1)*D1 - M_(1,1)*D2
R2 := M_(0,2)*D1 - M_(1,2)*D2
R3 := M_(0,3)*D1 - M_(1,3)*D2
-- map P^3xP^3->(P^3)* is [R0:R1:R2:R3], then composition with (P^3)*->P^2 is [R1:R2:R3] <-> {x=R1y+R2z+R3w=0}

F := rationalMap {R0,R1,R2,R3}
C1 = projectiveVariety(F B1)
C2 = projectiveVariety(F B2)
C3 = projectiveVariety(F B3)
C4 = projectiveVariety(F B4)
X22 = C1+C2+C3+C4
