-- embed elliptic curve in P^5 via O(6p)

loadPackage "MultiprojectiveVarieties"

R = QQ[x,y,z]
I = ideal(y^2*z-x*(x-z)*(x-2*z))

F = rationalMap {z^2,x*z,y*z,x^2,x*y,y^2}
J = F I
X = projectiveVariety(J)
