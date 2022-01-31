loadPackage "QuillenSuslin"

R = QQ[a1, a2, a3, a4, a5, b1, b2, y, z] -- P^5 x P^2 x P^2 (conic, line, point)
S = QQ[A1, A2, A3, A4, A5, B1, B2]

M = matrix{
    {2+a1*y+a3*z, a1+2*a2*y+a4*z, a3+a4*y+2*a5*z},
    {1, b1, b2}}
I = ideal(1+a1*y+a2*y*y+a3*z+a4*y*z+a5*z*z, 1+b1*y+b2*z) + maxMinors(M)

F = map(R/I, S, {a1, a2, a3, a4, a5, b1, b2}) -- S inside R/I
