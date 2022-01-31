loadPackage "MultiprojectiveVarieties"
loadPackage "QuillenSuslin"

-- q in E, r in E
R = QQ[q0,q1,q2,r0,r1,r2]

I1 = ideal(q0*q0+2*q0*q1+3*q0*q2+4*q0+q1*q1+q2*q2-1, r0*r0+2*r0*r1+3*r0*r2+4*r0+r1*r1+r2*r2-1)
I2 = ideal(5*q0*q0+11*q0*q1+7*q0*q2+8*q0+4*q1*q1+q2*q2+6*q1+2,
           5*r0*r0+11*r0*r1+7*r0*r2+8*r0+4*r1*r1+r2*r2+6*r1+2)

M1 = matrix{  {2*q0+2*q1+3*q2+4, 2*q0+2*q1, 3*q0+2*q2, 4*q0-2},
    {10*q0+11*q1+7*q2+8, 6*q0+8*q1+6, 7*q0+2*q2, 8*q0+6*q1+4},
    {2*r0+2*r1+3*r2+4, 2*r0+2*r1, 3*r0+2*r2, 4*r0-2},
    {10*r0+11*r1+7*r2+8, 6*r0+8*r1+6, 7*r0+2*r2, 8*r0+6*r1+4} }

I3 = ideal(det(M1))
