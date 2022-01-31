loadPackage "MultiprojectiveVarieties"
loadPackage "QuillenSuslin"

-- q in E, r in E
R = QQ[q0,q1,q2,q3,r0,r1,r2,r3,Degrees=>{4:{1,0},4:{0,1}}]

I1 = ideal(q0*q0+2*q0*q1+3*q0*q2+4*q0*q3+q1*q1+q2*q2-q3*q3, r0*r0+2*r0*r1+3*r0*r2+4*r0*r3+r1*r1+r2*r2-r3*r3)
I2 = ideal(5*q0*q0+11*q0*q1+7*q0*q2+8*q0*q3+4*q1*q1+q2*q2+6*q1*q3+2*q3*q3,
           5*r0*r0+11*r0*r1+7*r0*r2+8*r0*r3+4*r1*r1+r2*r2+6*r1*r3+2*r3*r3)

M1 = matrix{  {2*q0+2*q1+3*q2+4*q3, 2*q0+2*q1, 3*q0+2*q2, 4*q0-2*q3},
    {10*q0+11*q1+7*q2+8*q3, 6*q0+8*q1+6*q3, 7*q0+2*q2, 8*q0+6*q1+4*q3},
    {2*r0+2*r1+3*r2+4*r3, 2*r0+2*r1, 3*r0+2*r2, 4*r0-2*r3},
    {10*r0+11*r1+7*r2+8*r3, 6*r0+8*r1+6*r3, 7*r0+2*r2, 8*r0+6*r1+4*r3} }

I3 = ideal(det(M1))
