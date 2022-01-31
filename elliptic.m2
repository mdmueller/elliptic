loadPackage "QuillenSuslin"
loadPackage "MultiprojectiveVarieties"

R=QQ[b,c,d,x1,y1,z1,w1,x2,y2,z2,w2,x3,y3,z3,w3,x4,y4,z4,w4,Degrees=>{3:{1,0,0,0,0},4:{0,1,0,0,0},4:{0,0,1,0,0},4:{0,0,0,1,0},4:{0,0,0,0,1}}]

-- p1,p2,p3,p4 in E
I = ideal(x1*x1+2*x1*y1+3*x1*z1+4*x1*w1+y1*y1+z1*z1-w1*w1,
    x2*x2+2*x2*y2+3*x2*z2+4*x2*w2+y2*y2+z2*z2-w2*w2,
    x3*x3+2*x3*y3+3*x3*z3+4*x3*w3+y3*y3+z3*z3-w3*w3,
    x4*x4+2*x4*y4+3*x4*z4+4*x4*w4+y4*y4+z4*z4-w4*w4,
    5*x1*x1+11*x1*y1+7*x1*z1+8*x1*w1+4*y1*y1+z1*z1+6*y1*w1+2*w1*w1,
    5*x2*x2+11*x2*y2+7*x2*z2+8*x2*w2+4*y2*y2+z2*z2+6*y2*w2+2*w2*w2,
    5*x3*x3+11*x3*y3+7*x3*z3+8*x3*w3+4*y3*y3+z3*z3+6*y3*w3+2*w3*w3,
    5*x4*x4+11*x4*y4+7*x4*z4+8*x4*w4+4*y4*y4+z4*z4+6*y4*w4+2*w4*w4)
