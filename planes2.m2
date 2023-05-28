loadPackage "QuillenSuslin"
loadPackage "MultiprojectiveVarieties"

R = QQ[alpha,beta,gamma,delta,q0,q1,q2,q3,Degrees=>{4:{1,0},4:{0,1}}]

-- q in Q_1
I1 = ideal(q0*q0+2*q0*q1+3*q0*q2+4*q0*q3+q1*q1+q2*q2-q3*q3)

-- q in Q_2
I2 = ideal(5*q0*q0+11*q0*q1+7*q0*q2+8*q0*q3+4*q1*q1+q2*q2+6*q1*q3+2*q3*q3)

-- T_qQ_1 \cap T_qQ_2 \subset H
M = matrix{ {alpha, beta, gamma, delta},
     	    {2*q0+2*q1+3*q2+4*q3, 2*q0+2*q1, 3*q0+2*q2, 4*q0-2*q3},
	    {10*q0+11*q1+7*q2+8*q3, 6*q0+8*q1+6*q3, 7*q0+2*q2, 8*q0+6*q1+4*q3} }
I3 = maxMinors(M)

I = I1+I2+I3

-- q = p
M2 = matrix{ {0, -1, 0, 1}, {q0, q1, q2, q3} }
I4 = maxMinors(M2)
