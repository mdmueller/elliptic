import numpy as np
from sympy import *
from random import random
from cmath import sqrt
from scipy.linalg import null_space

def M(y0,y1,y2,y3,y4,y5):
    return np.array([[0, -y4, y3, y2, -y1, 0],
                     [42*y3, -13*y3 - y5, 0, 42*y0 - 13*y1 + 2*y3, 0, -y1],
                     [-y5, 0, 2*y2, 0, 0, -y0],
                     [0, 42*y2 - 13*y4, 42*y1-y5, y4, -13*y1+y3, -y2]])

def F(x1,y1,z1,x2,y2,z2):
    M1 = M(z1**2,x1*z1,y1*z1,x1**2,x1*y1,y1**2)
    M2 = M(z2**2,x2*z2,y2*z2,x2**2,x2*y2,y2**2)
    m = (y2*z1-y1*z2)/(x2*z1-x1*z2)
    x3 = (13*z1*z2 - x2*z1 - x1*z2) + z1*z2*m**2
    y3 = y1*z2 + m*(x3-x1*z2)
    z3 = z1*z2
    M3 = M(z3**2,x3*z3,y3*z3,x3**2,x3*y3,y3**2)
    A = np.concatenate((M1,M2))
    B = A[:-1]
    D1 = np.linalg.det(np.concatenate((B[:0],B[1:])))
    D2 = np.linalg.det(np.concatenate((B[:1],B[2:])))
    D3 = np.linalg.det(np.concatenate((B[:2],B[3:])))
    D4 = np.linalg.det(np.concatenate((B[:3],B[4:])))
    D5 = np.linalg.det(np.concatenate((B[:4],B[5:])))
    D6 = np.linalg.det(np.concatenate((B[:5],B[6:])))
    D7 = np.linalg.det(np.concatenate((B[:6],B[7:])))
    C = A[1:]
    E1 = np.linalg.det(np.concatenate((C[:0],C[1:])))
    E2 = np.linalg.det(np.concatenate((C[:1],C[2:])))
    E3 = np.linalg.det(np.concatenate((C[:2],C[3:])))
    E4 = np.linalg.det(np.concatenate((C[:3],C[4:])))
    E5 = np.linalg.det(np.concatenate((C[:4],C[5:])))
    E6 = np.linalg.det(np.concatenate((C[:5],C[6:])))
    E7 = np.linalg.det(np.concatenate((C[:6],C[7:])))
    P = np.concatenate((np.array([D1*A[0]-D2*A[1]+D3*A[2]-D4*A[3],
                                  E1*A[1]-E2*A[2]+E3*A[3]]), M3))
    # remove the 1st column to get 6x5; should hopefully be rank 5 still
    P2 = P[:,1:]
    F1 = np.linalg.det(np.concatenate((P2[:0],P2[1:])))
    F2 = np.linalg.det(np.concatenate((P2[:1],P2[2:])))
    F3 = np.linalg.det(np.concatenate((P2[:2],P2[3:])))
    F4 = np.linalg.det(np.concatenate((P2[:3],P2[4:])))
    F5 = np.linalg.det(np.concatenate((P2[:4],P2[5:])))
    F6 = np.linalg.det(np.concatenate((P2[:5],P2[6:])))
    # F1M1 - F2M2 + ... - F6M6 = 0
    H = F1*P[0] - F2*P[1]
    return list([a/H[0] for a in H])

def rand_elliptic_pt():
    # return a random point on E \subset P^5
    # choose random x, then y determined from that
    x = random()*20
    z = 1.
    y = sqrt(x*(x-6.)*(x-7.))
    return (x,y,z)

def to_p5(x,y,z):
    return (z**2,x*z,y*z,x**2,x*y,y**2)

out = []
for i in range(20):
    x1,y1,z1 = rand_elliptic_pt()
    x2,y2,z2 = rand_elliptic_pt()
    w0,w1,w2,w3,w4,w5=tuple(F(x1,y1,z1,x2,y2,z2))
    '''out.append([w0**2,w0*w1,w0*w2,w0*w3,w0*w4,w0*w5,
                w1**2,w1*w2,w1*w3,w1*w4,w1*w5,
                w2**2,w2*w3,w2*w4,w2*w5,
                w3**2,w3*w4,w3*w5,
                w4**2,w4*w5,
                w5**2])'''
    out.append([w0,w1,w2,w3,w4,w5])
out = np.array(out)
for sol in out:
    print('{' + ','.join([str(x).replace('j','*ii') for x in sol]) + '}')
#A = np.array([[symbols('a{}{}'.format(i,j)) for i in range(6) for j in range(6)]])



'''
l = [0,
     0,
     0,
 (-1180.29383+14525.029j),
 (78299.5993+23275.0463j),
 (-13010.5105+27216.6508j),
 (295.073477-3631.25724j),
 (-39149.7996-11637.5231j),
     0,
 (-2298.20002-2219.00616j),
 (-11622.3036-14006.2762j),
 (3252.6277-6804.1627j),
 (4596.40013+4438.01233j),
 (5811.15185+7003.13811j),
     0,
     0,
     0,
 (-762.359324-996.302424j),
 (190.589828+249.075606j),
     0,
     0]
E = [sum([x*y for x,y in zip(l, A)]) for A in out]
#import pdb;pdb.set_trace()
print(max([abs(y) for y in E]))
#print(np.linalg.matrix_rank(out, tol=.000001))
#print(100000*np.transpose(null_space(out, rcond=.00000000000000001)))
'''
