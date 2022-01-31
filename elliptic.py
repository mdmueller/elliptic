from sympy import *
x = symbols('x')
y = symbols('y')
z = symbols('z')
w = symbols('w')

f1 = 5*x**2+11*x*y+7*x*z+8*x*w+4*y**2+z**2+6*y*w+2*w**2
f2 = x**2+2*x*y+3*x*z+4*x*w+y**2+z**2-w**2

def find_planes(L):
    # given a line L, find all planes H containing L and tangent to E
    # return list of (x, H) where H is tangent to E at x

    # ok first, solve for x where T_xE intersects L
    a,b,c,d = L
    M = Matrix([[10*x+11*y+7*z+8*w, 11*x+8*y+6*w, 7*x+2*z, 8*x+6*y+4*w],
                [2*x+2*y+3*z+4*w,   2*x+2*y,      3*x+2*z, 4*x-2*w],
                [a,                 b,            c,       d],
                [1,                 0,            0,       0]])
    solns = nonlinsolve([M.det(), f1, f2, w-1], [x,y,z,w]) # TODO: work projectively or something...
    print(solns)
    import pdb;pdb.set_trace()
    #import pdb;pdb.set_trace()


# lines are in the form (a,b,c,d) corresponding to {ax+by+cz+dw=x=0}
find_planes((0.,1.,0.,0.))
