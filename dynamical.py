import matplotlib.pyplot as plt
from math import cos, sin, pi

VERTICAL = None
eps = .0000000000001

def perp_bisector(p, q):
    v = ((p[0]+q[0])/2., (p[1]+q[1])/2.)
    slope = (q[1]-p[1])/(q[0]-p[0]) if abs(q[0]-p[0])>eps else VERTICAL
    a = -1./slope if abs(slope)>eps else VERTICAL
    if slope == VERTICAL:
        a = 0.
    elif a == VERTICAL:
        return (VERTICAL, v[0]) # v[0] is x-value
    # y = y_0 + a(x-x_0) = ax + b
    b = v[1]-a*v[0]
    return (a,b)

def intersection(line1, line2):
    if line2[0] == VERTICAL:
        line1, line2 = line2, line1
    a1, b1 = line1
    a2, b2 = line2
    if a1 == VERTICAL:
        x = b1
        y = a2*x+b2
    else:
        x = (b2-b1)/(a1-a2)
        y = a1*x+b1
    return (x,y)
    
def circumcenter(p1, p2, p3):
    line1 = perp_bisector(p1, p2)
    line2 = perp_bisector(p1, p3)
    return intersection(line1, line2)

def f(p1, p2, p3, N):
    points = [p1, p2, p3]
    for i in range(N):
        points.append(circumcenter(points[-3], points[-2], points[-1]))
    return points

p3 = (1.,0.)
p1 = (0.,0.)
p2 = (cos(15*pi/180), sin(15*pi/180))
points = f(p1, p2, p3, 20)

plt.plot([p[0] for p in points], [p[1] for p in points], 'ro')
plt.show()
