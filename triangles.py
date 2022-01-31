from matplotlib import pyplot as plt
from matplotlib import interactive
from math import cos, sin, pi
from random import random

N = 2000 # number of tries
EPS = .000

def make_polygon(k):
    # k = number of sides
    return [(cos(2*pi*i/k), sin(2*pi*i/k)) for i in range(k)]

def make_point(poly):
    k = len(poly)
    # we want k numbers in [0,1] which add to 1
    deltas = [random()+.1 for i in range(k)]
    # now normalize
    size = sum(deltas)
    deltas = [x/size for x in deltas]
    x = sum(deltas[i]*poly[i][0] for i in range(k))
    y = sum(deltas[i]*poly[i][1] for i in range(k))
    return (x, y)

def make_points(k, n):
    # choose n random points in a regular k-gon
    poly = make_polygon(k)
    points = [make_point(poly) for i in range(n)]
    return (poly, points)

def between(a, b, c):
    # see if a is between b and c (inclusive)
    if b<=c:
        return b-EPS<=a and a<=c+EPS
    else:
        return c-EPS<=a and a<=b+EPS

def segments_intersect(seg1, seg2):
    # segments are tuples ((x1, y1), (x2, y2))
    m1 = (seg1[1][1]-seg1[0][1])/(seg1[1][0]-seg1[0][0])
    m2 = (seg2[1][1]-seg2[0][1])/(seg2[1][0]-seg2[0][0])
    if m1 == m2 and seg1 != seg2:
        return False
    # x corresponds to intersection of lines
    x = (seg1[0][1]-seg2[0][1]+m2*seg2[0][0]-m1*seg1[0][0])/(m2-m1)
    y = m1*(x-seg1[0][0])+seg1[0][1]
    return (between(x,seg1[0][0],seg1[1][0]) and between(y,seg1[0][1],seg1[1][1])
            and between(x,seg2[0][0],seg2[1][0]) and between(y,seg2[0][1],seg2[1][1]))

def intersection_pt(seg1, seg2):
    # just find the intersection pt
    m1 = (seg1[1][1]-seg1[0][1])/(seg1[1][0]-seg1[0][0])
    m2 = (seg2[1][1]-seg2[0][1])/(seg2[1][0]-seg2[0][0])
    if m1 == m2 and seg1 != seg2:
        return False
    # x corresponds to intersection of lines
    x = (seg1[0][1]-seg2[0][1]+m2*seg2[0][0]-m1*seg1[0][0])/(m2-m1)
    y = m1*(x-seg1[0][0])+seg1[0][1]
    return (x,y)


def num_regions(k, n):
    poly, ps = make_points(k, n)
    segments = [(p, q) for p in poly for q in ps]
    points = []
    for i in range(N):
        p = make_point(poly)
        for q in points:
            if not any([segments_intersect((p, q), seg) for seg in segments]):
                break
        else:
            points.append(p)
    #print(len(points))
    '''
    plt.plot([q[0] for q in ps], [q[1] for q in ps], 'bo')
    plt.plot([q[0] for q in poly]+[poly[0][0]], [q[1] for q in poly]+[poly[0][1]], 'go-')
    for seg in segments:
        plt.plot([seg[0][0],seg[1][0]],[seg[0][1],seg[1][1]], 'go-')
    for q in points:
        plt.plot([q[0]],[q[1]],'ro')
    plt.show()
    '''
    return len(points)

if __name__=="__main__":
    k = 5 # triangle
    n=4
    m=0
    for i in range(100):
        m=max(m,num_regions(k,n))
        print(m)
    #print(3, num_regions(k, 3))
    #for n in range(1, 10):
    #    print(n, num_regions(k, n))


