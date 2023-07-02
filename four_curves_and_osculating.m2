-- look at intersection of images of (2,2) and (3,1) loci under projection away from H_0

load "osculating.m2"
Dv2 = D
load "four_curves_linear.m2"
D = D1+D2+D3+D4

x=symbol x
y=symbol y
z=symbol z
R = QQ[x,y,z]

Dosc = (map(R,ring(ideal(Dv2)),{x,y,z})) (ideal(Dv2))
Dfour = (map(R,ring(ideal(D)),{x,y,z})) (ideal(D))
Z = projectiveVariety(radical(Dosc+Dfour))
