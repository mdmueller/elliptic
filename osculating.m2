loadPackage "MultiprojectiveVarieties"
loadPackage "RationalPoints2"

R = QQ[x,y,z,w]
F = x^2+y^2+z^2-x*y-w^2
G = 3*x^2-10*x*y-5*x*w+4*y^2+z^2+6*y*w+2*w^2
E = projectiveVariety(ideal(F,G))

fx = 2*x-y
fy = 2*y-x
fz = 2*z
fxx = 2
fxy = -1
fxz = 0
fyy = 2
fyz = 0
fzz = 2
gx = 6*x-10*y-5
gy = 8*y-10*x+6
gz = 2*z
gxx = 6
gxy = -10
gxz = 0
gyy = 8
gyz = 0
gzz = 2

yx = (gz*fx-gx*fz)/(gy*fz-gz*fy)
zx = (gy*fx-gx*fy)/(gz*fy-gy*fz)

u = fxx + 2*fxy*yx + 2*fxz*zx + fyy*(yx)^2+2*fyz*yx*zx+fzz*(zx)^2
v = gxx + 2*gxy*yx + 2*gxz*zx + gyy*(yx)^2+2*gyz*yx*zx+gzz*(zx)^2

fw = -2
gw = -5*x+6*y+4

h = method()
h(frac R) := x -> homogenize(numerator(x),w)/homogenize(denominator(x),w)

H = multirationalMap{ rationalMap { h(u*fx - v*gx), h(u*fy - v*gy), h(u*fz - v*gz), h(u*fw - v*gw), 0 }} -- why is the 0 necessary?
C = H(E)
y = gens(ring(C))
F = multirationalMap {rationalMap {y_0,y_1,y_2,y_3}}
C = image(F)
