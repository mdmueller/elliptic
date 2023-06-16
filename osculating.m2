loadPackage "MultiprojectiveVarieties"
loadPackage "RationalPoints2"

R = QQ[x,y,z,w]
F = x^2+y^2+z^2-x*y-w^2
G = 3*x^2-10*x*y-5*x*w+4*y^2+z^2+6*y*w+2*w^2
E = projectiveVariety(ideal(F,G))

fx = 2*x-y
fy = 2*y-x
fz = 2*z
fw = -2*w
fxx = 2
fxy = -1
fxz = 0
fxw = 0
fyy = 2
fyz = 0
fyw = 0
fzz = 2
fzw = 0
fww = -2
gx = 6*x-10*y-5*w
gy = 8*y-10*x+6*w
gz = 2*z
gw = -5*x+6*y+4*w
gxx = 6
gxy = -10
gxz = 0
gxw = -5
gyy = 8
gyz = 0
gyw = 6
gzz = 2
gzw = 0
gww = 4

M = matrix{{fx, fy, fz, fw}, {gx, gy, gz, gw}, {1, 0, 0, 0}}
D1 = det(submatrix'(M, {0}))
D2 = det(submatrix'(M, {1}))
D3 = det(submatrix'(M, {2}))
D4 = det(submatrix'(M, {3}))
vprime = matrix{{D1, -D2, D3, -D4}} -- v'(0)

-- Hessian
Hf = matrix{{ fxx, fxy, fxz, fxw }, { fxy, fyy, fyz, fyw }, { fxz, fyz, fzz, fzw }, { fxw, fyw, fzw, fww }}
Hg = matrix{{ gxx, gxy, gxz, gxw }, { gxy, gyy, gyz, gyw }, { gxz, gyz, gzz, gzw }, { gxw, gyw, gzw, gww }}
u = (vprime*Hg*transpose(vprime))_(0,0)
v = (vprime*Hf*transpose(vprime))_(0,0)

H = multirationalMap{ rationalMap { u*fx - v*gx, u*fy - v*gy, u*fz - v*gz, u*fw - v*gw, 0 }} -- why is the 0 necessary?
C = H(E)
y = gens(ring(C))
F = multirationalMap {rationalMap {y_0,y_1,y_2,y_3}}
C = image(F)
