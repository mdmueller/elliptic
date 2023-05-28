-- find hyperplanes in P^5 tangent to E (see p5.m2) at a_4 = q1 with order 5
-- y^2 = x(x-6)(x-7)

R = QQ[h0,h1,h2,h3,h4,h5]
S = R[t]
phi = map(R,S)
x = t+3 -- t = x-3
y = 6 - 3/4*t - 73/192*t^2 + 55/1536*t^3 - 3349/442368*t^4 -- degree 4 taylor poly of sqrt(x(x-6)(x-7)) centered at x=3
f = h0 + h1*x + h2*y + h3*x^2 + h4*x*y + h5*y^2
M = entries(((coefficients(f))_1)_0) -- coeffs of f
Ord1 = ideal(phi(M_(-1))) -- hyperplanes with order >=1 contact with E at a (i.e., they contain a)
Ord2 = Ord1 + ideal(phi(M_(-2))) -- hyperplanes tangent to E at a
Ord3 = Ord2 + ideal(phi(M_(-3))) -- ord >= 3 tangency
Ord4 = Ord3 + ideal(phi(M_(-4))) -- ord >= 4

