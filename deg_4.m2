load "four_curves_linear.m2"
load "osculating.m2"

--k = toField(QQ[sqrt10]/(sqrt10^2-10))
k = QQ; sqrt10= 10 --TODO: fix?

R = k[x,y,z,w] -- P^3
F = x^2+y^2+z^2-x*y-w^2
G = 3*x^2-10*x*y-5*x*w+4*y^2+z^2+6*y*w+2*w^2
E = projectiveVariety(ideal(F,G))

S = k[X,Y,Z,W] -- (P^3)*

-- import loci into (P^3)*
import = locus -> projectiveVariety((map(S, ring(ideal(locus)), {X,Y,Z,W})) ideal(locus))
X22 = import X22
X31 = import X31

-- take plane curve, return cusps
getcusps = C -> (
    F := (ideal(C))_0;
    t := gens(ring(F));
    f := sub(F, {t_2=>1}); -- dehomogenize
    lin1 := diff(t_0,f);
    lin2 := diff(t_1,f); -- these must vanish at a singular point
    a := 1/2*diff(t_0^2,f);
    b := diff(t_0*t_1,f);
    c := 1/2*diff(t_1^2,f);
    delta := b^2-4*a*c; -- vanishes if there's a cusp
    cusps := minimalPrimes(ideal(f,lin1,lin2,delta)))

project := (P,X) -> (multirationalMap {P})(X)

pi4 = rationalMap {Y,Z,W} -- project from H_0 = {x=0}
Y22 = project(pi4, X22)
Y31 = project(pi4, X31)

-- use the plane H_0 = {x+y+z+w=0} (meets E transversely) to define another map (P^3)*->P^2
pi1111 = rationalMap {X-W, Y-W, Z-W}

-- use the plane H_0 = {-2/3*x+5*y+w=0} (meets E at 2 distinct tangent points)
pi22 = rationalMap {X+2/3*W, Y-5*W, Z}

-- use the plane H_0 = {27x-6y+8sqrt(10)*z-38w=0} (meets E in (2,1,1))
pi211 = rationalMap {X+27/38*W, Y-6/38*W, Z+8*sqrt10/38*W}

-- use the plane H_0 = {-6451/7264*x+45/454*y+30*sqrt(10)/227*z+w=0} (meets E at a,b, one with order 3)
pi31 = rationalMap {X+6451/7264*W, Y-45/454*W, Z-30*sqrt10/(227)*W}
