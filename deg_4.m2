load "four_curves_linear.m2"
load "osculating.m2"

-- takes projective 0-dim ideal I, returns lengths of points
getLengths = I -> (
    K := primaryDecomposition(I);
    L := {};
    for j from 1 to #K do (
	J := K_(j-1);
	D1 := degree(projectiveVariety(J));
	D2 := degree(projectiveVariety(radical(J)));
	L := append(L, {D1,D2});
	);
    L)

-- take plane curve, return cusps --TODO: make homogeneous so we don't miss points...
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
    cusps := projectiveVariety(homogenize(ideal(f,lin1,lin2,delta), t_2)))

-- take reduction of variety
reduce = variety -> projectiveVariety(radical(ideal(variety)))

-- given a map P and a dim 0 subvariety V, find the # of points in V with fiber having at least 2 points
doublefibers := (P,V) -> (
    U := reduce V;
    L := minimalPrimes ideal U;
    result := 0;
    for i from 1 to #L do (
	p := L_(i-1); -- point(s) of V
	len := degree(projectiveVariety(p));
	pre := reduce(P ^* (projectiveVariety(p)));
	--print(len); print(degree(pre));
	if degree(pre) > len then result = result+len;
	);
    result)

k = ZZ/107
sqrt10 = 44 -- 44^2 = 10 mod 107

R = k[x,y,z,w] -- P^3
F = x^2+y^2+z^2-x*y-w^2
G = 3*x^2-10*x*y-5*x*w+4*y^2+z^2+6*y*w+2*w^2
E = projectiveVariety(ideal(F,G))

S = k[X,Y,Z,W] -- (P^3)*

-- import loci into (P^3)*
import = locus -> projectiveVariety((map(S, ring(ideal(locus)), {X,Y,Z,W})) ideal(locus))
X22 = import X22
X31 = import X31

pi4 = multirationalMap {rationalMap {Y,Z,W}} -- project from H_0 = {x=0}
Y22 = pi4(X22) -- three conics and a line
Y31 = pi4(X31) -- deg 10 curve with 15 cusps, 20 nodes
fourtorsion = reduce((pi4 ^* (getcusps(Y31)))*X31) -- 4-torsion points
-- n + c = length(reduced singular locus)
-- n + 2c = length(singular locus)
-- c = length(singular locus) - length(reduced singular locus)
-- n = length(reduced singular locus) - c = 2*length(reduced singular locus) - length(singular locus)
sing = singularLocus(Y31)
--print(2*degree(reduce(sing)) - degree(sing))
--print doublefibers((multirationalMap{pi4})|X31, sing)
print(degree(reduce(Y22*Y31))-16) -- subtract 16 for the images of 4-torsion points

-- use the plane H_0 = {13x+9y+7z+11w=0} (meets E transversely) to define another map (P^3)*->P^2
use S
pi1111 = multirationalMap {rationalMap {X+Y-2*W, Y+5*Z-4*W, Z+2*X-3*W}}
Z22 = pi1111(X22) -- deg 8, 24 nodes
Z31 = pi1111(X31) -- deg 12, 16 cusps (corresponding to 4-torsion), 38 nodes
sing = singularLocus(Z31)
--print doublefibers((multirationalMap{pi1111})|X31, sing)
print(degree(reduce(Z22*Z31))-16)

-- use the plane H_0 = {-2/3*x+5*y+w=0} (meets E at 2 distinct tangent points)
use S
pi22 = multirationalMap {rationalMap {X+2/3*W, Y-5*W, Z}}
W22 = pi22(X22)
W31 = pi22(X31) -- cusps are: TODO
sing = singularLocus(W31)
--print doublefibers((multirationalMap{pi22})|X31, sing)
print(degree(reduce(W22*W31))-16)

-- use the plane H_0 = {27x-6y+8sqrt(10)*z-38w=0} (meets E in (2,1,1))
use S
pi211 = multirationalMap {rationalMap {X+27/38*W, Y-6/38*W, Z+8*sqrt10/38*W}}
J22 = pi211(X22)
J31 = pi211(X31)
sing = singularLocus(J31)
--print doublefibers((multirationalMap{pi211})|X31, sing)
print(degree(reduce(J22*J31))-16)

-- use the plane H_0 = {-6451/7264*x+45/454*y+30*sqrt(10)/227*z+w=0} (meets E at a,b, one with order 3)
use S
pi31 = multirationalMap {rationalMap {X+6451/7264*W, Y-45/454*W, Z-30*sqrt10/(227)*W}}
K22 = pi31(X22)
K31 = pi31(X31)
sing = singularLocus(K31)
--print doublefibers((multirationalMap{pi31})|X31, sing)
print(degree(reduce(K22*K31))-16)
