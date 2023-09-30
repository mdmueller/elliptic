load "four_curves_linear.m2"
load "osculating.m2"

-- get tangent line to C at p
tangentLine = (C, p) -> (
    R := ring(ideal(C));
    M := sub(jacobian(ideal(C)), p);
    projectiveVariety(ideal(matrix{gens(R)} * M)))

-- get tangent cone to C at the point [1:0:0:0]
tCone = C -> (
    R := ring(ideal(C));
    X := (gens(R))_0;
    Q := k[Y,Z,W];
    phi := map(Q, R, {1,Y,Z,W});
    phi2 := map(R, Q);
    I := phi ideal C;
    reduce projectiveVariety homogenize(phi2 tangentCone I, X))

-- takes projective 0-dim ideal I, returns lengths of points
getLengths = I -> (
    K := primaryDecomposition(I);
    L := {};
    D1 := 0;
    D2 := 0;
    J := 0;
    for j from 1 to #K do (
	J = K_(j-1);
	D1 = degree(projectiveVariety(J));
	D2 = degree(projectiveVariety(radical(J)));
	L = append(L, {D1,D2});
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
reduce = V -> projectiveVariety(radical(ideal(V)))

-- given a map P and a dim 0 subvariety V, find the # of points in V with fiber (excluding H0) having at least 2 points
doublefibers = (P,V) -> (
    U := reduce V;
    L := minimalPrimes ideal U;
    H0 := baseLocus(P);
    result := 0;
    for i from 1 to #L do (
	p := projectiveVariety(L_(i-1)); -- point(s) of V
	len := degree(p);
	pre := reduce(P ^* p); -- fiber of p
	if degree(pre) > len then result = result+len;
	);
    result)

k = ZZ/6151
sqrt10 = 2024 -- 2024^2 = 10 mod 6151

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
use S

print "(4)"
H0 = {X=>1, Y=>0, Z=>0, W=>0}
-- H0 \cap E = {[0:-1:0:1]}
bad = pi4(projectiveVariety(ideal(-Y+W))) -- undefined maps
Q31 = pi4(tCone(X31))
Q22 = pi4(tangentLine(X22, H0)) -- we see Q22==Q31
fourtorsion = reduce((pi4 ^* (getcusps(Y31)))*X31) -- 4-torsion points
sing = singularLocus(Y31)
print doublefibers(pi4|X31, sing-bad)
print(degree(reduce(Y22*Y31)-reduce(pi4(fourtorsion))-Q22-reduce(bad)))

-- use a plane H_0 that meets E transversely to define another map (P^3)*->P^2
print "(1,1,1,1)"
use S
pi1111 = multirationalMap {rationalMap {X+2*Y-3*Z+4*W, 7*X+5*Y-11*Z+9*W, 13*X+14*Y-15*Z+17*W}}
Z22 = pi1111(X22) -- deg 8, 24 nodes
Z31 = pi1111(X31) -- deg 12, 16 cusps (corresponding to 4-torsion), 38 nodes
sing = singularLocus(Z31)
print doublefibers(pi1111|X31, sing)
print(degree(reduce(Z22*Z31)-reduce(pi1111(fourtorsion))))

-- use the plane H_0 = {2265*x+2398*y-2709*z+w=0} (meets E at two distinct tangent points)
-- (random point in X22 obtained through four_curves_linear.m2)
-- H_0\cap E = {[1044:302:316:1], [-1609:2738:2829:1]}
print "(2,2)"
use S
H0 = {X=>2265, Y=>2398, Z=>-2709, W=>1}
pi22 = multirationalMap {rationalMap {X-2265*W, Y-2398*W, Z+2709*W}}
bad = pi22(projectiveVariety(ideal(1044*X+302*Y+316*Z+W)) + projectiveVariety(ideal(-1609*X+2738*Y+2829*Z+W))) -- undefined maps
W22 = pi22(X22)
Q22 = pi22(tangentLine(X22, H0))
W31 = pi22(X31) -- cusps are: TODO
sing = singularLocus(W31)
print doublefibers(pi22|X31, sing-bad)
print(degree(reduce(W22*W31)-reduce(pi22(fourtorsion))-reduce(bad)))

-- use the plane H_0 = {27x-6y+8sqrt(10)*z-38w=0} (meets E in (2,1,1))
-- H_0\cap E = {[1937:-2587:604:1], [2267:324:3016:1], [1904:-2229:-627:1]}
print "(2,1,1)"
use S
pi211 = multirationalMap {rationalMap {X+27/38*W, Y-6/38*W, Z+8*sqrt10/38*W}}
bad = pi211(projectiveVariety(ideal(1937*X-2587*Y+604*Z+W))+projectiveVariety(ideal(2267*X+324*Y+3016*Z+W))+projectiveVariety(ideal(1904*X-2229*Y-627*Z+W)))
J22 = pi211(X22)
J31 = pi211(X31)
sing = singularLocus(J31)
print doublefibers(pi211|X31, sing-bad)
print(degree(reduce(J22*J31)-reduce(pi211(fourtorsion))-reduce(bad)))

-- use the plane H_0 = {-1134*x-2478*y-2455*z+w=0} (meets E at a,b, one with order 3)
-- obtained from osculating.m2
-- H_0\cap E = {[-1737:327:2394:1], [2391:-2088:-1104:1]}
print "(3,1)"
use S
pi31 = multirationalMap {rationalMap {X+1134*W, Y+2478*W, Z+2455*W}}
bad = pi31(projectiveVariety(ideal(-1737*X+327*Y+2394*Z+W))+projectiveVariety(ideal(2391*X-2088*Y-1104*Z+W)))
K22 = pi31(X22)
K31 = pi31(X31)
sing = singularLocus(K31)
print doublefibers(pi31|X31, sing-bad)
print(degree(reduce(K22*K31)-reduce(pi31(fourtorsion))-reduce(bad)))
