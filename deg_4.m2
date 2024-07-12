load "four_curves_linear.m2"
load "osculating.m2"

-- get point from ideal in R
getPoint = I -> (
    use R;
    J := sub(I,{w=>1});
    {lift(x%J,k), lift(y%J,k), lift(z%J,k), 1_k})

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
T = k[x,y,z,w,X,Y,Z,W,Degrees=>{4:{1,0},4:{0,1}}] -- P^3 x (P^3)*
phiR = map(T,R)
phiS = map(T,S)
use S

-- import loci into (P^3)*
import = locus -> projectiveVariety((map(S, ring(ideal(locus)), {X,Y,Z,W})) ideal(locus))
X22 = import X22
X31 = import X31
-- using tangentLine(E,{x=>-2726, y=>-47, z=>453, w=>1})
-- which contains {x=>1162, y=>-1578, z=>2, w=>1}
X211 = projectiveVariety(ideal(-2726*X-47*Y+453*Z+W, 1162*X-1578*Y+2*Z+W))
-- using tangentLine(E,{x=>-347, y=>-1991, z=>-2918, w=>1})
-- which contains {x=>294, y=>1278, z=>5, w=>3}
X211p = projectiveVariety(ideal(-347*X-1991*Y-2918*Z+W, 294*X+1278*Y+5*Z+3*W))

labels = {"(4)", "(2,2)", "(3,1)", "(2,1,1)", "(1,1,1,1)"}
-- planes H_0
planes = {{Y,Z,W}, {X-2265*W, Y-2398*W, Z+2709*W}, {X+1134*W, Y+2478*W, Z+2455*W},
    {X+27/38*W, Y-6/38*W, Z+8*sqrt10/38*W}, {X+2*Y-3*Z+4*W, 7*X+5*Y-11*Z+9*W, 13*X+14*Y-15*Z+17*W}}
maps = {}

for i from 0 to 4 do (
    H0 = planes_i;
    use T;
    << "First partition: " << labels_i << endl;
    -- find pairs (H0, p) where p in E\cap H0
    A = projectiveVariety(phiR(ideal(F,G))) * projectiveVariety(phiS(ideal(H0))) * projectiveVariety(ideal(X*x+Y*y+Z*z+W*w));
    points = minimalPrimes(preimage(phiR, ideal(A))); -- points of E\cap H0
    piH0 = multirationalMap {rationalMap H0};
    maps = append(maps, piH0);
    bad = projectiveVariety(ideal(1_S));
    if i < 4 then (
    	for j from 0 to (#points - 1) do (
	    P = getPoint(points_j);
	    use S;
	    bad = bad + projectiveVariety(ideal(P_0*X+P_1*Y+P_2*Z+P_3*W));
	    );
	);
    bad = reduce(piH0(bad));
    Y22 = piH0(X22);
    Y31 = piH0(X31);
    Y211 = piH0(X211);
    Y211p = piH0(X211p);

    if i==0 then (
	fourtorsion = reduce((piH0 ^* (getcusps(Y31)))*X31); -- 4-torsion points
	Q31 = piH0(tCone(X31));
	use S;
	Q22 = piH0(tangentLine(X22, {X=>1, Y=>0, Z=>0, W=>0}));
	assert(Q31 == Q22);
	bad = bad+Q22; -- throw out this image of tangent lines at H_0
	);
    << "(3,1),(3,1): " << doublefibers(piH0|X31, reduce(singularLocus(Y31))-bad) << endl;
    << "(3,1),(2,2): " << degree(reduce(Y22*Y31)-bad-piH0(fourtorsion)) << endl;
    << "(2,2),(2,2): " << doublefibers(piH0|X22, reduce(singularLocus(Y22))-bad) << endl;
    << "(2,1,1)fixed,(3,1): " << degree(reduce(Y31*Y211)-bad-piH0(reduce(X211*X31)))<<endl;
    << "(2,1,1)fixed,(2,1,1)fixed: " << degree(reduce(Y211*Y211p)-bad-piH0(reduce(X211*X211p)))<<endl;	       
    )
