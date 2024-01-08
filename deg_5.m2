loadPackage "MultiprojectiveVarieties"

X = Grassmannian(1, 4, CoefficientRing=>QQ) -- dim 6 in P^9
setRandomSeed 1234
R = QQ[t0, t1, t2, t3, t4] -- P^4
M = random(QQ^10,QQ^5)
T = transpose(matrix{{t0,t1,t2,t3,t4}})
F = rationalMap {transpose(M*T)} -- generic inclusion P^4 -> P^9
t = gens target F
phi = map(target F, ring X, t)
X = phi(X)
E = projectiveVariety(F ^* X) -- E = X \cap P^4, cut out by 5 equations
t = matrix {gens ring ideal E}
F = (gens(ideal E))_{0,1,2}  -- choose 3 of 5 equations

--M = diff(t, transpose F) || random(ZZ^1,ZZ^5)
M = diff(t, transpose F) || matrix{{1,0,0,0,0}}
D1 := det(submatrix'(M, {0}))
D2 := det(submatrix'(M, {1}))
D3 := det(submatrix'(M, {2}))
D4 := det(submatrix'(M, {3}))
D5 := det(submatrix'(M, {4}))
vprime = matrix{{D1, -D2, D3, -D4, D5}} -- v'(0)

hess := G -> vprime*diff(transpose(t)*t, G)*transpose(vprime)
--M = M || random(ZZ^1,ZZ^5)
M = M || matrix{{0,1,0,0,0}}
--N = -(transpose(matrix{{hess(F_(0,0)), hess(F_(0,1)), hess(F_(0,2))}})) || random(ZZ^2,ZZ^1) -- N is 5x1
N = -(transpose(matrix{{hess(F_(0,0)), hess(F_(0,1)), hess(F_(0,2))}})) || matrix{{0},{0}} -- N is 5x1
-- M*v''(0) = N
D = det(M)
L = {}
for i from 0 to 4 do (
    L = append(L, det(M_{0..(i-1)} | N | M_{(i+1)..4}))/D)
vprime2 = matrix{L} -- v''(0)
print("hi")

hess2 := G -> (
    A := 3*vprime2*diff(transpose(t)*t, G)*transpose(vprime);
    for j from 0 to 4 do (
	for k from 0 to 4 do (
	    print(k);
	    for l from 0 to 4 do (
		A = A + diff(t_(0,j)*t_(0,k)*t_(0,l), G)*vprime_(0,j)*vprime_(0,k)*vprime_(0,l))));
    A)
--N = -(transpose(matrix{{hess2(F_(0,0)), hess2(F_(0,1)), hess2(F_(0,2))}})) || random(ZZ^2,ZZ^1) -- N is 5x1
H1 = hess2(F_(0,0))
print("done 1")
H2 = hess2(F_(0,1))
print("done 2")
H3 = hess2(F_(0,2))
print("done 3")
N = -(transpose(matrix{{H1, H2, H3}})) || matrix{{0},{0}} -- N is 5x1
-- M*v'''(0) = N
D = det(M)
L = {}
for i from 0 to 4 do (
    L = append(L, det(M_{0..(i-1)} | N | M_{(i+1)..4}))/D)
vprime3 = matrix{L} -- v''(0)
