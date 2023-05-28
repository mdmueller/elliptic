-- like p5.m2, but for any embedding in P^n

loadPackage "Divisor"
loadPackage "MultiprojectiveVarieties"

d = 4 -- embed E in P^(2d-1) by O(2d*p)
S = QQ[X,y,z]
R = S/ideal(y^2*z-X*(X-6*z)*(X-7*z))
p = divisor(ideal(z,X)) -- use capital X to avoid clash with x in ring U later
F = mapToProjectiveSpace(d*p, Variable => w) -- R <- Q[w_1,...,w_d]

v2 = veronese(d-1, 2) -- P^d -> P^{(d+1 choose 2)-1}
T = target v2 -- P^d
w = gens(source F)
phi = map(source F, T, w)
G = F*phi*v2 -- E -> P^d -> P^{(d+1 choose 2)-1} = P^N
I = trim(kernel G) -- represents image of E in P^N, contained in some linear space
D = degrees I
assert(D == sort(D))
linforms = {}	 
for i from 1 to #D do (
    if D_(i-1)=={1} then linforms = append(linforms, I_(i-1)))
H = ideal(linforms) -- the largest linear space in P^N containing E; should be of dimension 2d-1
assert(dim(H) == 2*d)

-- S_1 = {H \in (P^(2d-1))*: H\cap E has even multiplicities}
-- we can get S_1 as the image of T_1 under projection P^N -> P^(2d-1)
-- T_1 is the image of Veronese w2: (P^d)* -> (P^N)*
U = source v2 -- P^N
x = gens U
forms = {}
for i from 1 to #x do (
    forms = append(forms, 2*x_(i-1))
)
j = 0
for i from 1 to d do (
    forms = replace(j, 1/2*forms_j, forms);
    j = j+(d+1-i)
    )
w2 = v2*map(U,U,forms)
T1 = kernel w2
