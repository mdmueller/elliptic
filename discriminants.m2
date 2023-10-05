-- count certain maps P^1->P^1

loadPackage "Resultants"

K = QQ--ZZ/6151
r = 21

discriminants = (g, type) -> (
    coeffs := ((coefficients(g))_1)_0;
    D1 := affineDiscriminant(g);
    a := coeffs_0;
    b := coeffs_1;
    c := coeffs_2;
    d := coeffs_3;
    e := coeffs_4;
    D2 := c^2 - 3*b*d + 12*a*e;
    D3 := 64*a^3*e-16*a^2*c^2+16*a*b^2*c-16*a^2*b*d-3*b^4;
    -- see https://en.wikipedia.org/wiki/Quartic_function, quantities Delta_0 and D
    if type=="(3,1)" then (D1,D2) else if type=="(2,2)" then (D1,D3) else "??")

-- count maps P^1->P^1 with ramification (3,1) over infinity, (2,1,1) over 0, (3,1) over 1
-- where fiber of infinity is 3*[infinity]+[w], fiber of 0 is 2[0]+[1]+[r], r is fixed and w is variable
R = K[k,w][t]
-- f(t) = t^2*(t-1)*(t-r)/(k*(t-w))
g = t^2*(t-1)*(t-r)-k*(t-w)
I = ideal(discriminants(g, "(3,1)"))
<< "(3,1),(2,1,1),(3,1): " << degree(radical I) << endl

-- count maps P^1->P^1 with ramification (3,1) over infinity, (2,1,1) over 0, (2,2) over 1
-- where fiber of infinity is 3*[infinity]+[w], fiber of 0 is 2[0]+[1]+[r], r is fixed and w is variable
R = K[k,w][t]
-- f(t) = t^2*(t-1)*(t-r)/(k*(t-w))
g = t^2*(t-1)*(t-r)-k*(t-w)
I = ideal(discriminants(g, "(2,2)"))
I0=I
<< "(3,1),(2,1,1),(2,2): " << degree(radical I) << endl

-- count maps P^1->P^1 with ramification (2,2) over infinity, (2,1,1) over 0, (3,1) over 1
-- where fiber of infinity is 2[infinity]+2[w], fiber of 0 is 2[0]+[1]+[r], r is fixed and w is variable
R = K[k,w][t]
-- f(t) = t^2*(t-1)*(t-r)/(k*(t-w)^2)
g = t^2*(t-1)*(t-r)-k*(t-w)^2
I = ideal(discriminants(g, "(3,1)"))
<< "(2,2),(2,1,1),(3,1): " << degree(radical I) << endl

-- count maps P^1->P^1 with ramification (2,1,1) over infinity, (2,1,1) over 0, (3,1) over 1 and (3,1) over another unspecified point
-- where fiber of infinity is 2[infinity]+[a]+[b], fiber of 0 is 2[0]+[1]+[r], r is fixed and a,b are variable
S = K[k,l,a,b][t]
-- f(t) = t^2*(t-1)*(t-r)/(k(t-a)(t-b))
g = t^2*(t-1)*(t-r)-k*(t-a)*(t-b)
h = t^2*(t-1)*(t-r)-l*(t-a)*(t-b) -- want g,h to both have (3,1) roots
I = ideal(discriminants(g,"(3,1)"))+ideal(discriminants(h,"(3,1)"))
P = minimalPrimes(I)
-- then look at P: P_0 includes l=k (bad), P_1 and P_2 are points with a,b swapped
-- so 1 solution


-- count maps P^1->P^1 with ramification (2,1,1) over infinity, (2,1,1) over 0, (3,1) over 1 and (2,2) over another unspecified point
-- where fiber of infinity is 2[infinity]+[a]+[b], fiber of 0 is 2[0]+[1]+[r], r is fixed and a,b are variable
S = K[k,l,a,b][t]
-- f(t) = t^2*(t-1)*(t-r)/(k(t-a)(t-b))
g = t^2*(t-1)*(t-r)-k*(t-a)*(t-b)
h = t^2*(t-1)*(t-r)-l*(t-a)*(t-b) -- want g,h to both have (3,1) roots
I = ideal(discriminants(g,"(3,1)"))+ideal(discriminants(h,"(2,2)"))
--Q = minimalPrimes(I)

