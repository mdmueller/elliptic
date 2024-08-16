from sage.all import SymmetricGroup
import functools
import pickle

def transitive(elts, d):
    # return True if elts generate a transitive subgroup of S(d) (acting on {1,...,d})
    orbit = set([1]) # orbit containing 1; build it iteratively
    while True:
        old_orbit = set(orbit)
        for g in elts:
            orbit |= {g(x) for x in orbit}
        if len(orbit) == len(old_orbit):
            break
    return True if len(orbit)==d else False

def hurwitz_helper(perms, d, G, result, perms_so_far, connected):
    # count lists of permutations in G whose product is result
    # perms should be a tuple of tuples
    if len(perms) == 0:
        if result != G.identity():
            return (0,0)
        auts = len([h for h in G if all([h*g*h.inverse()==g for g in perms_so_far])])
        if connected:
            if transitive(perms_so_far, d):
                print(perms_so_far)
            return (1,auts) if transitive(perms_so_far, d) else (0,0) # all([any([g(n)!=n for g in perms_so_far]) for n in range(1,d+1)]) else (0,0)
        return (1,auts)
    sigma = perms[-1]
    count_total = 0
    count_nonorbifold = 0
    for g in G:
        if tuple(g.cycle_type())==sigma:
            t, n = hurwitz_helper(perms[:-1], d, G, result*g.inverse(), perms_so_far+[g], connected)
            count_total += t
            count_nonorbifold += n
    return (count_total, count_nonorbifold)

@functools.cache
def hurwitz1(perms, d, connected):
    G = SymmetricGroup(d)
    with open('hurwitzknown.pickle', 'rb') as f:
        known = pickle.load(f)
    if (perms, d, connected) in known:
        return known[(perms, d, connected)]
    total, nonorbifold = hurwitz_helper(perms, d, G, G.identity(), [], connected)
    x = fact(d)
    total /= x
    nonorbifold /= x
    known[(perms, d, connected)] = (total, nonorbifold)
    with open('hurwitzknown.pickle', 'wb') as f:
        pickle.dump(known, f, pickle.HIGHEST_PROTOCOL)
    return (total, nonorbifold)

def hurwitz_count(perms, d, connected=True):
    # perms is a list of partitions sigma_1, ..., sigma_n
    # weighted count of x_1,...,x_n in S_d such that x_1*...*x_n=1 and x_i has cycle type sigma_i
    # if connected=True, only count connected covers
    # return (total, nonorbifold)
    if d==1:
        return (1,1)
    perms = tuple([tuple(sorted(x)[::-1]) for x in sorted(perms) if len(x)!=sum(x)]) # remove (1,...,1)
    return hurwitz1(perms, d, connected)

@functools.cache
def fact(n):
    result = 1
    for i in range(2,n+1):
        result *= i
    return result

def auts(perm):
    # return number of automorphisms of permutation
    factor = 1
    for x in set(perm):
        factor *= fact(perm.count(x))
    return factor

def marked_hurwitz(perms, d):
    # return marked orbifold count
    x = hurwitz_count(perms, d)[0]
    for perm in perms:
        x *= auts(perm)
    return x

def f(a,n,m,d):
    # return H((n,m,p),(a,1,...),(b,1,...)) where n+m+p=d and a+b-3=d, with a<=b and n<=m<=p
    assert(3<=a and a<=d)
    b = d-a+3
    a,b = min([a,b]),max([a,b])
    assert(1<=n and 1<=m)
    p = d-m-n
    n,m,p = min([n,m,p]),n+m+p-min([n,m,p])-max([n,m,p]),max([n,m,p])
    X = 2/auts([n,m,p])*len([(n1,m1) for n1 in range(1,n+1) for m1 in range(1,m+1) if a-p-n1<=m1 and m1<=a-1-n1])
    return X

def S(a1,a2,a3):
    # return S_{{a1,a2,a3},{}}(d) where d=a1+a2+a3-4
    A = (a1,a2,a3)
    d = a1+a2+a3-4
    val = 0
    for i,ai in enumerate(A):
        aj, ak = tuple(A[:i]+A[i+1:])
        for e in range(1, aj+ak-2):
            f = aj+ak-2-e # e+f = d-a1+2 = a2+a3-2
            if e>f:
                continue

            #Xp = hurwitz_count([[d],[ai]+[1]*(d-ai),[e,f]+[1]*(d-e-f)],d)[0]
            if [e,f]==[1,1]:
                X = 1/d
            elif 1 in [e,f]:
                X = 1
            elif e==f:
                X = (ai-1)/2
            else:
                X = ai-1
            #Y = hurwitz_count([[e,f],[aj]+[1]*(e+f-aj),[ak]+[1]*(e+f-ak)],e+f)[0]
            Y = min([e,f,aj-1,ak-1])/auts([e,f])
            '''
            X = marked_hurwitz([[d],[ai]+[1]*(d-ai),[e,f]+[1]*(d-e-f)],d)
            Y = marked_hurwitz([[e,f],[aj]+[1]*(e+f-aj),[ak]+[1]*(e+f-ak)],e+f)
            '''
            N = (aj+ak-2)*X*Y*fact(d-ai)*auts([e,f]+[1]*(d-e-f))*fact(d-aj)*fact(d-ak)/fact(ai-2)
            #print(i,e,f,X,Y,N)
            val += N
    return val*2/(auts(A)*fact(d-a1)*fact(d-a2)*fact(d-a3))

if __name__ == '__main__':
    first=[int(x) for x in input('mu1? ').split(',')]
    second=[int(x) for x in input('mu2? ').split(',')]
    third=[int(x) for x in input('mu3? ').split(',')]
    deg = sum(first)
    assert(deg==sum(second))
    assert(deg==sum(third))
    print(hurwitz_count([first,second,third],deg))
