# Sab(mu) = N_{mu, (x,1,...) for x in a, (x,1,...) for x in b} with b's fixed

from hurwitz import S, f, auts
import functools

def unique_enum(l):
    return {y:x for x,y in enumerate(l)}

def special(n,m,z1,z2,a,mu):
    if z1+z2-n-m-3<=0:
        return 0
    X = f(z1,n,m,z1+z2-3) * Sab(a,[z1+z2-n-m-3],mu)
    X *= (auts([n,m,z1+z2-n-m-3])/auts([z1,z2]))
    return X

def threes(mu):
    return ns(mu, 3)

def fours(mu):
    return ns(mu, 4)

def ns(mu, n):
    # family where ramification is mostly (n,1,...,1)
    d = sum(mu)
    assert(all([x>0 for x in mu]))
    assert((d-2)%(n-2)==0)
    k = int((d-2)/(n-2)) # number of appearances of (n,1,...,1)
    twos = int(d+len(mu)-k*(n-1))
    if twos<0:
        return 0
    return Sab([n]*k + [2]*twos, [], mu)

def Sab(a, b, mu, log=False):
    return Sab_helper(tuple(a), tuple(b), tuple(mu), log)

@functools.cache
def Sab_helper(a, b, mu, log=False):
    output = (lambda x: print(x)) if log else (lambda x: None)
    d = sum(mu)
    output((a,b,mu))
    if any([x<1 for x in b]):
        return 0
    b = [x for x in b if x>1]
    a = list(a)
    mu = list(mu)
    if any([x<=0 for x in mu]):
        return 0
    elif a==[2,2] and (not b) and (not mu):
        return 1
    assert len(a) == len(mu)+2, '|a|=|mu|+2 is necessary for dimension reasons'
    assert sum(a)-len(a)+sum(b)-len(b)+sum(mu)-len(mu)==2*d, 'Riemann-Hurwitz failed'

    if b:
        if not mu:
            return 0
        # apply Claim 1
        n = mu[-1]
        mu2 = mu[:-1]
        y = b[0]
        b2 = b[1:]
        output('applying claim 1')
        X = Sab(a, b2, mu2+[n-y+1])

        for z, i in unique_enum(a).items():
            a2 = a[:i]+a[i+1:]
            output('applying claim 1')
            Y = min([n,y+z-n-2,z-1,y-1]) * Sab(a2,b2+[y+z-n-2],mu2)
            X += Y
        return X
    else: # b is empty
        if len(mu)==1: # base case
            output('base case: {}'.format(a))
            assert(len(a)==3)
            return S(*a)
        # apply Claim 2
        n = mu[-2]
        m = mu[-1]
        mu2 = mu[:-2]
        X = 0

        for z, i in unique_enum(a).items():
            a2 = a[:i]+a[i+1:]
            output('adding term for z={}: Sab({},[],{})'.format(z,a2,mu2+[n+m-z+2]))
            Y = min([n,m,z-1,n+m-z+1])*Sab(a2,[],mu2+[n+m-z+2])
            X += Y

        for z1, i in unique_enum(a).items():
            if z1<3:
                continue
            # could we have z1=z2?
            j = a.index(z1)
            if j != i: # multiple z1s in a
                output('applying Claim 2 with z1=z2={} for {},{},{}'.format(z1,a,b,mu))
                a2 = a[:j]+a[j+1:i]+a[i+1:]
                X += special(n,m,z1,z1,a2,mu2)
                output(special(n,m,z1,z1,a2,mu2))
            for z2, j in unique_enum(a).items():
                if i>=j or z2<3:
                    continue
                output('applying Claim 2 with z1={},z2={} for {},{},{}'.format(z1,z2,a,b,mu))
                a2 = a[:i]+a[i+1:j]+a[j+1:]
                X += special(n,m,z1,z2,a2,mu2)
                output(special(n,m,z1,z2,a2,mu2))
        return X

if __name__ == '__main__':
    for k in range(3,12):
        t = 1

        d = k+t+1
        a = 5
        b = d-a
        l = a+2*b+2-t-2*k
        print(Sab([3]*k+[2]*l, [t], [a]+[1]*b))
    '''
    k = 30
    for n in range(3,10):
        A = ns([1]*((n-2)*(k-1)+2),n)
        B = ns([1]*((n-2)*(k)+2),n)
        print((n,B/A))
    '''
