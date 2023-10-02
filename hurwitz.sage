def hurwitz_helper(perms, d, G, result, perms_so_far):
    # count lists of permutations in G whose product is result
    if len(perms) == 0:
        if result != G.identity():
            return (0,0,0)
        auts = len([h for h in G if all([h*g*h.inverse()==g for g in perms_so_far])])
        return (1,1,auts) if all([any([g(n)!=n for g in perms_so_far]) for n in range(1,d+1)]) else (1,0,auts)
    sigma = perms[-1]
    count_total = 0
    count_connected = 0
    count_nonorbifold = 0
    for g in G:
        if g.cycle_type()==sigma:
            t, c, n = hurwitz_helper(perms[:-1], d, G, result*g.inverse(), perms_so_far+[g])
            count_total += t
            count_connected += c
            count_nonorbifold += n
    return (count_total, count_connected, count_nonorbifold)

def hurwitz_count(perms, d):
    # perms is a list of partitions sigma_1, ..., sigma_n
    # weighted count of x_1,...,x_n in S_d such that x_1*...*x_n=1 and x_i has cycle type sigma_i
    # return (total, connected)
    G = SymmetricGroup(d)
    total, connected, nonorbifold = hurwitz_helper(perms, d, G, G.identity(), [])
    return (total/factorial(d), connected/factorial(d), nonorbifold/factorial(d))

print(hurwitz_count([[2,2],[3,1],[3,1]],4))
