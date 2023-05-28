S = ['random({1,0,0},R)' for i in range(3)] + ['random({0,1,0},R)' for i in range(3)] + ['random({0,0,1},R)' for i in range(3)]
possibilities = list(set([', '.join(S[:i]+S[i+1:j]+S[j+1:]) for i in range(9) for j in range(i+1, 9)]))
for i, x in enumerate(possibilities):
    print('X{0} = projectiveVariety(I3 + ideal({1}))'.format(i, x))
