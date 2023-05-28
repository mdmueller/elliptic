f = open('elliptic_output2.txt')
s = f.read()
f.close()

lines = s.split('\n')
data2 = {}
for i, line in enumerate(lines[:-1]):
    x = line.replace('[', '').replace(']', '').split(',')
    data2[((i/(32*32*32))%32, (i/(32*32))%32, (i/32)%32, i%32)] = float(x[-1])

data = data2 #{x:data2[x] for x in data2 if x[0]!= 13}


def neighbors(I):
    a,b,c,d = I
    L = []
    for x1 in [-1,0,1]:
        A = a+x1
        for x2 in [-1,0,1]:
            B = b+x2
            for x3 in [-1,0,1]:
                C = c+x3
                for x4 in [-1,0,1]:
                    D = d+x4
                    if any([(y < 0 or y>=32) for y in (A,B,C,D)]) or (A,B,C,D)==(a,b,c,d) or (A,B,C,D) not in data:
                        continue
                    L.append((A,B,C,D))
    return L

T = [(a,b,c,d) for (a,b,c,d) in data if all([data[neighbor]<data[(a,b,c,d)] for neighbor in neighbors((a,b,c,d))])]
import pdb;pdb.set_trace()
