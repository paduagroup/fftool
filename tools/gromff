#!/usr/bin/env python

import sys

if len(sys.argv) < 4:
    print("Get OPLS force field coefficients from Gromacs database")
    print("usage: gromff ffbonded.itp i j [k [l]]")
    quit()

i = sys.argv[2]
j = sys.argv[3]
bond = True
angle = dihed = False
if len(sys.argv) >= 5:
    bond = False
    angle = True
    k = sys.argv[4]
if len(sys.argv) == 6:
    angle = False
    dihed = True
    l = sys.argv[5]

eps = 1.0e-9

f = open(sys.argv[1], 'r')
for line in f:
    tok = line.strip().split()

    if bond and len(tok) >= 5 and tok[2] == '1' and float(tok[3]):
        kr = r0 = 0.0
        if i == tok[0] and j == tok[1] or j == tok[0] and i == tok[1]:
            print(line)
            kr = float(tok[4]) / 100.0
            r0 = float(tok[3]) * 10.0
            print("{0:<3s} {1:<3s}  harm  {2:6.3f}  "\
                  "{3:7.1f}".format(i, j, r0, kr))
        
    elif angle and len(tok) >= 6 and tok[3] == '1' and float(tok[4]):
        kth = th0 = 0.0
        if i == tok[0] and j == tok[1] and k == tok[2] or \
                k == tok[0] and j == tok[1] and i == tok[2]:
            print(line)
            kth = float(tok[5])
            th0 = float(tok[4])
            print("{0:<3s} {1:<3s} {2:<3s}  harm  {3:6.1f} "\
                  "{4:7.1f}".format(i, j, k, th0, kth))

    elif dihed and len(tok) >= 11:
        v1 = v2 = v3 = v4 = 0.0
        if (i == tok[0] or tok[0] == 'X') and j == tok[1] and k == tok[2] and \
            (l == tok[3] or tok[3] == 'X') or \
            (l == tok[0] or tok[0] == 'X') and k == tok[1] and j == tok[2] and \
            (i == tok[3] or tok[3] == 'X'):
            if tok[4] == '3':
                print(line)
                c1, c2, c3, c4 = [ float(c) for c in tok[6:10] ]
                v4 = -c4/4.0 + eps
                v3 = -c3/2.0 + eps
                v2 = -c2 - c4 + eps
                v1 = -2.0*c1 - 1.5*c3 + eps
            elif tok[4] == '5':
                print(line)
                v1, v2, v3, v4 = [ float(v)+1.0e-6 for v in tok[5:9] ]
            print("{0:<3s} {1:<3s} {2:<3s} {3:<3s}  opls {4:9.4f} {5:9.4f} "\
                  "{6:9.4f} {7:9.4f}".format(i, j, k, l, v1, v2, v3, v4))
