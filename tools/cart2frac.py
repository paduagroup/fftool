#!/usr/bin/env python

import sys

if len(sys.argv) < 4:
    print('Convert cartesian to factional coordinates (orthogonal lattice vectors)')
    print('usage: cart2frac a b c file.xyz')
    sys.exit(1)

a, b, c = [ float(d) for d in sys.argv[1:4] ]

with open(sys.argv[4], 'r') as f:
    print(f.readline().strip())
    print(f.readline().strip())
    line = f.readline()
    while line:
        tok = line.split()
        x, y, z = [ float(q) for q in tok[1:] ]
        print(" %-4s %12.7f %12.7f %12.7f" % (tok[0], x / a, y / b, z / c))
        line = f.readline()
