#!/usr/bin/env python

import sys

if len(sys.argv) != 5:
    print "Convert Gromacs RB dihedral coefficients to OPLS-AA"
    print "usage: gromdihe C1 C2 C3 C4"
    quit()

c1 = c2 = c3 = c4 = 0.0
c1, c2, c3, c4 = [ float(a) for a in sys.argv[1:] ]

v4 = -c4/4.0
v3 = -c3/2.0
v2 = -c2 - c4
v1 = -2.0*c1 - 1.5*c3

print v1, v2, v3, v4

