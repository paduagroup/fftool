#!/usr/bin/env python

import sys, math

if len(sys.argv) < 2:
    print('Halgren estimation of angle force constants')
    print('usage: ktheta.py Z_i C_j Z_k r_ij/A r_jk/A theta_ijk/deg')
    print('result in kJ mol-1 rad-2')
    sys.exit(0)

zi, cj, zk, rij, rjk, th = [float(arg) for arg in sys.argv[1:]]

th *= math.pi / 180.0
d = ((rij - rjk)/(rij + rjk))**2
k = 1053.85 * zi * cj * zk / (rij + rjk) * th**(-2) * math.exp(-2.0*d)

print(k)

