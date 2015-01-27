#!/usr/bin/env python

import sys, pprint

c = [ [0.0 for i in range(6)] for j in range(6) ]

readeng = False
for line in sys.stdin:
    if readeng:
        tok = line.strip().split()
        eini = float(tok[0])
        efin = float(tok[2])
        readeng = False        
    if line.startswith('Elastic Constant'):
        tok = line.strip().split()
        i = int(tok[2][1]) - 1
        j = int(tok[2][2]) - 1
        c[i][j] = c[j][i] = float(tok[4])
        units = tok[5]
    elif '  Energy initial, next-to-last, final' in line:
        readeng = True

dashes = "".ljust(76, '-') 
print dashes
print "  Total lattice energy (initial) = %10.6f eV" % eini
print "  Total lattice energy (final)   = %10.6f eV" % efin
print dashes
print
print "  Elastic Constant Matrix (%s):" % units
print
print dashes
print "  ij",
for j in range(6):
    print "     %d     " % (j + 1),
print 
print dashes
for i in range(6):
    print "  %d " % (i + 1),
    for j in range(6):
        print " %10.4f" % c[i][j], 
    print
print dashes
