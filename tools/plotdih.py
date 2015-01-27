#!/usr/bin/env python

import sys
import numpy as np
import matplotlib.pyplot as plt

if len(sys.argv) != 5:
    print 'Plot OPLS cosine dihedral function'
    print 'usage: plotdih.py V1 V2 V3 V4'
    sys.exit(0)

v1, v2, v3, v4 = [ float(v) for v in sys.argv[1:5] ]

phi = np.linspace(0, 360)

u = v1/2.*(1.0+np.cos(phi*np.pi/180.0)) + \
    v2/2.*(1.0-np.cos(2*phi*np.pi/180.0)) + \
    v3/2.*(1.0+np.cos(3*phi*np.pi/180.0)) + \
    v4/2.*(1.0-np.cos(4*phi*np.pi/180.0))

plt.plot(phi,u)

plt.xlim(0, 360)
plt.xticks(np.linspace(0, 360, 13))
plt.grid(True)

plt.show()
