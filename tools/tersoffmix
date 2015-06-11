#!/usr/bin/env python

import sys, math

class tersoff:
    '''Tersoff potential parameters'''

    def __init__(self, line):
        '''init from a string following LAMMPS convention'''

        tok = line.strip().split()
        if len(tok) < 17:
            print('Error: Tersoff parameter lines should have 17 fields.\n')
            print('See documentation of LAMMPS pair_style tersoff.')
            sys.exit(1)

        self.param = {
            'atom1': tok[0],
            'atom2': tok[1],
            'atom3': tok[2],
            'm': float(tok[3]),
            'gamma': float(tok[4]),
            'lambda3': float(tok[5]),
            'c': float(tok[6]),
            'd': float(tok[7]),
            'h': float(tok[8]),
            'n': float(tok[9]),
            'beta': float(tok[10]),
            'lambda2': float(tok[11]),
            'B': float(tok[12]),
            'R': float(tok[13]),
            'D': float(tok[14]),
            'lambda1': float(tok[15]),
            'A': float(tok[16])
            }

    def __setitem__(self, key, value):
        self.param[key] = value

    def __getitem__(self, key):
        return self.param[key]
    
    def __str__(self):
        return "{0:<2s} {1:<2s} {2:<2s} {3:3.1f} {4:3.1f} {5:3.1f} {6:8.1f} "\
               "{7:8.4f} {8:8.4f} {9:8.5f} {10:11.4e} {11:8.4f} {12:8.2f} "\
               "{13:6.3f} {14:6.3f} {15:8.4f} {16:8.2f}".format(\
                self['atom1'], self['atom2'], self['atom3'],
                self['m'], self['gamma'], self['lambda3'],
                self['c'], self['d'], self['h'],
                self['n'], self['beta'], self['lambda2'],
                self['B'], self['R'], self['D'], self['lambda1'], self['A'])
    
    def mix(self, other, chi = 1.0, atom3 = 'X'):
        new = tersoff(str(self))
        new['atom2'] = other['atom1']
        new['atom3'] = atom3
        new['lambda1'] = 0.5 * (self['lambda1'] + other['lambda1'])
        new['lambda2'] = 0.5 * (self['lambda2'] + other['lambda2']) # mu
        new['A'] = math.sqrt(self['A'] * other['A'])
        new['B'] = math.sqrt(self['B'] * other['B']) * chi
        Ri = self['R'] - self['D']
        Si = self['R'] + self['D']
        Rj = other['R'] - other['D']
        Sj = other['R'] + other['D']
        Rij = math.sqrt(Ri * Rj)
        Sij = math.sqrt(Si * Sj)
        new['R'] = 0.5 * (Sij + Rij)
        new['D'] = 0.5 * (Sij - Rij)
        return new

    def header(self):
        return "# i j k  m   gam lam3     c     d        h=cos0  n        beta"\
          "         lam2       B     R      D       lam1       A"

    
class chi_ij:
    '''binary interaction chi_ij'''

    def __init__(self, line):
        tok = line.strip().split()
        if len(tok) < 4:
            print('Error: chi_ij lines should have 4 fields.')
            sys.exit(1)

        self.data = {'atom1': tok[1], 'atom2': tok[2], 'chi': float(tok[3])}

    def __setitem__(self, key, value):
        self.data[key] = value

    def __getitem__(self, key):
        return self.data[key]

    def __str__(self):
        return("{0:<2s} {1:<-2s} {2:6.4f}".format(self['atom1'],
                                                  self['atom2'], self['chi']))


def main():

    if len(sys.argv) < 2:
        print("Mixing rules for Tersoff potential parameters.")
        print("J Tersoff, Phys Rev B 39 (1989) 5566.")
        print("usage: tersoffmix.py [-p] file")
        print("  The format for each line for a given atom type is described")
        print("  in the documentation of LAMMPS pair_style tersoff.")
        print("  Unlike interaction parameters are given in additional lines:")
        print("  chi atom_i atom_j value")
        sys.exit(1)

    if sys.argv[1] == '-p':
        justprint = True
        infile = sys.argv[2]
    else:
        justprint = False
        infile = sys.argv[1]

    atom = []
    atomlist = []
    chi = []
    with open(infile, 'r') as f:
        line = f.readline()
        while line:
            tok = line.strip().split()
            if len(tok) and tok[0] != '#':
                if tok[0].lower() == 'chi':
                    chi.append(chi_ij(line))
                else:
                    if justprint:
                        atom.append(tersoff(line))
                    else:
                        if tok[0] not in atomlist:
                            atomlist.append(tok[0])
                            atom.append(tersoff(line))
            line = f.readline()

    if justprint:
        print(atom[0].header())
        for at in atom:
            print(str(at))
        sys.exit(0)

    for x in chi:
        print('# chi ' + str(x))

    print atom[0].header()
    for ati in atom:
        for atj in atom:
            xij = 1.0
            for x in chi:
                if x['atom1'] == ati['atom1'] and \
                    x['atom2'] == atj['atom1'] or \
                    x['atom2'] == ati['atom1'] and \
                    x['atom1'] == atj['atom1']:
                    xij = x['chi']
            for atk in atom:
                at = ati.mix(atj, chi = xij, atom3 = atk['atom3'])
                print(str(at))

if __name__ == "__main__":
    main()
