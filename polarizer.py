#!/usr/bin/env python3
# polarizer.py - add Drude oscillators to LAMMPS data file.
# Agilio Padua <agilio.padua@univ-bpclermont.fr>, version 2015/04/04
# http://tim.univ-bpclermont.fr/apadua

"""
Add Drude oscillators to LAMMPS data file.

The format of file containing specification of Drude oscillators is:

# type  dm/u  dq/e  k/(kJ/molA2)  alpha/A3  thole
C3H     1.0   0.0   4184.0        2.051     2.6
...

where:
dm is the mass to put on the Drude particle (taken from its core),
dq is the charge to put on the Drude particle (taken from its core),
k is the garmonic force constant of the bond between Drude and core,
alpha is the polarizability,
thole is a Thole damping parameter.

A Drude particle is created for each atom in the LAMMPS data file
that corresponds to an atom type given in the Drude file.
Since LAMMPS uses numbers for atom types in the data file, a comment
after each line in the Masses section has to be introduced to allow
identification of the atom types within the force field database:

Masses

   1   12.011  # C3H
   2   12.011  # CTO
...


This script will add new atom types, new atoms, new bond types and
new bonds to the data file. It will also add a new seciton called
"Drudes" containing information about core-drude particle relations.

"""

import sys
import argparse
import math
from copy import deepcopy

# keywords of header and main sections (from data.py in Pizza.py)

hkeywords = ["atoms", "ellipsoids", "lines", "triangles", "bodies",
             "bonds", "angles", "dihedrals", "impropers",
             "atom types", "bond types", "angle types", "dihedral types",
             "improper types", "xlo xhi", "ylo yhi", "zlo zhi", "xy xz yz"]

skeywords = [["Masses", "atom types"],
             ["Pair Coeffs", "atom types"],
             ["Bond Coeffs", "bond types"], ["Angle Coeffs", "angle types"],
             ["Dihedral Coeffs", "dihedral types"],
             ["Improper Coeffs", "improper types"],
             ["BondBond Coeffs", "angle types"],
             ["BondAngle Coeffs", "angle types"],
             ["MiddleBondTorsion Coeffs", "dihedral types"],
             ["EndBondTorsion Coeffs", "dihedral types"],
             ["AngleTorsion Coeffs", "dihedral types"],
             ["AngleAngleTorsion Coeffs", "dihedral types"],
             ["BondBond13 Coeffs", "dihedral types"],
             ["AngleAngle Coeffs", "improper types"],
             ["Atoms", "atoms"], ["Ellipsoids", "ellipsoids"],
             ["Lines", "lines"], ["Triangles", "triangles"],
             ["Bodies", "bodies"],
             ["Bonds", "bonds"],
             ["Angles", "angles"], ["Dihedrals", "dihedrals"],
             ["Impropers", "impropers"], ["Velocities", "atoms"],
             ["Molecules", "atoms"], ["Drudes", "atoms"]]

class Data:
    
    def __init__(self, datafile):
        '''read LAMMPS data file (from data.py in Pizza.py)'''
        
        self.nselect = 1

        f = open(datafile, "r")

        self.title = f.readline()
        self.names = {}

        headers = {}
        while 1:
            line = f.readline().strip()
            if len(line) == 0:
                continue
            found = 0
            for keyword in hkeywords:
                if keyword in line:
                    found = 1
                    words = line.split()
                    if keyword == "xlo xhi" or keyword == "ylo yhi" or \
                      keyword == "zlo zhi":
                        headers[keyword] = (float(words[0]), float(words[1]))
                    elif keyword == "xy xz yz":
                        headers[keyword] = \
                          (float(words[0]), float(words[1]), float(words[2]))
                    else:
                        headers[keyword] = int(words[0])
            if not found:
                break
    
        sections = {}
        while 1:
            found = 0
            for pair in skeywords:
                keyword, length = pair[0], pair[1]
                if keyword == line:
                    found = 1
                    if length not in headers:
                        raise RuntimeError("data section {} has no matching"\
                                           " header value".format(line))
                    f.readline()
                    list = []
                    for i in range(headers[length]):
                        list.append(f.readline())
                    sections[keyword] = list
            if not found:
                raise RuntimeError("invalid section {} in data"\
                                   " file".format(line))
            f.readline()
            line = f.readline()
            if not line:
                break
            line = line.strip()
        
        f.close()
        self.headers = headers
        self.sections = sections
        
    def write(self, filename):
        '''write out a LAMMPS data file (from data.py in Pizza.py)'''

        with open(filename, "w") as f:
            f.write(self.title + '\n')
            for keyword in hkeywords:
                if keyword in self.headers:
                    if keyword == "xlo xhi" or keyword == "ylo yhi" or \
                       keyword == "zlo zhi":
                        pair = [ str(p) for p in self.headers[keyword] ]
                        f.write(pair[0] + ' ' + pair[1] + ' ' + keyword + '\n')
                    elif keyword == "xy xz yz":
                        triple = [ str(t) for t in self.headers[keyword] ]
                        f.write(triple[0] + ' ' + triple[1] + ' ' + triple[2] +
                                ' ' + keyword + '\n')
                    else:
                        f.write(str(self.headers[keyword]) + ' ' +
                                keyword + '\n')
            for pair in skeywords:
                keyword = pair[0]
                if keyword in self.sections:
                    f.write("\n{}\n".format(keyword))
                    for line in self.sections[keyword]:
                        f.write(line)

    def atomtypes(self):
        """extract atom IDs from data"""
        atomtypes = []
        for line in self.sections['Masses']:
            tok = line.split()
            if len(tok) < 4:
                print("warning: missing type for atom ID " + tok[0])
                continue
            atomtype = {}
            atomtype['id'] = int(tok[0])
            atomtype['m'] = float(tok[1])
            atomtype['type'] = tok[3]
            atomtypes.append(atomtype)
        return atomtypes

    def bondtypes(self):
        """extract bond types from data"""
        bondtypes = []
        for line in self.sections['Bond Coeffs']:
            tok = line.split()
            bondtype = {}
            bondtype['id'] = int(tok[0])
            bondtype['k'] = float(tok[1])
            bondtype['r0'] = float(tok[2])
            note = ''
            for i in range(3, len(tok)):
                note += tok[i] + ' '
            bondtype['note'] = note.strip()
            bondtypes.append(bondtype)
        return bondtypes

    def atoms(self):
        """extract atom registers and atom IDs from data"""
        atoms = []
        for line in self.sections['Atoms']:
            tok = line.split()
            atom = {}
            atom['n'] = int(tok[0])
            atom['mol'] = int(tok[1])
            atom['id'] = int(tok[2])
            atom['q'] = float(tok[3])
            atom['x'] = float(tok[4])
            atom['y'] = float(tok[5])
            atom['z'] = float(tok[6])
            note = ''
            for i in range(7, len(tok)):
                note += tok[i] + ' '
            atom['note'] = note.strip()
            atoms.append(atom)
        return atoms
        
# --------------------------------------

class Drude:
    """specification of drude oscillator types"""

    def __init__(self, drudefile):
        self.types = []
        with open(drudefile, "r") as f:
            for line in f:
                line = line.strip()
                if line.startswith('#') or len(line) == 0:
                    continue
                tok = line.split()
                drude = {}
                drude['type'] = tok[0]
                drude['dm'] = float(tok[1])
                drude['dq'] = float(tok[2])
                drude['k'] = float(tok[3])
                drude['alpha'] = float(tok[4])
                drude['thole'] = float(tok[5])
                self.types.append(drude)

# --------------------------------------

def atomline(at):
    return "{0:7d} {1:7d} {2:4d} {3:6.3f} {4:13.6e} {5:13.6e} {6:13.6e} "\
           "{7}\n".format(at['n'], at['mol'], at['id'], at['q'],
                          at['x'], at['y'], at['z'], at['note'])

def massline(att):
    return "{0:4d} {1:8.3f}  # {2}\n".format(att['id'], att['m'], att['type'])

def bdtline(bdt):
    return "{0:4d} {1:12.6f} {2:12.6f}  {3}\n".format(bdt['id'], bdt['k'],
                                                      bdt['r0'], bdt['note'])

def bondline(bd):
    return "{0:7d} {1:4d} {2:7d} {3:7d}  # {4}\n".format(bd['n'], bd['id'],
                                            bd['i'], bd['j'], bd['note'])

def drudeline(at):
    return "{0:7d} {1:1d} {2:7d}\n".format(at['n'], at['dflag'], at['dd'])

# --------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description = 'Add Drude dipole polarization to LAMMPS data file')
    parser.add_argument('-d', '--drudeff', default = 'drude.ff',
                        help = 'Drude parameter file (default: drude.ff)')
    parser.add_argument('infile', help = 'input LAMMPS data file')
    parser.add_argument('outfile', help = 'output LAMMPS data file')
    args = parser.parse_args()
    
    data = Data(args.infile)
    drude = Drude(args.drudeff)

    natom = data.headers['atoms']
    nbond = data.headers['bonds']
    nattype = data.headers['atom types']
    nbdtype = data.headers['bond types']

    # create new atom types (IDs) for Drude particles
    atomtypes = data.atomtypes()
    newids = []
    for att in atomtypes:
        for ddt in drude.types:
            if ddt['type'] == att['type']:
                nattype += 1
                newid = {}
                newid['id'] = ddt['id'] = nattype
                newid['m'] = ddt['dm']
                att['m'] -= ddt['dm']
                newid['type'] = att['type'] + ' DD'
                att['type'] += ' DC'
                ddt['type'] += ' DC'
                newids.append(newid)
                break

    data.headers['atom types'] += len(newids)
    for att in newids:
        data.sections['Masses'].append(massline(att))

    # create new bond types for core-Drude bonds
    bondtypes = data.bondtypes()
    newbdtypes = []
    for att in atomtypes:
        for ddt in drude.types:
            if ddt['type'] == att['type']:
                nbdtype += 1
                newbdtype = {}
                newbdtype['id'] = ddt['bdid'] = nbdtype
                newbdtype['k'] = ddt['k']
                newbdtype['r0'] = 0.0
                newbdtype['note'] = '# ' + ddt['type'] + '-DD'
                newbdtypes.append(newbdtype)
                break

    data.headers['bond types'] += len(newbdtypes)
    for bdt in newbdtypes:
        data.sections['Bond Coeffs'].append(bdtline(bdt))

    # create new atoms for Drude particles and bonds with their cores
    atoms = data.atoms()    
    newatoms = []
    newbonds = []
    for atom in atoms:
        atom['dflag'] = 0                 # 1: core, 2: drude, 0: other
        atom['dd'] = 0                    # partner drude or core
        for att in atomtypes:
            if att['id'] == atom['id']:
                break
        for ddt in drude.types:
            if ddt['type'] == att['type']:
                natom += 1
                newatom = deepcopy(atom)
                newatom['n'] = natom
                newatom['id'] = ddt['id']
                newatom['q'] = ddt['dq']
                newatom['note'] = atom['note'] + ' DD'
                newatom['dflag'] = 2
                newatom['dd'] = atom['n']
                newatoms.append(newatom)
                atom['q'] -= ddt['dq']
                atom['dflag'] = 1
                atom['dd'] = natom
                atom['note'] += ' DC'
                                
                nbond += 1
                newbond = {}
                newbond['n'] = nbond
                newbond['id'] = ddt['bdid']
                newbond['i'] = atom['n']
                newbond['j'] = newatom['n']
                newbond['note'] = ddt['type'] + '-DD'
                newbonds.append(newbond) 
                break
            
    data.headers['atoms'] += len(newatoms)
    data.headers['bonds'] += len(newbonds)
    data.sections['Atoms'] = []
    for atom in atoms + newatoms:
        data.sections['Atoms'].append(atomline(atom))        
    for bond in newbonds:                  
        data.sections['Bonds'].append(bondline(bond))

    data.sections['Drudes'] = []
    for atom in atoms + newatoms:
        data.sections['Drudes'].append(drudeline(atom))
        
    data.write(args.outfile)
        
                                    
if __name__ == '__main__':
    main()
