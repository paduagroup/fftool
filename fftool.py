#!/usr/bin/env python
# fftool.py - generate force field parameters for molecular system
# Agilio Padua <agilio.padua@univ-bpclermont.fr>, version 2013/11/29
# http://tim.univ-bpclermont.fr/apadua

# Copyright (C) 2013 Agilio A.H. Padua
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import sys, argparse, math

kcal = 4.184                            # kJ


# --------------------------------------

atomic_weights = {'H':    1.008,
                  'Li':   6.941,
                  'B':   10.811,
                  'C':   12.011,
                  'N':   14.006,
                  'O':   15.999,
                  'F':   18.998,
                  'Ne':  20.180,
                  'Na':  22.990,
                  'Mg':  24.305,
                  'Al':  26.982,
                  'Si':  28.086,
                  'P':   30.974,
                  'S':   32.065,
                  'Cl':  35.453,
                  'Ar':  39.948,
                  'K':   39.098,
                  'Ca':  40.078,
                  'Br':  79.904,
                  'I':  126.904}

def atomic_weight(name):
    if name[:2] in atomic_weights:
        return atomic_weights[name[:2]]
    elif name[0] in atomic_weights:
        return atomic_weights[name[0]]
    else:
        print 'warning: unknown atomic weight for %s atom' % (name)
        return 0.0

def atomic_symbol(name):
    if name[:2] in atomic_weights:
        return name[:2]
    elif name[0] in atomic_weights:
        return name[0]
    else:
        print 'warning: unknown symbol for %s atom' % (name)
        return ''


# --------------------------------------

class vector:
    '''minimal 3D vector'''

    def __init__(self, x = 0.0, y = 0.0, z = 0.0):
        if isinstance(x, tuple) or isinstance(x, list):
            self.x, self.y, self.z = x
        else:
            self.x = x
            self.y = y
            self.z = z

    def __getitem__(self, index):
        if index == 0:
            return self.x
        elif index == 1:
            return self.y
        elif index == 2:
            return self.z
        else:
            raise IndexError('vector index out of range')

    def __abs__(self):
        return math.sqrt(self.x**2 + self.y**2 + self.z**2)
    
    def __add__(self, other):
        if isinstance(other, vector):
            return vector(self.x + other.x, self.y + other.y, self.z + other.z)
        else:
            raise TypeError('wrong type in vector addition')

    def __sub__(self, other):
        if isinstance(other, vector):
            return vector(self.x - other.x, self.y - other.y, self.z - other.z)
        else:
            raise TypeError('wrong type in vector subtraction')

    def __mul__(self, other): 
        if isinstance(other, vector): # dot product
            return self.x*other.x + self.y*other.y + self.z*other.z
        else:
            return vector(self.x*other, self.y*other, self.z*other)

    def __div__(self, other):
        return vector(self.x/other, self.y/other, self.z/other)

    def __neg__(self):
        return vector(-self.x, -self.y, -self.z)

    def __str__(self):
        return '( ' + ', '.join([str(val) for val in (self.x, self.y, self.z)]) + ' )'

    def __repr__(self):
        return str(self) + ' instance at 0x' + str(hex(id(self))[2:].upper())    

    def cross(self, other):
        return vector(self.y * other.z - self.z * other.y,  
                      self.z * other.x - self.x * other.z,  
                      self.x * other.y - self.y * other.x)

    def unit(self):
        return self / abs(self)


# --------------------------------------

class zmat:
    '''z-matrix representing a molecule, read from file'''

    def __init__(self, filename, nmols = 1):
        self.nmols = nmols
        self.zatom = []
        self.connect = []
        self.improper = []
        
        # TODO: add code to handle z-matrix with tokens
        
        with open(filename, 'r') as f:

            # read molecule name
            line = f.readline().strip() 
            while line.startswith('#') or line == '':
                line = f.readline().strip()
            self.name = line

            #read z-matrix
            line = f.readline().strip()
            while line.startswith('#') or line == '':
                line = f.readline().strip()
            
            tok = line.split()
            if len(tok) > 1:   # there can be line numbers
                shift = 1
            else:
                shift = 0

            variables = False
            while line != '' and not line.lower().startswith('var'):
                tok = line.split()
                name = tok[shift]
                ir = ia = id = 0
                r = a = d = 0.0
                rvar = avar = dvar = ''
                if (len(tok) - shift) > 1:
                    ir = int(tok[shift+1])
                    if tok[shift+2][0].isalpha():
                        rvar = tok[shift+2]
                        variables = True
                    else:
                        r = float(tok[shift+2])
                    if (len(tok) - shift) > 3:
                        ia = int(tok[shift+3])
                        if tok[shift+4][0].isalpha():
                            avar = tok[shift+4]
                            variables = True
                        else:
                            a = float(tok[shift+4])
                        if (len(tok) - shift) > 5:
                            id = int(tok[shift+5])
                            if tok[shift+6][0].isalpha():
                                dvar = tok[shift+6]
                                variables = True
                            else:
                                d = float(tok[shift+6])
                zatom = {'name': name,
                        'ir': ir, 'rvar': rvar, 'r': r,
                        'ia': ia, 'avar': avar, 'a': a,
                        'id': id, 'dvar': dvar, 'd': d}
                self.zatom.append(zatom)
                line = f.readline().strip()
                
            # read variables
            if variables:
                if line.lower().startswith('var') or line == '':
                    line = f.readline().strip()
                while line != '':
                    tok = line.split('=')
                    key = tok[0].strip()
                    val = float(tok[1])
                    for rec in self.zatom:
                        if rec['rvar'] == key:
                            rec['r'] = val
                        if rec['avar'] == key:
                            rec['a'] = val
                        if rec['dvar'] == key:
                            rec['d'] = val
                    line = f.readline().strip()
                        
            # read connects
            while line.startswith('#') or line == '':
                line = f.readline().strip()
            while line.startswith('connect'):
                tok = line.split()
                atomi = int(tok[1])
                atomj = int(tok[2])
                self.connect.append([atomi, atomj])
                line = f.readline().strip()

            # read improper dihedrals
            while line.startswith('#') or line == '':
                line = f.readline().strip()            
            while line.startswith('improper'):
                tok = line.split()
                atomi = int(tok[1])
                atomj = int(tok[2])
                atomk = int(tok[3])
                atoml = int(tok[4])
                self.improper.append([atomi, atomj, atomk, atoml])
                line = f.readline().strip()
                
            # read force field file
            while line.startswith('#') or line == '':
                line = f.readline().strip()
            self.ff = line.split()[0]
            
    def show(self):
        print self.name
        print
        i = 0
        for rec in self.zatom:
            i += 1
            if rec['ir'] == 0:
                print '%-3d %-5s' % (i, rec['name'])
            elif rec['ia'] == 0:
                print '%-3d %-5s %3d %6.3f' % (i, rec['name'], rec['ir'], rec['r'])
            elif rec['id'] == 0:
                print '%-3d %-5s %3d %6.3f %3d %6.1f' % \
                    (i, rec['name'], rec['ir'], rec['r'], rec['ia'], rec['a'])
            else:
                print '%-3d %-5s %3d %6.3f %3d %6.1f %3d %6.1f' % \
                    (i, rec['name'], rec['ir'], rec['r'], rec['ia'], rec['a'],
                     rec['id'], rec['d'])
        print
        if len(self.connect) > 0:
            print 'connects'
            for c in self.connect:
                print '%3d (%5s) -- %3d (%5s)' % \
                    (c[0], self.zatom[c[0]-1]['name'],
                     c[1], self.zatom[c[1]-1]['name'])
            print
        print self.ff


# --------------------------------------
 
class atom:
    '''atom in a molecule or in a force field'''

    def __init__(self, name, m = 0.0):
        self.name = name
        if (m == 0.0):
            self.m = atomic_weight(self.name)
        else:
            self.m = m
        self.ityp = -1                    # index of atom type for this atom
        self.x = 0.0
        self.y = 0.0
        self.z = 0.0

    def __str__(self):
        if hasattr(self, 'type'):
            return 'atom %-5s %-2s  m = %7.3f  q = %+6.3f  %s %s' % \
              (self.name, self.type, self.m, self.q, self.pot, str(self.par))
        else:
            return 'atom %-5s  m = %7.3f' % (self.name, self.m)

    def setpar(self, attp, q, pot, par):
        self.type = attp
        self.q = q
        self.pot = pot
        self.par = par


class bond:
    '''covalent bond in a molecule or in a force field'''

    def __init__(self, i = -1, j = -1, r = 0.0):
        self.i = i
        self.j = j
        self.r = r
        self.ityp = -1

    def __str__(self):
        if hasattr(self, 'name'):
            if self.i != -1:
                return 'bond %5d %5d  %s  %s %s' % \
                  (self.i + 1, self.j + 1, self.name, self.pot, str(self.par))
            else:
                return 'bond %s  %s %s' % (self.name, self.pot, str(self.par))                
        else:
            return 'bond %5d %5d' % (self.i + 1, self.j + 1)

    def setpar(self, iatp, jatp, pot, par):
        self.name = '%s-%s' % (iatp, jatp)
        self.iatp = iatp
        self.jatp = jatp
        self.pot = pot
        self.par = par


class angle:
    '''valence angle'''

    def __init__(self, i = -1, j = -1, k = -1, theta = 0.0):
        self.i = i
        self.j = j
        self.k = k
        self.theta = theta
        self.ityp = -1

    def __str__(self):
        if hasattr(self, 'name'):
            if self.i != -1:
                return 'angle %5d %5d %5d  %s  %s %s' % \
                  (self.i + 1, self.j + 1, self.k + 1,
                   self.name, self.pot, str(self.par))
            else:
                return 'angle %s  %s %s' % (self.name, self.pot, str(self.par))                
        else:
            return 'angle %5d %5d %5d' % (self.i + 1, self.j + 1, self.k + 1)

    def setpar(self, iatp, jatp, katp, pot, par):
        self.name = '%s-%s-%s' % (iatp, jatp, katp)
        self.iatp = iatp
        self.jatp = jatp
        self.katp = katp
        self.pot = pot
        self.par = par


class dihed:
    '''dihedral angle (torsion)'''

    def __init__(self, i = -1, j = -1, k = -1, l = -1, phi = 0.0):
        self.i = i
        self.j = j
        self.k = k
        self.l = l
        self.phi = phi
        self.ityp = -1

    def __str__(self):
        if hasattr(self, 'name'):
            if self.i != -1:
                return 'dihedral %5d %5d %5d %5d  %s  %s %s' % \
                  (self.i + 1, self.j + 1, self.k + 1, self.l + 1, \
                   self.name, self.pot, str(self.par))
            else:
                return 'dihedral %s  %s %s' % \
                  (self.name, self.pot, str(self.par))                
        else:
            return 'dihedral %5d %5d %5d %5d' % \
              (self.i + 1, self.j + 1, self.k + 1, self.l + 1)

    def setpar(self, iatp, jatp, katp, latp, pot, par):
        self.name = '%s-%s-%s-%s' % (iatp, jatp, katp, latp)
        self.iatp = iatp
        self.jatp = jatp
        self.katp = katp
        self.latp = latp
        self.pot = pot
        self.par = par


class dimpr:
    '''dihedral angle (improper)'''

    def __init__(self, i = -1, j = -1, k = -1, l = -1, phi = 0.0):
        self.i = i
        self.j = j
        self.k = k
        self.l = l
        self.phi = phi
        self.ityp = -1

    def __str__(self):
        if hasattr(self, 'name'):
            if self.i != -1:
                return 'improper %5d %5d %5d %5d  %s  %s %s' % \
                  (self.i + 1, self.j + 1, self.k + 1, self.l + 1, \
                   self.name, self.pot, str(self.par))
            else:
                return 'improper %s  %s %s' % \
                  (self.name, self.pot, str(self.par))                
        else:
            return 'improper %5d %5d %5d %5d' % \
              (self.i + 1, self.j + 1, self.k + 1, self.l + 1)

    def setpar(self, iatp, jatp, katp, latp, pot, par):
        self.name = '%s-%s-%s-%s' % (iatp, jatp, katp, latp)
        self.iatp = iatp
        self.jatp = jatp
        self.katp = katp
        self.latp = latp
        self.pot = pot
        self.par = par


class mol:
    '''molecule'''
    
    def zmat2cart(self, z):
        natoms = len(self.atom)    
        if natoms != len(z.zatom):
            print 'error in mol.zmat2cart(): different values of natoms'
            sys.exit(2)

        if natoms == 0:
            return self

        # first atom at origin
        self.atom[0].x = 0.0
        self.atom[0].y = 0.0
        self.atom[0].z = 0.0
        if natoms == 1:
            return self

        # second atom at distance r from first along xx
        self.atom[1].x = z.zatom[1]['r']
        self.atom[1].y = 0.0
        self.atom[1].z = 0.0
        if natoms == 2:
            return self

        # third atom at distance r from ir forms angle a 3-ir-ia in plane xy
        r = z.zatom[2]['r']
        ir = z.zatom[2]['ir'] - 1
        ang = z.zatom[2]['a'] * math.pi / 180.0
        ia = z.zatom[2]['ia'] - 1

        # for this construction, the new atom is at point (x, y), atom
        # ir is at point (xr, yr) and atom ia is at point (xa, ya).
        # Theta is the angle between the vector joining ir to ia and
        # the x-axis, a' (= theta - a) is is the angle between r and
        # the x-axis. x = xa + r cos a', y = ya + r sin a'.  From the
        # dot product of a unitary vector along x with the vector from
        # ir to ia, theta can be calculated: cos theta = (xa - xr) /
        # sqrt((xa - xr)^2 + (ya - yr)^2).  If atom ia is in third or
        # forth quadrant relative to atom ir, ya - yr < 0, then theta
        # = 2 pi - theta. */
        delx = self.atom[ia].x - self.atom[ir].x
        dely = self.atom[ia].y - self.atom[ir].y
        theta = math.acos(delx / math.sqrt(delx*delx + dely*dely))
        if dely < 0.0:
            theta = 2 * math.pi - theta
        ang = theta - ang
        self.atom[2].x = self.atom[ir].x + r * math.cos(ang)
        self.atom[2].y = self.atom[ir].y + r * math.sin(ang)
        self.atom[2].z = 0.0
        if natoms == 3:
            return self
        
        # nth atom at distance r from atom ir forms angle a at 3-ir-ia
        # and dihedral angle between planes 3-ir-ia and ir-ia-id
        for i in range(3, natoms):
            r = z.zatom[i]['r']
            ir = z.zatom[i]['ir'] - 1
            ang = z.zatom[i]['a'] * math.pi / 180.0
            ia = z.zatom[i]['ia'] - 1
            dih = z.zatom[i]['d'] * math.pi / 180.0
            id = z.zatom[i]['id'] - 1

            # for this construction the new atom is at point A, atom ir is
            # at B, atom ia at C and atom id at D.  Point a is the
            # projection of A onto the plane BCD.  Point b is the
            # projection of A along the direction BC (the line defining
            # the dihedral angle between planes ABC and BCD). n = CD x BC
            # / |CD x BC| is the unit vector normal to the plane BCD. m =
            # BC x n / |BC x n| is the unit vector on the plane BCD normal
            # to the direction BC.
            #                               
            #                               .'A
            #                 ------------.' /.-----------------
            #                /           b /  .               /
            #               /           ./    .              /
            #              /           B......a      ^      /
            #             /           /              |n    /
            #            /           /                    /
            #           /           C                    /
            #          /             \                  /
            #         /               \                /
            #        /plane BCD        D              /
            #       ----------------------------------
            #
            #                    A              C------B...b
            #                   /.             /        .  .
            #                  / .            /    |m    . .
            #                 /  .           /     V      ..
            #         C------B...b          D              a
            #
            
            BA = r
            vB = vector(self.atom[ir].x, self.atom[ir].y, self.atom[ir].z)
            vC = vector(self.atom[ia].x, self.atom[ia].y, self.atom[ia].z)
            vD = vector(self.atom[id].x, self.atom[id].y, self.atom[id].z)

            vBC = vC - vB
            vCD = vD - vC
            
            BC = abs(vBC)
            bB = BA * math.cos(ang)
            bA = BA * math.sin(ang)
            aA = bA * math.sin(dih)
            ba = bA * math.cos(dih)

            vb = vC - vBC * ((BC - bB) / BC)
            vn = (vCD.cross(vBC)).unit()
            vm = (vBC.cross(vn)).unit()
            va = vb + vm * ba
            vA = va + vn * aA

            self.atom[i].x = vA.x
            self.atom[i].y = vA.y
            self.atom[i].z = vA.z

        return self
    
    def fromzmat(self, z):
        self.name = z.name

        for zat in z.zatom:
            self.atom.append(atom(zat['name']))
            self.m += atomic_weight(zat['name'])
        i = 1
        while i < len(z.zatom):    
            self.bond.append(bond(i, z.zatom[i]['ir'] - 1))
            i += 1
        for cn in z.connect:
            self.bond.append(bond(cn[0] - 1, cn[1] - 1))
            
        self.zmat2cart(z)

        natoms = len(self.atom)
        nbonds = len(self.bond)

        # identify valence angles
        i = 0
        while i < natoms:        # find neighbour atoms to each atom i
            nb = 0
            neib = []
            for bd in self.bond:          
                if i == bd.i:
                    neib.append(bd.j)
                    nb += 1
                elif i == bd.j:
                    neib.append(bd.i)
                    nb += 1
            k = 0
            while k < nb - 1:
                l = k + 1
                while l < nb:
                    self.angle.append(angle(neib[k], i, neib[l]))
                    l += 1
                k += 1
            i += 1

        # identify dihedral angles
        k = 0
        while k < nbonds:       # find bonds around non-terminal bonds
            l = 0
            while l < nbonds:
                if k == l:
                    l += 1
                    continue
                if self.bond[k].i == self.bond[l].i:
                    j = 0
                    while j < nbonds:
                        if j == k or j == l:
                            j += 1
                            continue
                        if self.bond[k].j == self.bond[j].i:
                            self.dihed.append(dihed(self.bond[l].j, self.bond[k].i,
                                                    self.bond[k].j, self.bond[j].j))
                        elif self.bond[k].j == self.bond[j].j:
                            self.dihed.append(dihed(self.bond[l].j, self.bond[k].i,
                                                    self.bond[k].j, self.bond[j].i))
                        j += 1
                elif self.bond[k].i == self.bond[l].j:
                    j = 0
                    while j < nbonds:
                        if j == k or j == l:
                            j += 1
                            continue
                        if self.bond[k].j == self.bond[j].i:
                            self.dihed.append(dihed(self.bond[l].i, self.bond[k].i,
                                                    self.bond[k].j, self.bond[j].j))
                        elif self.bond[k].j == self.bond[j].j:
                            self.dihed.append(dihed(self.bond[l].i, self.bond[k].i,
                                                    self.bond[k].j, self.bond[j].i))
                        j += 1
                l += 1
            k += 1

        # add improper dihedrals
        for di in z.improper:                 
            self.dimpr.append(dihed(di[0] - 1, di[1] - 1, di[2] - 1, di[3] - 1))

        self.ff = z.ff
        return self

    def __init__(self, name, nmols = 1):
        self.atom = []
        self.bond = []
        self.angle = []
        self.dihed = []
        self.dimpr = []
        self.m = 0
        
        try:                              # read from zmat file
            with open(name, 'r'):
                self.filename = name
            self.fromzmat(zmat(name))                
        except IOError:
            self.filename = ''
            self.name = name

        self.nmols = nmols
        
    def __str__(self):
        return 'molecule %s  %d atoms  m = %8.4f' % \
            (self.name, len(self.atom), self.m)
            
    def charge(self):
        q = 0.0
        for at in self.atom:
            q += at.q
        return q

    def show(self):
        print '%s: %d molecules' % (self.name, self.nmols)
        print '%d atoms' % len(self.atom)
        for at in self.atom:
            print at
        print '%d bonds' % len(self.bond)
        for bd in self.bond:
            print bd
        print '%d angles' % len(self.angle)
        for an in self.angle:
            print an
        print '%d dihedrals' % len(self.dihed)
        for dh in self.dihed:
            print dh
        print '%d improper' % len(self.dimpr)
        for di in self.dimpr:
            print di

    def showxyz(self, symbol = False):
        print len(self.atom)
        print self.name
        for a in self.atom:
            if symbol:
                atname = atomic_symbol(a.name)
            else:
                atname = a.name
            print '%-5s %15.6f %15.6f %15.6f' % (atname, a.x, a.y, a.z)

    def writexyz(self, symbol = False):
        outfile = (self.filename).rsplit('.', 1)[0] + '.xyz'
        with open(outfile, 'w') as f:
            f.write(str(len(self.atom)) + '\n')
            f.write(self.name + '\n')
            for a in self.atom:
                if symbol:
                    atname = atomic_symbol(a.name)
                else:
                    atname = a.name
                f.write('%-5s %15.6f %15.6f %15.6f\n' % (atname, a.x, a.y, a.z))

            
# --------------------------------------

class forcefield:
    '''force field parameter database'''

    def __init__(self, filename):
        self.filename = filename
        self.atom = []
        self.bond = []
        self.angle = []
        self.dihed = []
        self.dimpr = []
        
        with open(self.filename, 'r') as f:
            i = ib = ia = ih = im = 0
            for line in f:
                if line.startswith('#') or line.strip() == '':
                    continue
                
                if line.lower().startswith('atom'):
                    section = 'atoms'
                    continue
                elif line.lower().startswith('bond'):
                    section = 'bonds'
                    continue
                elif line.lower().startswith('angl'):
                    section = 'angles'
                    continue
                elif line.lower().startswith('dihe'):
                    section = 'dihedrals'
                    continue
                elif line.lower().startswith('impro'):
                    section = 'improper'
                    continue

                tok = line.strip().split()

                if section == 'atoms':
                    name = tok[0]
                    attp = tok[1]
                    m = float(tok[2])
                    q = float(tok[3])
                    pot = tok[4]
                    par = [float(p) for p in tok[5:]]
                    self.atom.append(atom(name, m))
                    self.atom[i].setpar(attp, q, pot, par)
                    i += 1

                elif section == 'bonds':
                    iatp = tok[0]
                    jatp = tok[1]
                    pot = tok[2]
                    par = [float(p) for p in tok[3:]]
                    self.bond.append(bond())
                    self.bond[ib].setpar(iatp, jatp, pot, par)
                    ib += 1

                elif section == 'angles':
                    iatp = tok[0]
                    jatp = tok[1]
                    katp = tok[2]
                    pot = tok[3]
                    par = [float(p) for p in tok[4:]]
                    self.angle.append(angle())
                    self.angle[ia].setpar(iatp, jatp, katp, pot, par)
                    ia += 1

                elif section == 'dihedrals':
                    iatp = tok[0]
                    jatp = tok[1]
                    katp = tok[2]
                    latp = tok[3]
                    pot = tok[4]
                    par = [float(p) for p in tok[5:]]
                    self.dihed.append(dihed())
                    self.dihed[ih].setpar(iatp, jatp, katp, latp, pot, par)
                    ih += 1

                elif section == 'improper':
                    iatp = tok[0]
                    jatp = tok[1]
                    katp = tok[2]
                    latp = tok[3]
                    pot = tok[4]
                    par = [float(p) for p in tok[5:]]
                    self.dimpr.append(dimpr())
                    self.dimpr[im].setpar(iatp, jatp, katp, latp, pot, par)
                    im += 1

    def show(self):
        for at in self.atom:
            print at
        for bd in self.bond:
            print bd
        for an in self.angle:
            print an
        for dh in self.dihed:
            print dh
        for di in self.dimpr:
            print di


# --------------------------------------


class vdw:
    '''van der Waals interaction'''        
    
    def __init__(self, iat, jat, mix = 'g'):
        self.i = iat.name
        self.j = jat.name
        self.ityp = iat.ityp
        self.jtyp = jat.ityp
        
        if iat.pot != jat.pot:
            print 'error in vdw object: incompatible potential types'
            sys.exit(2)

        self.pot = iat.pot

        if len(iat.par) != len(jat.par):
            print 'error in vdw object: different lengths in parameter lists'
            sys.exit(2)

        if self.pot == 'lj':
            if iat.name == jat.name:
                self.par = iat.par
            else:
                self.par = [0.0, 0.0]
                if mix == 'g':
                    self.par[0] = math.sqrt(iat.par[0] * jat.par[0])
                else:
                    self.par[0] = (iat.par[0] + jat.par[0]) / 2.
                self.par[1] = math.sqrt(iat.par[1] * jat.par[1])
                
    def __str__(self):
        return 'vdw %2s %2s  %s %s' % (self.i, self.j, self.pot, str(self.par))


# --------------------------------------
    
def build_type_list(term, termtype):
    # build a list of atom or bonded term types based on the name attribute
    # term is an input list of atoms or bonded terms (bonds, angles, dihedrals)
    # termtype is a list containing the different types found
    for a in term:
        found = False
        for b in termtype:
            if a.name == b.name:
                found = True
        if not found:
            termtype.append(a)

def assign_type_index(term, termtype):
    # assign index numbers to the ityp attribute in atoms or bonded terms
    ntypes = len(termtype)
    for a in term:
        i = 0
        while i < ntypes:
            if a.name == termtype[i].name:
                a.ityp = termtype[i].ityp = i
                break       
            i += 1


class system:
    '''Molecular system to be simulated'''
                
    def __init__(self, molecules, mix):
        self.mol = molecules                       # list of molecules
        self.attype = []                           # atom types
        self.bdtype = []                           # bond types
        self.antype = []                           # angle types
        self.dhtype = []                           # dihedral types
        self.ditype = []                           # improper types
        self.vdw = []

        # set force field parameters
        atomerror = error = False
        for m in self.mol:
            ff = forcefield(m.ff)

            # identify atom types and set parameters
            for at in m.atom:
                found = False
                for ffat in ff.atom:     
                    if at.name == ffat.name:
                        at.setpar(ffat.type, ffat.q, ffat.pot, ffat.par)
                        found = True
                if not found:
                    print 'error in molecule %s: no parameters for atom %s' % \
                      (m.name, at.name)
                    atomerror = True
                    continue    # no point in looking for bonded terms
            
            # identify bonded terms and set parameters
            for bd in m.bond:
                bd.name = '%s-%s' % (m.atom[bd.i].type, m.atom[bd.j].type)
                found = False
                for ffbd in ff.bond:
                    namestr = '%s-%s' % (ffbd.iatp, ffbd.jatp)
                    namerev = '%s-%s' % (ffbd.jatp, ffbd.iatp)
                    if bd.name == namestr or bd.name == namerev: 
                        bd.setpar(ffbd.iatp, ffbd.jatp, ffbd.pot, ffbd.par)
                        found = True
                if not found:
                    print 'error in molecule %s: no parameters for bond %s' % \
                      (m.name, bd.name)
                    error = True
                                        
            for an in m.angle:
                an.name = '%s-%s-%s' % \
                  (m.atom[an.i].type, m.atom[an.j].type, m.atom[an.k].type)
                found = False
                for ffan in ff.angle:
                    namestr = '%s-%s-%s' % (ffan.iatp, ffan.jatp, ffan.katp)
                    namerev = '%s-%s-%s' % (ffan.katp, ffan.jatp, ffan.iatp)
                    if an.name == namestr or an.name == namerev:
                        an.setpar(ffan.iatp, ffan.jatp, ffan.katp,
                                  ffan.pot, ffan.par)
                        found = True
                if not found:
                    print 'error in molecule %s: no parameters for angle %s' % \
                      (m.name, an.name)
                    error = True
                        
            for dh in m.dihed:
                dh.name = '%s-%s-%s-%s' % (m.atom[dh.i].type, m.atom[dh.j].type,
                                           m.atom[dh.k].type, m.atom[dh.l].type)
                found = False
                for ffdh in ff.dihed:
                    namestr = '%s-%s-%s-%s' % \
                      (ffdh.iatp, ffdh.jatp, ffdh.katp, ffdh.latp)
                    namerev = '%s-%s-%s-%s' % \
                      (ffdh.latp, ffdh.katp, ffdh.jatp, ffdh.iatp)
                    if dh.name == namestr or dh.name == namerev:
                        dh.setpar(ffdh.iatp, ffdh.jatp, ffdh.katp, ffdh.latp,
                                  ffdh.pot, ffdh.par)
                        found = True
                if not found:
                    print 'error in molecule %s: no parameters for dihedral %s' % \
                      (m.name, dh.name)
                    error = True

            for di in m.dimpr:
                di.name = '%s-%s-%s-%s' % (m.atom[di.i].type, m.atom[di.j].type,
                                           m.atom[di.k].type, m.atom[di.l].type)
                found = False
                for ffdi in ff.dimpr:
                    namestr = '%s-%s-%s-%s' % \
                      (ffdi.iatp, ffdi.jatp, ffdi.katp, ffdi.latp)
                    namerev = '%s-%s-%s-%s' % \
                      (ffdi.latp, ffdi.katp, ffdi.jatp, ffdi.iatp)
                    if di.name == namestr or di.name == namerev:
                        di.setpar(ffdi.iatp, ffdi.jatp, ffdi.katp, ffdi.latp,
                                  ffdi.pot, ffdi.par)
                        found = True
                if not found:
                    print 'error in molecule %s: no parameters for improper %s' % \
                      (m.name, di.name)
                    error = True

        if atomerror or error:
            sys.exit(2)
                    
        # at this point force field infomation was read for all molecules
        
        # build lists of different atom and bonded term types in the system
        for m in self.mol:
            build_type_list(m.atom, self.attype)
            build_type_list(m.bond, self.bdtype)
            build_type_list(m.angle, self.antype)
            build_type_list(m.dihed, self.dhtype)
            build_type_list(m.dimpr, self.ditype)

        nattypes = len(self.attype)
        nbdtypes = len(self.bdtype)
        nantypes = len(self.antype)
        ndhtypes = len(self.dhtype)
        nditypes = len(self.ditype)

        # assign the type index for all atoms and bonded terms in the system
        for m in self.mol:
            assign_type_index(m.atom, self.attype)
            assign_type_index(m.bond, self.bdtype)
            assign_type_index(m.angle, self.antype)
            assign_type_index(m.dihed, self.dhtype)
            assign_type_index(m.dimpr, self.ditype)

        # set non-bonded parameters for all i-j pairs
        i = 0
        while i < nattypes:
            j = i
            while j < nattypes:
                self.vdw.append(vdw(self.attype[i], self.attype[j], mix))
                j += 1
            i += 1

    def show(self):
        for m in self.mol:
            print '%s  %d molecules force field %s' % (m.name, m.nmols, m.ff)
            for at in m.atom:
                print at
            for bd in m.bond:
                print bd
            for an in m.angle:
                print an
            for dh in m.dihed:
                print dh
            for di in m.dimpr:
                print di
        for nb in self.vdw:
            print nb

    def writepackmol(self, boxlen):
        with open('pack.inp', 'w') as f:
            f.write('# created by fftool\n')
            f.write('tolerance 3.0\n')
            f.write('filetype xyz\n')
            f.write('output simbox.xyz\n')
            for m in self.mol:
                xyzfile = (m.filename).rsplit('.', 1)[0] + '.xyz'
                f.write('\nstructure %s\n' % xyzfile)
                f.write('  number %s\n' % m.nmols)
                f.write('  inside box %.1f %.1f %.1f %.1f %.1f %.1f\n' % \
                        (-boxlen/2., -boxlen/2., -boxlen/2.,
                         boxlen/2., boxlen/2., boxlen/2.))
                f.write('end structure\n')
        return boxlen
                
    def writelmp(self, boxlen, mix = 'g', allpairs = False):
        try:
            with open('simbox.xyz', 'r') as fx:
                pass
        except IOError:
            print 'warning: cannot open simbox.xyz, lammps files not created'
            return

        natoms = nbonds = nangles = ndiheds = 0
        for m in self.mol:
            natoms += m.nmols * len(m.atom)
            nbonds += m.nmols * len(m.bond)
            nangles += m.nmols * len(m.angle)
            ndiheds += m.nmols * (len(m.dihed) + len(m.dimpr))
        
        with open('in.lmp', 'w') as fi:
            fi.write('# created by fftool\n\n')
            fi.write('units real\n')
            fi.write('boundary p p p\n\n')

            fi.write('atom_style full\n')
            if nbonds > 0:
                fi.write('bond_style harmonic\n')
            if nangles > 0:
                fi.write('angle_style harmonic\n')
            if ndiheds > 0:
                fi.write('dihedral_style opls\n')
            fi.write('special_bonds lj/coul 0.0 0.0 0.5\n\n')

            fi.write('read_data data.lmp\n')
            fi.write('# read_restart restart.*.lmp\n')
            fi.write('# reset_timestep 0\n\n')

            fi.write('pair_style hybrid lj/cut/coul/long 10.0 10.0\n')
            if not allpairs:
                if (mix == 'g'):
                    fi.write('pair_modify mix geometric tail yes\n')
                else:
                    fi.write('pair_modify mix arithmetic tail yes\n')
                fi.write('kspace_style pppm 1.0e-4\n\n')
                for att in self.attype:
                    fi.write('pair_coeff %4d %4d  %s  %8.4f %8.4f  # %s %s\n' % \
                             (att.ityp + 1, att.ityp + 1, 'lj/cut/coul/long',
                              att.par[1] / kcal, att.par[0], att.name, att.name))
            else:
                fi.write('pair_modify tail yes\n')
                fi.write('kspace_style pppm 1.0e-4\n\n')
                for nb in self.vdw:
                    fi.write('pair_coeff %4d %4d  %s  %8.4f %8.4f  # %s %s\n' % \
                             (nb.ityp + 1, nb.jtyp + 1, 'lj/cut/coul/long',
                              nb.par[1] / kcal, nb.par[0], nb.i, nb.j))
            fi.write('\n')

            fi.write('variable nsteps equal 10000\n')
            fi.write('variable nprint equal ${nsteps}/100\n')
            fi.write('variable ndump equal ${nsteps}/100\n')
            fi.write('# variable nrestart equal ${nsteps}/10\n\n')

            fi.write('variable temp equal 300.0\n')
            fi.write('variable press equal 1.0\n\n')

            fi.write('neighbor 2.0 bin\n\n')

            fi.write('timestep 1.0\n\n')

            fi.write('velocity all create ${temp} 12345\n\n')            

            fi.write('thermo_style multi\n')
            fi.write('thermo ${nprint}\n\n')

            shakebd = shakean = False
            for bdt in self.bdtype:
                if bdt.pot == 'cons':
                    shakebd = True
            for ant in self.antype:
                if ant.pot == 'cons':
                    shakean = True
            if shakebd or shakean:
                fi.write('fix fSHAKE all shake 0.0001 20 ${nprint}')
                if shakebd:
                    fi.write(' b')
                    for bdt in self.bdtype:
                        if bdt.pot == 'cons':
                            fi.write(' %d' % (bdt.ityp + 1))
                if shakean:
                    fi.write(' a')
                    for ant in self.antype:
                        if ant.pot == 'cons':
                            fi.write(' %d' % (ant.ityp + 1))
                fi.write('\n\n')

            fi.write('fix fNPT all npt temp ${temp} ${temp} 100 '\
                     'iso ${press} ${press} 500\n\n')

            fi.write('# compute cRDF all rdf 100 1 1\n')
            fi.write('# fix fRDF all ave/time 20 100 ${nsteps} '\
                     'c_cRDF file rdf.lammps mode vector\n\n')
            
            fi.write('# compute cMSD all msd\n')
            fi.write('# fix fMSD all ave/time 1 1 ${ndump} '\
                     'c_cMSD[1] c_cMSD[2] c_cMSD[3] c_cMSD[4] file msd.lammps\n\n')

            fi.write('dump dCONF all custom ${ndump} dump.lammpstrj '\
                     'id mol type element x y z ix iy iz\n')
            fi.write('dump_modify dCONF element')
            for att in self.attype:
                fi.write(' %s' % atomic_symbol(att.name))
            fi.write('\n\n')

            fi.write('# restart ${nrestart} restart.*.lmp\n\n')

            fi.write('run ${nsteps}\n\n')

            fi.write('write_restart restart.*.lmp\n')
            fi.write('write_data data.*.lmp\n')

        with open('data.lmp', 'w') as fd:
            fd.write('created by fftool\n\n')
            fd.write('%d atoms\n' % natoms)
            if nbonds > 0:
                fd.write('%d bonds\n' % nbonds)
            if nangles > 0:
                fd.write('%d angles\n' % nangles)
            if ndiheds > 0:
                fd.write('%d dihedrals\n' % ndiheds)
            fd.write('\n')
                
            fd.write('%d atom types\n' % len(self.attype))
            if nbonds > 0:
                fd.write('%d bond types\n' % len(self.bdtype))
            if nangles > 0:
                fd.write('%d angle types\n' % len(self.antype))
            if ndiheds > 0:
                ndht = len(self.dhtype)     # needed later
                fd.write('%d dihedral types\n' % (ndht + len(self.ditype)))
            fd.write('\n')
            
            x = []
            y = []
            z = []
            with open('simbox.xyz', 'r') as fx:
                natoms = int(fx.readline().strip())
                title = fx.readline().strip()
                i = 0
                while i < natoms:
                    tok = fx.readline().strip().split()
                    x.append(float(tok[1]))
                    y.append(float(tok[2]))
                    z.append(float(tok[3]))
                    i += 1
        
            tol = 2.0
            fd.write('%f %f xlo xhi\n' % (-boxlen/2. - tol, boxlen/2. + tol))
            fd.write('%f %f ylo yhi\n' % (-boxlen/2. - tol, boxlen/2. + tol))
            fd.write('%f %f zlo zhi\n' % (-boxlen/2. - tol, boxlen/2. + tol))

            fd.write('\nMasses\n\n')
            for att in self.attype:
                fd.write('%4d %8.3f  # %s\n' % (att.ityp + 1, att.m, att.name))

            if nbonds > 0:
                fd.write('\nBond Coeffs\n\n')
                for bdt in self.bdtype:
                    fd.write('%4d %7.1f %6.3f  # %s\n' % \
                             (bdt.ityp + 1, bdt.par[1] / (2.0 * kcal),
                              bdt.par[0], bdt.name))

            if nangles > 0:
                fd.write('\nAngle Coeffs\n\n')
                for ant in self.antype:
                    fd.write('%4d %7.2f %7.2f  # %s\n' % \
                             (ant.ityp + 1, ant.par[1] / (2.0 * kcal),
                              ant.par[0], ant.name))

            if ndiheds > 0:
                fd.write('\nDihedral Coeffs\n\n')
                for dht in self.dhtype:
                    fd.write('%4d %9.4f %9.4f %9.4f %9.4f  # %s\n' % \
                             (dht.ityp + 1,
                              dht.par[0] / kcal, dht.par[1] / kcal,
                              dht.par[2] / kcal, dht.par[3] / kcal, dht.name))
                for dit in self.ditype:
                    fd.write('%4d %9.4f %9.4f %9.4f %9.4f  # %s\n' % \
                             (ndht + dit.ityp + 1,
                              dit.par[0] / kcal, dit.par[1] / kcal,
                              dit.par[2] / kcal, dit.par[3] / kcal, dit.name))

            fd.write('\nAtoms\n\n')
            i = 0
            nmol = 0
            for m in self.mol:
                im = 0
                while im < m.nmols:
                    for at in m.atom:
                        fd.write('%7d %7d %4d %6.3f %13.6e %13.6e %13.6e  '\
                                 '# %-6s %s\n' % \
                                 (i + 1, nmol + 1, at.ityp + 1, at.q, 
                                  x[i], y[i], z[i], at.name, m.name))
                        i += 1
                    im += 1
                    nmol += 1

            if nbonds > 0:
                fd.write('\nBonds\n\n')
                i = shift = 1
                for m in self.mol:
                    natoms = len(m.atom)
                    im = 0
                    while im < m.nmols:
                        for bd in m.bond:
                            fd.write('%7d %4d %7d %7d  # %s\n' % \
                                     (i, bd.ityp + 1, bd.i + shift, bd.j + shift,
                                      bd.name))
                            i += 1
                        shift += natoms
                        im += 1

            if nangles > 0:
                fd.write('\nAngles\n\n')
                i = shift = 1
                for m in self.mol:
                    natoms = len(m.atom)
                    im = 0
                    while im < m.nmols:
                        for an in m.angle:
                            fd.write('%7d %4d %7d %7d %7d  # %s\n' % \
                                     (i, an.ityp + 1, an.i + shift,
                                      an.j + shift, an.k + shift, an.name))
                            i += 1
                        shift += natoms
                        im += 1

            if ndiheds > 0:
                fd.write('\nDihedrals\n\n')
                i = shift = 1
                for m in self.mol:
                    natoms = len(m.atom)
                    im = 0
                    while im < m.nmols:
                        for dh in m.dihed:
                            fd.write('%7d %4d %7d %7d %7d %7d  # %s\n' % \
                                     (i, dh.ityp + 1, dh.i + shift, dh.j + shift,
                                      dh.k + shift, dh.l + shift, dh.name))
                            i += 1
                        for di in m.dimpr:
                            fd.write('%7d %4d %7d %7d %7d %7d  # %s\n' % \
                                     (i, ndht + di.ityp + 1, di.i + shift, di.j + shift,
                                      di.k + shift, di.l + shift, di.name))
                            i += 1
                        shift += natoms
                        im += 1
                    
            fd.write('\n')
                    

    def writedlp(self, boxlen, cos4 = False):
        with open('FIELD', 'w') as f:
            f.write('created by fftool\n')
            f.write('units kJ\n\n')
            
            f.write('molecular types %d\n' % (len(self.mol)))
            for mol in self.mol:
                f.write('%s\n' % mol.name)
                f.write('nummols %d\n' % mol.nmols)

                f.write('atoms %d\n' % len(mol.atom))
                for at in mol.atom:
                    f.write('%-5s %8.4f %6.3f 1  # %s\n' % \
                            (at.name, at.m, at.q, at.type))

                ncons = 0
                for bd in mol.bond: 
                    if bd.pot == 'cons':
                        ncons += 1
                f.write('constraints %d\n' % ncons)
                for bd in mol.bond:
                    if bd.pot == 'cons':
                        f.write('%4d %4d %6.3f  # %s\n' % \
                                (bd.i + 1, bd.j + 1, bd.par[0], bd.name))
                f.write('bonds %d\n' % (len(mol.bond) - ncons))
                for bd in mol.bond:
                    if bd.pot != 'cons':
                        f.write('%4s %4d %4d %7.1f %6.3f  # %s\n' % \
                                (bd.pot, bd.i + 1, bd.j + 1, bd.par[1], bd.par[0],
                                 bd.name))
                                                                  
                f.write('angles %d\n' % len(mol.angle))
                for an in mol.angle:
                    f.write('%4s %4d %4d %4d %7.2f %7.2f  # %s\n' % \
                            (an.pot, an.i + 1, an.j + 1, an.k + 1,
                             an.par[1], an.par[0], an.name))
                             
                f.write('dihedrals %d\n' % (len(mol.dihed) + len(mol.dimpr)))
                for dh in mol.dihed:
                    if cos4:
                        pot = 'cos4'
                        f.write('%4s %4d %4d %4d %4d %9.4f %9.4f %9.4f %9.4f'\
                                ' %6.3f %6.3f  # %s\n' % \
                                (pot, dh.i + 1, dh.j + 1, dh.k + 1, dh.l + 1,
                                 dh.par[0], dh.par[1], dh.par[2], dh.par[3],
                                 0.5, 0.5, dh.name))
                    else:
                        pot = 'cos3'
                        f.write('%4s %4d %4d %4d %4d %9.4f %9.4f %9.4f'\
                                ' %6.3f %6.3f  # %s\n' % \
                                (pot, dh.i + 1, dh.j + 1, dh.k + 1, dh.l + 1,
                                 dh.par[0], dh.par[1], dh.par[2], 0.5, 0.5, dh.name))
                for di in mol.dimpr:
                    if cos4:
                        pot = 'cos4'
                        f.write('%4s %4d %4d %4d %4d %9.4f %9.4f %9.4f %9.4f'\
                                ' %6.3f %6.3f  # %s\n' % \
                                (pot, di.i + 1, di.j + 1, di.k + 1, di.l + 1,
                                 di.par[0], di.par[1], di.par[2], di.par[3], 0.5, 0.5,
                                 di.name))
                    else:
                        pot = 'cos3'
                        f.write('%4s %4d %4d %4d %4d %9.4f %9.4f %9.4f'\
                                ' %6.3f %6.3f  # %s\n' % \
                                (pot, di.i + 1, di.j + 1, di.k + 1, di.l + 1,
                                 di.par[0], di.par[1], di.par[2], 0.5, 0.5, di.name))
                f.write('finish\n')

            f.write('vdw %d\n' % len(self.vdw))
            for nb in self.vdw:
                if nb.pot == 'lj':
                    f.write('%-5s %-5s %4s %10.6f %8.4f\n' % \
                            (nb.i, nb.j, nb.pot, nb.par[1], nb.par[0]))
                
            f.write('close\n')

        name = []
        x = []
        y = []
        z = []
        try:
            with open('simbox.xyz', 'r') as fx:
                natoms = int(fx.readline().strip())
                title = fx.readline().strip()
                i = 0
                while i < natoms:
                    tok = fx.readline().strip().split()
                    name.append(tok[0])
                    x.append(float(tok[1]))
                    y.append(float(tok[2]))
                    z.append(float(tok[3]))
                    i += 1
        
            with open('CONFIG', 'w') as fc:
                fc.write(title + '\n')
                fc.write(' %9d %9d %9d\n' % (0, 1, natoms))
                tol = 4.0
                fc.write(' %19.9f %19.9f %19.9f\n' % (boxlen + tol, 0.0, 0.0))
                fc.write(' %19.9f %19.9f %19.9f\n' % (0.0, boxlen + tol, 0.0))
                fc.write(' %19.9f %19.9f %19.9f\n' % (0.0, 0.0, boxlen + tol))
                i = 0
                while i < natoms:
                    fc.write('%-8s %9d\n' % (name[i], i + 1))
                    fc.write(' %19.9f %19.9f %19.9f\n' % (x[i], y[i], z[i]))
                    i += 1

        except IOError:
            print 'warning: cannot open simbox.xyz, CONFIG not created'

            
# --------------------------------------


def main():
    parser = argparse.ArgumentParser(
        description = 'Force-field parameters and atomic coordinates for molecules'\
        ' described by z-matrices. '\
        'Simulation box can be built using packmol with the file produced')
    parser.add_argument('-r', '--rho', type=float, default = 5.0,
                        help = 'density in mol/L for packmol (default: 5.0)')
    parser.add_argument('-x', '--mix', default = 'g',
                        help = '[a]rithmetic or [g]eometric sigma_ij (default: g)')
    parser.add_argument('-l', '--lammps', action = 'store_true', 
                        help = 'save in lammps format '\
                        '(needs simbox.xyz built using packmol)')
    parser.add_argument('-a', '--allpairs', action = 'store_true', 
                        help = 'write all I J pairs to lammps input files')
    parser.add_argument('-d', '--dlpoly', action = 'store_true',
                        help = 'save in dlpoly format '\
                        '(needs simbox.xyz built using packmol)')
    parser.add_argument('-c', '--cos4', action = 'store_true', 
                        help = 'use cos4 dihedrals in dlpoly FIELD')
    parser.add_argument('-s', '--symbol', action = 'store_true',
                        help = 'write atomic symbols to xyz files '\
                        '(useful for rasmol)')
    parser.add_argument('-q', '--quiet', action = 'store_true')
    parser.add_argument('infile', nargs='+',
                        help = 'file(s) with z-matrix: n1 infile1 [n2 infile2 '\
                        '...], where ni are the numbers of molecules in '\
                        'the system')
    args = parser.parse_args()

    if len(args.infile) == 1:
        nmols = [1]
        files = args.infile
    else:
        nmols = args.infile[::2]   # even elements are numbers of molecules
        files = args.infile[1::2]  # odd elements are zmat files

    m = []
    i = 0
    nmol = 0
    if not args.quiet:
        print 'atomic coordinates'
    for zfile in files:
        m.append(mol(zfile))
        m[i].nmols = int(nmols[i])
        nmol += m[i].nmols
        if not args.quiet:
            print '  ' + zfile.rsplit('.', 1)[0] + '.xyz'
        m[i].writexyz(args.symbol)
        i += 1

    boxlen = math.pow(nmol / (args.rho * 6.022e+23 * 1.0e-27), 1./3.) 

    s = system(m, args.mix)
    if not args.quiet:
        print 'charges'
        for spec in m:
            print '  %+.3f' % spec.charge() 
        print 'packmol input\n  pack.inp'
        s.writepackmol(boxlen)

    if args.lammps:
        if not args.quiet:
            print 'force field and coordinates\n  in.lmp\n  data.lmp'
        s.writelmp(boxlen, args.mix, args.allpairs)
    elif args.dlpoly:
        if not args.quiet:
            print 'force field and coordinates\n  FIELD\n  CONFIG'
        s.writedlp(boxlen, args.cos4)


if __name__ == '__main__':
    main()

