#!/usr/bin/env python

import argparse
import xml.etree.ElementTree as ET


class forcefield(object):
    '''force field parameter database'''

    def __init__(self, filename):

        self.ftree = ET.ElementTree(ET.Element('ForceField'))
        root = self.ftree.getroot()

        for line in open(filename, 'r'):
            if line.startswith('#') or line.strip() == '':
                continue
        
            if line.lower().startswith('atom'):
                section = 'atoms'
                atoms = ET.SubElement(root, 'Atoms')
                continue
            elif line.lower().startswith('bond'):
                section = 'bonds'
                bonds = ET.SubElement(root, 'Bonds')
                continue
            elif line.lower().startswith('angl'):
                section = 'angles'
                angles = ET.SubElement(root, 'Angles')
                continue
            elif line.lower().startswith('dihe'):
                section = 'dihedrals'
                if not root.find('Torsions'):
                    torsions = ET.SubElement(root, 'Torsions')
                continue
            elif line.lower().startswith('impro'):
                section = 'improper'
                if not root.find('Torsions'):
                    torsions = ET.SubElement(root, 'Torsions')
                continue
            tok = line.strip().split()
            if section == 'atoms':
                at = ET.SubElement(atoms, 'Atom')
                at.set('type', tok[0])
                at.set('class', tok[1])
                at.set('mass', tok[2])
                at.set('charge', tok[3])
                at.set('potential', tok[4])
                at.set('sigma', tok[5])
                at.set('epsilon', tok[6])
            elif section == 'bonds':
                bd = ET.SubElement(bonds, 'Bond')
                bd.set('class1', tok[0])
                bd.set('class2', tok[1])
                bd.set('potential', tok[2])
                bd.set('length', tok[3])
                bd.set('k', tok[4])
            elif section == 'angles':
                an = ET.SubElement(angles, 'Angle')
                an.set('class1', tok[0])
                an.set('class2', tok[1])
                an.set('class3', tok[2])
                an.set('potential', tok[3])
                an.set('angle', tok[4])
                an.set('k', tok[5])
            elif section == 'dihedrals':
                di = ET.SubElement(torsions, 'Proper')
                di.set('class1', tok[0])
                di.set('class2', tok[1])
                di.set('class3', tok[2])
                di.set('class4', tok[3])
                di.set('potential', tok[4])
                di.set('v1', tok[5])
                di.set('v2', tok[6])
                di.set('v3', tok[7])
                di.set('v4', tok[8])
            elif section == 'improper':
                im = ET.SubElement(torsions, 'Improper')
                im.set('class1', tok[0])
                im.set('class2', tok[1])
                im.set('class3', tok[2])
                im.set('class4', tok[3])
                im.set('potential', tok[4])
                im.set('v1', tok[5])
                im.set('v2', tok[6])
                im.set('v3', tok[7])
                im.set('v4', tok[8])

    def write(self, outfile):
        '''write force field to xml file'''

        ET.indent(self.ftree.getroot(), space=' ')
        self.ftree.write(outfile)


def main():
    parser = argparse.ArgumentParser(description =
        'convert force field file to xml',
        formatter_class = argparse.RawTextHelpFormatter)
    parser.add_argument('infile', help = 'force field file')
    parser.add_argument('outfile', help = 'output OpenMM xml file')
    args = parser.parse_args()

    ff = forcefield(args.infile)
    ff.write(args.outfile)


if __name__ == '__main__':
    main()
