#!/usr/bin/env python

import argparse
import json

class forcefield(object):
    '''force field parameter database'''

    def __init__(self, filename):

        self.ff = {}

        for line in open(filename, 'r'):
            if line.startswith('#') or line.strip() == '':
                continue
            if line.lower().startswith('atom'):
                section = 'atoms'
                self.ff['Atoms'] = []
                continue
            elif line.lower().startswith('bond'):
                section = 'bonds'
                self.ff['Bonds'] = []
                continue
            elif line.lower().startswith('angl'):
                section = 'angles'
                self.ff['Angles'] = []
                continue
            elif line.lower().startswith('dihe'):
                section = 'dihedrals'
                self.ff['Dihedrals'] = []
                continue
            elif line.lower().startswith('impro'):
                section = 'improper'
                self.ff['Improper'] = []
                continue
            tok = line.strip().split()
            if section == 'atoms':
                at = {}
                at['type'] = tok[0]
                at['class'] = tok[1]
                at['mass'] = float(tok[2])
                at['charge'] = float(tok[3])
                at['potential'] = tok[4]
                at['sigma'] = float(tok[5])
                at['epsilon'] = float(tok[6])
                self.ff['Atoms'].append(at)
            elif section == 'bonds':
                bd = {}
                bd['class1'] = tok[0]
                bd['class2'] = tok[1]
                bd['potential'] = tok[2]
                bd['length'] = float(tok[3])
                bd['k'] = float(tok[4])
                self.ff['Bonds'].append(bd)
            elif section == 'angles':
                an = {}
                an['class1'] = tok[0]
                an['class2'] = tok[1]
                an['class3'] = tok[2]
                an['potential'] = tok[3]
                an['angle'] = float(tok[4])
                an['k'] = float(tok[5])
                self.ff['Angles'].append(an)
            elif section == 'dihedrals':
                di = {}
                di['class1'] = tok[0]
                di['class2'] = tok[1]
                di['class3'] = tok[2]
                di['class4'] = tok[3]
                di['potential'] = tok[4]
                di['v1'] = float(tok[5])
                di['v2'] = float(tok[6])
                di['v3'] = float(tok[7])
                di['v4'] = float(tok[8])
                self.ff['Dihedrals'].append(di)
            elif section == 'improper':
                im = {}
                im['class1'] = tok[0]
                im['class2'] = tok[1]
                im['class3'] = tok[2]
                im['class4'] = tok[3]
                im['potential'] = tok[4]
                im['v1'] = float(tok[5])
                im['v2'] = float(tok[6])
                im['v3'] = float(tok[7])
                im['v4'] = float(tok[8])
                self.ff['Improper'].append(im)

    def write(self, outfile):
        '''write force field dict to json file'''

        with open(outfile, 'w') as file: 
            json.dump(self.ff, file, indent=4)


def main():
    parser = argparse.ArgumentParser(description =
        'convert force field file to json',
        formatter_class = argparse.RawTextHelpFormatter)
    parser.add_argument('infile', help = 'force field file')
    parser.add_argument('outfile', help = 'output OpenMM xml file')
    args = parser.parse_args()

    ff = forcefield(args.infile)
    ff.write(args.outfile)


if __name__ == '__main__':
    main()
