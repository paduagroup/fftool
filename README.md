fftool
======

_[Agilio Padua](http://tim.univ-bpclermont.fr/apadua)_

This is a tool to build force field input files for molecular
dynamics.


Contents
--------

* `fftool.py`: builds a simulation box and the corresponding force
    field for systems containing molecules, ions or extended
    materials. It requires the
    [Packmol](http://www.ime.unicamp.br/~martinez/packmol/) software
    to create coordinates. The output are files in formats suitable
    for the [LAMMPS](http://lammps.sandia.gov/) or
    [DL_POLY](http://www.stfc.ac.uk/CSE/randd/ccg/software/DL_POLY/25526.aspx)
    molecular dynamics packages.

* `tools/`: several utility scripts.

* `examples/`: examples of molecule files and force field databases.


Requirements
------------

* [Python](http://www.python.org/) versions 2.7 or 3 should work.

* [Packmol](http://www.ime.unicamp.br/~martinez/packmol/) to pack
  molecules and materials in the simultion box.

* [PyPy](http://pypy.org) (optional) can bring enormous speed
  improvements for large systems.


Obtaining
---------

Download the files or else clone the repository (easier to stay updated):

    git clone https://github.com/agiliopadua/fftool.git


Tutorial
--------

These are instructions on how to build an initial configuration for a
system composed of molecules, ions or materials.

1. For each molecule, ion or fragment of a material prepare a file
   with atomic coordinates and/or connectivity (covalent bonds). The
   formats accepted by this tool are `.zmat`, `.mol` or `.xyz`.

    A `.zmat` file has the molecule name in the first line, followed
    by one empy line, then the z-matrix. See the `examples` directory
    and the Wikipedia entry for "Z-matrix (chemistry)". Variables can
    be used in place of distances, angles and dihedrals. Connectivity
    is inferred from the z-matrix by default. In this case cyclic
    molecules require additional `connect` records to close
    rings. Improper dihedrals can be indicated by `improper`
    records. If a `reconnect` record is present, then connectivity
    will be guessed based on bond distances from the force
    field. After the z-matrix the informations above, the name of a
    file with force field parameters can be supplied.

    MDL `.mol` format is widely used, containing a table with
    coordinates and also bonds. The name of a file with force field
    parameters can be given in the first line after the molecule name,
    or in the third line. If the keyword `reconnect` is present after
    the force field filename, then connectivity will be guessed based
    on bond distances from the force field.

    The `.xyz` is also widely used, containing atomic coordinates
    only. The name of a file with force field parameters can be given
    in the second line after the molecule name, and in this case
    connectivity is inferred from the bond lengths in the force field.

    There are many free tools to create MDL mol files, xyz files or
    z-matrices, which are common formats in computational chemistry
    ([Open Babel](http://openbabel.org/),
    [Avogadro](http://avogadro.cc/)). Manual editing of the files is
    usually necessary in order to match the atom names with those of
    the force field.

2. Use the `fftool.py` script to create `.xyz` files with atomic
   coordinates for the components of your system, plus an input file
   for `packmol`. For help type `fftool.py -h`. For example, to build
   a simulation box with 40 ethanol and 300 water molecules and a
   density of 38.0 mol/L do:

        fftool.py 40 ethanol.zmat 300 spce.zmat -r 38.0

    Alternatively the side length of the the simulation box (cubic) in
    angstroms can be supplied:

        fftool.py 40 ethanol.zmat 300 spce.zmat -b 20.0

3. Use `packmol` with the `pack.inp` file just created to generate the
   atomic coordinates in the simulation box (adjust the density/box
   size if necessary):

        packmol < pack.inp

    For more complex spatial arrangements of molecules and materials
    you can modify the `pack.inp` file to suit your needs (see the
    [Packmol](http://www.ime.unicamp.br/~martinez/packmol/)
    documentation).  Atomic coordinates for the full system will be
    written to `simbox.xyz`.

4. Use `fftool.py` to build the input files for LAMMPS or DL_POLY
   containing the force field parameters and the coordinates of all
   the atoms (from `simbox.xyz`):

        fftool.py 40 ethanol.zmat 300 spce.zmat --r 38.0 -l

    If no force field information was given explicitly in the molecule
    files, a default LJ potential with parameters zeroed will be
    assigned to atoms. No terms for bonds, angles or torsions will be
    created. This is suitable when working with non-additive,
    bond-order or other potentials often used for materials. The input
    files for MD simulations will have to be edited manually to
    include an interaction potential for the material.


Improper Dihedrals
------------------

Improper dihedrals are often used to increase the rigidity of planar
atoms (sp2), and differ from proper dihedrals in how they are
defined. A proper dihedral i-j-k-l is defined between bonded atoms
i-j, j-k, and k-l and corresponds to torsion around bond j-k, the
dihedral being the angle between planes ijk and jkl. An improper
dihedral i-j-k-l is defined between bonded atoms i-k, j-k and k-l,
therefore k is a central atom bonded to the three others. This central
atom of the improper dihedral is assumed to be the third in the
list. Often in the force field the same functional form is used both
for proper and improper torsions.

If `improper` records are supplied in a molecule file (in `.zmat`
format) then those improper dihedrals are assumed by
`fftool.py`. Otherwise, the script will search for candidate improper
dihedrals on all atoms with three bonds, with any of `.zmat`, `.mol`
or `.xyz` input formats. A number of warning messages may be printed
if there are atoms with three bonds, which can be ignored if the atoms
in question are not centers of improper torsions. The number and order
of the atoms in the true improper dihedrals should be verified in the
files created.


Periodic Boundary Conditions
----------------------------

For molecular systems the initial configuration will generaly not
contain molecules crossing the boundaries of the simulation box. If
the size of the box is indicated by just one value, `-b L`, or by
supplying the density, then the box will be cubic, centered on the
origin, and an extra space of 1 A is reserved in each dimension to
avoid overlaps in the initial configuration (as explained in the
`packmol` documentation).

For simulations with extended materials it is possible to create
chemical bonds across boundaries. The option `-p` allows specification
of periodic conditions along x, y, z or combinations thereof. It is
important in this case to supply precise dimensions for the simulation
box using the option `-b Lx,Ly,Lz`. An energy minimization step prior
to start of the MD simulation is highly recommended.

The coordinates of the atoms of the material have to be supplied in
`.xyz` format and prepared carefully so that distances across periodic
boundaries are within the tolerance to identify bonds. The number of
bonds in the files created should be verified.

It is important that only the material for which bonds are to be
established across boudaries is supplied in `.xyz` format. The initial
files for other molecules in the system should be in `.zmat` or `.mol`
formats, which containing connectivity information. This is to avoid
spurious bonds between atoms that happen to be positioned too close to
boudaries.

The `pack.inp` file will need manual editing in order to position the
atoms of the material precisely.



Force Field File Format
-----------------------

The `fftool.py` script reads a database of molecular force field terms
in the format described below. See the `examples` directory.

Blank lines and lines starting with `#` are ignored.

There are five sections, with headings `ATOMS`, `BONDS`, `ANGLES`,
`DIHEDRALS` and `IMPROPER`. Under each section heading, registers
concerning the different types of term in the force field are given.

`ATOMS` records describe, for each type of atom:
- the non-bonded atom type used for intermolecular interactions (these
  types may differ in the charges or intermolecular potential
  parameters)
- the bonded atom type used in intermolecular interactions (these
  types determine the intramolecular terms such as bonds, angles
  dihedrals)
- the mass in atomic units
- the electrostatic charge in elementary units
- the non-bonded potential type, e.g. _lj_
- potential parameters, such as Lennard-Jones _sigma_ and _epsilon_

        C3H   CT  12.011  -0.18   lj    3.50   0.27614

`BONDS` records describe covalent bonds between intramolecular atom types:
- two bonded atom types
- type of bond potential, e.g. _harm_ for harmonic potential or
  _cons_ for a constrained bond.
- bond potential parameters, with the force constant in the form k/2
  (x - x0)^2

        CT  CT   harm   1.529   2242.6

`ANGLES` records describe valence angles between intramolecular atom types:
- three bonded atom types, in which the central atom is bonded to the other
  two, e.g. i-j and j-k are bonded.
- type of angle potential, e.g. _harm`_ for harmonic potential or
  _cons_ for a constrained angle.
- angle potential parameters, with the force constant in the form k/2
  (x - x0)^2

        HC  CT  CT   harm   110.7   313.8

`DIHEDRALS` records describe torsion angles between intramolecular
atom types:
- four bonded atom types, in which atoms i-j, j-k, k-l are bonded.
- type of dihedral potential, e.g. _opls_ for OPLS cosine series with
  four terms.
- dihedral potential parameters, with the coefficients in the form V_n/2
  (1 +/- cos(n phi)).

        CT  CT  CT  CT   opls    5.4392   -0.2092    0.8368    0.0000

`IMPROPER` records describe improper dihedral angles between
intramolecular atom types:
- four bonded atom types, in which atoms i-k, j-k, k-l are bonded.
- type of dihedral potential, e.g. _opls_ for OPLS cosine series with
  four terms.
- dihedral potential parameters, with the coefficients in the form V_n/2
(1 +/- cos(n phi)).


References
----------

* [Packmol](http://www.ime.unicamp.br/~martinez/packmol/):
  L. Martinez et al. J Comp Chem 30 (2009) 2157, DOI:
  [10.1002/jcc.21224](http://dx.doi.org/10.1002/jcc.21224) 
  
* [LAMMPS](http://lammps.sandia.gov/): S. Plimton, J Comp Phys
  117 (1995) 1, DOI:
  [10.1006/jcph.1995.1039](http://dx.doi.org/10.1006/jcph.1995.1039)

* [DL_POLY](http://www.stfc.ac.uk/CSE/randd/ccg/software/DL_POLY/25526.aspx):
  I.T. Todorov and W. Smith, Daresbury Lab. 

