fftool
======

_[Agilio Padua](http://tim.univ-bpclermont.fr/apadua)_

Tool to build force field input files for molecular dynamics. 

Contents
--------

* `fftool.py`: python script to build simulation box of molecular or
    ionic liquids and their mixtures, as well as materials. It
    requires the
    [Packmol](http://www.ime.unicamp.br/~martinez/packmol/) software
    to create coordinates. Force field files are written in formats
    suitable for the [LAMMPS](http://lammps.sandia.gov/) or
    [DL_POLY](http://www.stfc.ac.uk/CSE/randd/ccg/software/DL_POLY/25526.aspx)
    molecular dynamics packages.

Requirements
------------

* [Python 2.7](http://www.python.org/)

* [Packmol](http://www.ime.unicamp.br/~martinez/packmol/)


Obtaining
---------

Download the files or else clone the repository (easier to stay updated):

    git clone https://github.com/agiliopadua/fftool.git


How to Use
----------

How to build an initial configuration of systems composed of
molecules, ions or materials.

1. For each molecule, ion or fragment of a material prepare a file
   with atomic coordinates and/or connectivity (covalent bonds). The
   formats accepted by this tool are `.zmat`, `.mol` or `.xyz`.

    A `.zmat` file has the molecule name in the first line, followed
    by one empy line, then the z-matrix. See the `examples` directory
    and the Wikipedia entry for "Z-matrix (chemistry)". Variables can
    be used for distances, angles and dihedrals. Connectivity is
    inferred from the z-matrix by default. In this case cyclic
    molecules require additional `connect` records to close rings. If
    a `reconnect` record is present, then connectivity is guessed
    based on bond distances from the force field. After the z-matrix
    the name of a file with force field parameters can be supplied.

    MDL `.mol` is a standard file format, containing a table with
    coordinates and also bonds. The name of a file with force field
    parameters can be given in the first line after the molecule name,
    or in the third line. If the keyword `reconnect` is present after
    the force field filename, then connectivity is guessed based on
    bond distances from the force field.

    The `.xyz` format is also standard, containing atomic coordinates
    only. The name of a file with force field parameters can be given
    in the second line after the molecule name, and in this case
    connectivity is inferred from the bond lengths.

    There are many free tools to create MDL mol files, xyz files or
    z-matrices, which are common formats in computational chemistry
    ([Open Babel](http://openbabel.org/),
    [Avogadro](http://avogadro.cc/)). Manual editing of the files is
    usually necessary in order to match the atom names with those of
    the force field.

2. Use the `fftool.py` script to create `.xyz` files with atomic
   coordinates for the components of your system, plus an input file for
   `packmol`. For help type `fftool.py -h`. To build a simulation box
   with 40 ethanol and 300 water molecules and a density of 38.0 mol/L
   do:

        fftool.py 40 ethanol.zmat 300 spce.zmat -r 38.0

    Alternatively the side length of the the simulation box (cubic) in
    angstroms can be supplied:

        fftool.py 40 ethanol.zmat 300 spce.zmat -b 20.0

3. Use `packmol` with the `pack.inp` file just created to build the
   simulation box (adjust the density/box size if necessary):

        packmol < pack.inp

    For more complex spatial arrangements of molecules and materials
    you can modify the `pack.inp` file to suit your needs (see the
    [Packmol](http://www.ime.unicamp.br/~martinez/packmol/)
    documentation).  Atomic coordinates for the full system will be
    written to `simbox.xyz`.

4. Use `fftool.py` to build the input files for LAMMPS or DL_POLY
   containing the force field parameters and the coordinates of all
   the atoms from `simbox.xyz`:

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
atoms (sp2). A proper dihedral i-j-k-l is defined between bonded atoms
i-j, j-k, and k-l and corresponds to torsion around bond j-k, the
dihedral being the angle between planes ijk and jkl. An improper
dihedral i-j-k-l is defined between bonded atoms i-k, j-k and k-l,
therefore k is a central atom bonded to the three others. Often the
same functional form used for proper torsions is also used for
improper dihedrals.

The `fftool.py` script searches for improper dihedrals on all atoms
with three bonds (which may be planar) so a number of warning messages
may be printed and can be ignored if the atoms in question are not
centers of improper torsions.

The number and order of the atoms in the true improper dihedrals
should be verified in the files created.


Periodic Boundary Conditions
----------------------------

For molecular systems the initial configuration will not contain
molecules crossing the boundaries of the simulation box. If the size
of the box is indicated by just one value, `-b L`, or by supplying the
density, then the box will be cubic and an extra space of 1 A is added
in each dimension to avoid overlaps in the initial configuration (as
explained in the `packmol` documentation).

For simulations with extended materials it is possible to create
chemical bonds across boundaries. The option `-p` allows specification
of periodic conditions in x, y, z or combinations thereof. It is
important in this case to supply precise dimensions for the simulation
box using the option `-b Lx,Ly,Lz`, even if the box is cubic. In this
manner no additional space will be added. The coordinates of the atoms
of the material have to be prepared carefully, so that distances
across periodic boundaries are within the tolerance to identify bonds.

Also, the `pack.inp` filed will likely need editing in order to
position the atoms of the material precisely.


References
----------

* [Packmol](http://www.ime.unicamp.br/~martinez/packmol/):
  L. Martinez et al. J Comp Chem 30 (2009) 2157, DOI:
  [10.1002/jcc.21224](http://dx.doi.org/10.1002/jcc.21224) 
  
* [LAMMPS](http://lammps.sandia.gov/): S. Plimton, J Comp Phys
  117 (1995) 1, DOI:
  [10.1006/jcph.1995.1039](http://dx.doi.org/10.1006/jcph.1995.1039)

* [DL_POLY](http://www.stfc.ac.uk/CSE/randd/ccg/software/DL_POLY/25526.aspx): I.T. Todorov and W. Smith, Daresbury Lab. 

