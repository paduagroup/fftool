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


How to use
----------

How to build an initial configuration of systems composed of
molecules, ions or materials.

1. For each molecule, ion or fragment of a material prepare a file
   with atomic coordinates and/or connectivity (covalent bonds). The
   formats accepted by this tool are `.zmat`, `.mol` or `.xyz`.

    A `.zmat` file has the molecule name in the first line, thane one
    empy line, then the z-matrix. See the `examples` directory and the
    Wikipedia entry for "Z-matrix (chemistry)". Variables can be used
    for distances, angles and dihedrals. Connectivity is inferred from
    the z-matrix. Cyclic molecules require additional `connect`
    records to close rings. Improper dihedrals must be indicated by
    additional `improper` records. After the z-matrix the name of a
    file with force field parameters can be supplied.

    A MDL `.mol` file contains a table with coordinates and also
    bonds. The name of a file with force field parameters can be
    given in the first line after the molecule name, or in the third
    line.

    A `.xyz` file contains atomic coordinates only. The name of a file
    with force field parameters can be given in the second line after
    the molecule name, and in this case connectivity is inferred from
    the bond lengths.

    There are many free tools to create MDL mol files, xyz files or
    z-matrices, which are common formats in computational chemistry
    ([Open Babel](http://openbabel.org/),
    [Avogadro](http://avogadro.cc/)). Manual editing of the files is
    usually necessary in order to match the atom names with those of
    the force field. (One disadvantage of using `.xyz` or `.mol` files
    is the difficulty to include improper dihedrals, which may be
    necessary for certain molecules.)

2. Use the `fftool.py` script to create `.xyz` files with atomic
   coordinates for the components of your system, plus an input file for
   `packmol`. For help type `fftool.py -h`. To build a simulation box
   with 40 ethanol and 300 water molecules and a density of 40.0 mol/L
   do:

        fftool.py 40 ethanol.zmat 300 spce.zmat --rho 40.0

3. Use `packmol` with the `pack.inp` file just created to buid the
   simulation box (adjust the density if necessary):

        packmol < pack.inp

    For more complex spatial arrangements of molecules and materials
    you can modify the `pack.inp` file to suit your needs (see the
    [Packmol](http://www.ime.unicamp.br/~martinez/packmol/)
    documentation).  Atomic coordinates for the full system will be
    written to `simbox.xyz`.

4. Use `fftool.py` to build the input files for LAMMPS or DL_POLY
   containing the force field parameters and the coordinates:

        fftool.py 40 ethanol.zmat 300 spce.zmat --rho 40.0 --lammps

    If no force field information is given explicitly in the molecule
    files, a default LJ potential with parameters zeroed will be
    assigned. No terms for bonds, angles or torsions will be
    created. This is suitable when working with non-additive,
    bond-order or other potentials often used for materials. The input
    files for MD simulations will have to be edited manually to
    include an interaction potential for the material.


References
----------

* [Packmol](http://www.ime.unicamp.br/~martinez/packmol/):
  L. Martinez et al. J Comp Chem 30 (2009) 2157, DOI:
  [10.1002/jcc.21224](http://dx.doi.org/10.1002/jcc.21224) 
  
* [LAMMPS](http://lammps.sandia.gov/): S. Plimton, J Comp Phys
  117 (1995) 1, DOI:
  [10.1006/jcph.1995.1039](http://dx.doi.org/10.1006/jcph.1995.1039)

* [DL_POLY](http://www.stfc.ac.uk/CSE/randd/ccg/software/DL_POLY/25526.aspx): I.T. Todorov and W. Smith, Daresbury Lab. 

