fftool
======

_[Agilio Padua](http://tim.univ-bpclermont.fr/apadua)_

Tool to build force field input files for molecular dynamics. 

Contents
--------

* `fftool.py`: python script to build simulation box of molecular or
    ionic liquids and their mixtures. Requires the
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

    git clone https://github.com/agiliopadua/ffool.git


How to use
----------

How to build an initial configuration of a molecular or ionic system.

1. For each molecule or ion prepare a file containing a z-matrix
   (`molecule.zmat`). See the `examples` directory and check the
   Wikipedia entry for "Z-matrix (chemistry)". The `fftool.py` script
   determines the connectivity (which atoms are linked by covalent
   bonds) from the z-matrix. Cyclic molecules require additional
   `connect` records to close rings. Improper dihedrals cannot be
   inferred from connectivity and must be indicated by `improper`
   records. After the z-matrix supply the name of a database of force
   field parameters (`database.ff`).

2. Use the `fftool.py` script to create `.xyz` files for the molecules
   in your system and an input file for `packmol`. For help type
   `fftool.py -h`. To build a simulation box with 40 ethanol and 300
   water molecules and a
   density of 40.0 mol/L do:

        fftool.py 40 ethanol.zmat 300 spce.zmat --rho 40.0

3. Use `packmol` with the `pack.inp` file just created to buid the
   simulation box (adjust the density if necessary):

        packmol < pack.inp

    Atom coordinates will be written to `simbox.xyz`. You can use a
    molecular viewer such as RasMol or VMD to look at the `.xyz` files
    (`fftool.py` has an option to write IUPAC atomic symbols instead
    of the atom names from the force field).

4. Use `fftool.py` to build the input files for LAMMPS or DL_POLY
   containing the force field parameters and the coordinates:

        fftool.py 40 ethanol.zmat 300 spce.zmat --rho 40.0 --lammps


References
----------

* [Packmol](http://www.ime.unicamp.br/~martinez/packmol/):
  L. Martinez et al. J Comp Chem 30 (2009) 2157, DOI:
  [10.1002/jcc.21224](http://dx.doi.org/10.1002/jcc.21224) 
  
* [LAMMPS](http://lammps.sandia.gov/): S. Plimton, J Comp Phys
  117 (1995) 1, DOI:
  [10.1006/jcph.1995.1039](http://dx.doi.org/10.1006/jcph.1995.1039)

* [DL_POLY](http://www.stfc.ac.uk/CSE/randd/ccg/software/DL_POLY/25526.aspx): I.T. Todorov and W. Smith, Daresbury Lab. 

