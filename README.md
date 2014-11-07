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

How to build an initial configuration of a molecular or ionic system.

1. For each molecule, ion or fragment of a material prepare a file
   with atomic coordinates and/or connectivity (covalent bonds). The
   formats accepted by this tool are `.zmat`, `.mol` or `.xyz`.

    A `.zmat` file has the molecule name in the first line, then one
    empty ine, then the z-matrix. See the `examples` directory and
    check the Wikipedia entry for "Z-matrix (chemistry)". Variables
    can be used for distances, angles and dihedrals. The connectivity
    is inferred from the z-matrix. Cyclic molecules require additional
    `connect` records to close rings. Improper dihedrals must be
    indicated by additional `improper` records. After the z-matrix the
    name of a database of force field parameters (`database.ff`) can
    be supplied.

    A MDL `.mol` file contains a table with coordinates and also
    bonds. The name of a database of force field parameters can be
    given in the first line after the molecule name, or in the third
    line.

    A `.xyz` file contains atomic coordinates only. The name of a
    database of force field parameters can be given in the second line
    after the molecule name. If a force field is given, then the
    connectivity is inferred from the equilibrium bond lengths in the
    force field.

    There are many freely available tools to create MDL mol files or
    xyz files, which are common formats in computational
    chemistry. Z-matrices can also be created by such tools. The
    `.zmat` format has the advantage of allowing specification of
    improper dihedrals. Once these are created, manual editing is
    necessary to match the names of atoms with those of the force
    field.

    If no force field database is provided, this tool will assign a
    default Lennard-Jones potential to each atom type, with parameters
    zeroed. This is useful to work with non-additive, bond-order or
    other models, often used for materials. The input files for the MD
    simulations will have to be edited afterwards to include such
    interaction models.

2. Use the `fftool.py` script to create `.xyz` files with atomic
   coordinates for the molecules in your system plus an input file for
   `packmol`. For help type `fftool.py -h`. To build a simulation box
   with 40 ethanol and 300 water molecules and a density of 40.0 mol/L
   do:

        fftool.py 40 ethanol.zmat 300 spce.zmat --rho 40.0

3. Use `packmol` with the `pack.inp` file just created to buid the
   simulation box (adjust the density if necessary):

        packmol < pack.inp

    For more complex spacial arrangements of molecules and materials,
    you can modify the `pack.inp` file to suit your needs (see the
    [Packmol](http://www.ime.unicamp.br/~martinez/packmol/)
    documentation).  Atom coordinates will be written to
    `simbox.xyz`. You can use a molecular viewer such as RasMol or VMD
    to look at the `.xyz` files (`fftool.py` has an option to write
    IUPAC atomic symbols instead of the atom names from the force
    field).

4. Use `fftool.py` to build the input files for LAMMPS or DL_POLY
   containing the force field parameters and the coordinates:

        fftool.py 40 ethanol.zmat 300 spce.zmat --rho 40.0 --lammps

    If no force field information is given explicitly in the molecule
    files, those atoms will ba assigned a default LJ potential with
    parameters zeroed. The input files for MD simulations will have to
    be edited manually to include the correct potential function and
    parameters for the material.


References
----------

* [Packmol](http://www.ime.unicamp.br/~martinez/packmol/):
  L. Martinez et al. J Comp Chem 30 (2009) 2157, DOI:
  [10.1002/jcc.21224](http://dx.doi.org/10.1002/jcc.21224) 
  
* [LAMMPS](http://lammps.sandia.gov/): S. Plimton, J Comp Phys
  117 (1995) 1, DOI:
  [10.1006/jcph.1995.1039](http://dx.doi.org/10.1006/jcph.1995.1039)

* [DL_POLY](http://www.stfc.ac.uk/CSE/randd/ccg/software/DL_POLY/25526.aspx): I.T. Todorov and W. Smith, Daresbury Lab. 

