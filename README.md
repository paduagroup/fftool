fftool
======

[![DOI](https://zenodo.org/badge/doi/10.5281/zenodo.18618.svg)](http://dx.doi.org/10.5281/zenodo.18618)

_[Agilio Padua](http://perso.ens-lyon.fr/agilio.padua)_

This is a Python tool to build force field input files for molecular dynamics.


Contents
--------

* `fftool`: builds a simulation box and the corresponding force
    field (non-polarizable or drude based polarizable) for systems containing molecules, ions or extended
    materials. It requires the
    [Packmol](http://www.ime.unicamp.br/~martinez/packmol/) software
    to generate coordinates in the box. The output are files in formats
    suitable for the [LAMMPS](http://lammps.sandia.gov/),
    [DL_POLY](http://www.stfc.ac.uk/CSE/randd/ccg/software/DL_POLY/25526.aspx), 
    [GROMACS](http://www.gromacs.org) and [NAMD](https://www.ks.uiuc.edu/Research/namd/)
    molecular dynamics packages.

* `tools/`: utility scripts.

* `examples/`: examples of molecule files and force field databases.


Requirements
------------

* [Python](http://www.python.org/)

* [Packmol](http://www.ime.unicamp.br/~martinez/packmol/) to pack
  molecules and materials in the simultion box.

* [PyPy](http://pypy.org) (optional) can bring speed
  improvements for large systems.


Obtaining
---------

Download the files or clone the repository:

git clone https://github.com/agiliopadua/fftool.git


Tutorial
--------

These are instructions on how to build an initial configuration for a
system composed of molecules, ions or materials.

1. For each molecule, ion or fragment of a material prepare a file with atomic
   coordinates and eventually connectivity (covalent bonds). The formats
   accepted by this tool are `.zmat`, `.mol`, `.pdb` or `.xyz`, which are
   common formats in computational chemistry.

   A `.zmat` file has the molecule name in the first line, followed by one
   empy line, then the z-matrix. See the `examples` directory and the
   Wikipedia entry for "Z-matrix (chemistry)". Variables can be used in place
   of distances, angles and dihedrals. Connectivity is inferred from the
   z-matrix by default. In this case cyclic molecules require additional
   `connect` records to close rings. Improper dihedrals can be indicated by
   `improper` records. If a `reconnect` record is present, then connectivity
   will be guessed based on bond distances from the force field (see
   below). After the z-matrix and the informations above, the name of a file
   with force field parameters can be supplied.

   The MDL Molfile `.mol` file format contains a table with coordinates and
   also bonds. The name of a file with force field parameters can be given in
   the first line after the molecule name or in the third line. If the keyword
   `reconnect` is present after the force field filename, then connectivity
   will be deduced based on bond distances from the force field.

   The PDB file format `.pdb` is widely used for proteins. The name of a file
   with force field parameters can be given on a `COMPND` record after the
   molecule name.  Connectivity is deduced from the bond lengths in the force
   field (`CONNECT` records are not read).

   The XYZ file format `.xyz` contains atomic coordinates only. The name of a
   file with force field parameters can be given in the second line after the
   molecule name and in this case connectivity is deduced from the bond
   lengths in the force field.

   There are many tools ([Open Babel](http://openbabel.org/),
   [Avogadro](http://avogadro.cc/), [VESTA](http://jp-minerals.org/vesta/en/))
   to create MDL mol files, xyz files or z-matrices. Manual editing of the
   files is usually necessary in order to match the atom names with those of
   the force field.

2. Use the `fftool` script to create an input file for `packmol`, which will
   use new `_pack.xyz` files with atomic coordinates for the components of
   your system. For help type `fftool -h`. For example, to build a simulation
   box with 40 ethanol and 300 water molecules and a density of 38.0 mol/L do:

        fftool 40 ethanol.zmat 300 spce.zmat -r 38.0

    Alternatively the side length of the the simulation box (here cubic) can
    be supplied in angstroms:

        fftool 40 ethanol.zmat 300 spce.zmat -b 20.0

3. Use `packmol` with the `pack.inp` file just created to generate the
   atomic coordinates in the simulation box:

        packmol < pack.inp

   Difficult convergence may indicate that density is too high, so adjust
   density or box size if necessary. For more complex spatial arrangements of
   molecules and materials you can modify the `pack.inp` to suit your needs
   (see the [Packmol](http://www.ime.unicamp.br/~martinez/packmol/)
   documentation).  Atomic coordinates for the full system will be written to
   `simbox.xyz`.

4. Use `fftool` to build the input files for LAMMPS (-l), DL_POLY (-d) or
   GROMACS (-g) containing the force field parameters and the coordinates of
   all the atoms (from `simbox.xyz`):

        fftool 40 ethanol.zmat 300 spce.zmat -r 38.0 -l

    If no force field information was given explicitly in the molecule files,
    a default LJ potential with parameters zeroed will be assigned to
    atoms. No terms for bonds, angles or torsions will be created. This is
    suitable when working with non-additive, bond-order or other potentials
    often used for materials. The input files for MD simulations will have to
    be edited manually to include an interaction potential for the material.


Deducing Bonds and Angles
-------------------------

When inferring connectivity from interatomic distances, distances in the
coordinates file are compared with equilibrium distances specified for bonds
in the force field and a tolerance of +/-0.25 angstrom is used to decide if a
bond should be present or not. So, the bond lengths in the conformation used
as input must be sufficiently close to those in the force field specification
for those bonds to be included in the potential energy fonction of the system.

Angles will be assigned to groups of three atoms i-j-k, with i-j and j-k
bonded, if the value of the angle in the conformation used as input is within
+/-15 degrees of the equilibrium angle in the force field specification. If
not, even if the atoms i-j-k are bonded, their angle will not be present in
the final potential energy function, although topologically the angle is
there. When running `fftool` to create a force field file (with `-l`, `-d` or
` -g` option) a warning message will show which such topological angles have
been "removed" because they deviate too much from the equilibrium angles in
the force field. This removal of angles avoids problems with atoms that have
more than four ligands, such as S or P atoms with five or six ligands. Around
these centers there are topological angles of 180 degrees to which no
potential energy of bending is attributed in force fields. For example, in the
octahedral PF6- anion there are two different values of F-P-F angles: twelve
90 degree angles between adjacent F atoms, and three 180 degree angles between
opposite F atoms; only the twelve 90 degree angles contribute with a harmonic
potential energy function in most force fields.

The tolerances for bond distances and angle values, 0.25 angstrom and 15
degrees, respectively, were chosen based on judgement. They can be set by
editing the `fftool` source, namely the global variables `BondTol` and
`AngleTol`. Use with care because spurious bonds and angles may be created if
the tolerances are too large.


Improper Dihedrals
------------------

Improper dihedrals are often used to increase the rigidity of planar atoms
(sp2) and differ from proper dihedrals in how they are defined. A proper
dihedral i-j-k-l is defined between bonded atoms i-j, j-k, and k-l and
corresponds to torsion around bond j-k, the dihedral being the angle between
planes i-j-k and j-k-l. An improper dihedral i-j-k-l is defined between bonded
atoms i-k, j-k and k-l, therefore k is a central atom bonded to the other
three. The central atom of the improper dihedral is assumed to be the third in
the list. Often in force fields the same potential energy function is used
both for proper and improper torsions.

If `improper` records are supplied in a molecule file (in `.zmat` format) then
those improper dihedrals are assumed by `fftool`. Otherwise, the script will
search for candidate improper dihedrals on all atoms with three bonds, with
any of `.zmat`, `.mol`, `.pdb` or `.xyz` input formats. A number of warning
messages will be printed if there are atoms with three bonds, which can be
ignored if the atoms in question are not centers of improper torsions. The
number and order of the atoms in the true improper dihedrals should be
checked in the files created.


Periodic Boundary Conditions
----------------------------

In molecular systems the initial configuration will generaly not contain
molecules crossing boundaries of the simulation box. A buffer distance of 1.5
angstrom is reserved at the box boundaries to avoid overlap of molecules from
periodic images in the initial configuration, as explained in the `packmol`
documentation (this empty space is added by `fftool` only for orthogonal
boxes). So the user should allow for this empty volume when supplying the size
of the box.

For simulations with extended materials it is possible to create chemical
bonds across boundaries. The option `-p` allows specification of periodic
conditions along x, y, z or combinations thereof. It is important in this case
to supply dimensions for the simulation box using the option `-b <l>` for a
cubic box, or `-b <lx,ly,lz>` for a general orthogonal box, or `-b
<a,b,c,alpha,beta,gamma>` for a general parallelepiped (triclinic box). An
energy minimization step prior to the start of the MD simulation is highly
recommended because no extra space is left near the boundaries and certain
molecules may overlap with those of neighboring images.

The coordinates of the atoms of the material have to be supplied in `.xyz`
format and prepared carefully so that distances across periodic boundaries are
within the tolerance to identify bonds. The number of bonds in the output
files created should be checked.

It is important that only the material for which bonds are to be established
across boudaries is supplied in `.xyz` format. The initial files for other
molecules in the system should be in `.zmat` or `.mol` formats, which contain
connectivity information. This is to avoid spurious bonds between atoms of the
molecular species that happen to be positioned too close to boundaries.

The `pack.inp` file will likely need manual editing in order to position the
atoms of the material precisely.


Force Field File Format
-----------------------

The `fftool` script reads a database of molecular force field terms in the
format described below. See the `examples` directory.

Blank lines and lines starting with `#` are ignored.

There are five sections, with headings `ATOMS`, `BONDS`, `ANGLES`, `DIHEDRALS`
and `IMPROPER`. Under each section heading, registers concerning the different
types of term in the force field are given.

`ATOMS` records describe, for each type of atom:
- the non-bonded atom type used for intermolecular interactions (these
  types may differ in the charges or intermolecular potential
  parameters)
- the bonded atom type used in intermolecular interactions (these
  types determine the intramolecular terms such as bonds, angles
  dihedrals)
- the mass in atomic units
- the electrostatic charge in elementary units
- the non-bonded potential type, e.g. `lj`
- potential parameters, namely Lennard-Jones `sigma` and `epsilon`

        C3H   CT  12.011  -0.18   lj    3.50   0.27614

`BONDS` records describe covalent bonds between intramolecular atom types:
- two bonded atom types
- type of bond potential, e.g. `harm` for harmonic potential or
  `cons` for a constrained bond.
- bond potential parameters, namely euqilibrium distance and force
  constant (the latter in the form k/2 (x - x0)^2)

        CT  CT   harm   1.529   2242.6

`ANGLES` records describe valence angles between intramolecular atom types:
- three bonded atom types, in which the central atom is bonded to the other
  two, e.g. i-j and j-k are bonded.
- type of angle potential, e.g. `harm` for harmonic potential or
  `cons` for a constrained angle.
- angle potential parameters, namely equilibrium angle and force
  constant (the latter in the form k/2 (x - x0)^2

        HC  CT  CT   harm   110.7   313.8

`DIHEDRALS` records describe torsion angles between intramolecular
atom types:
- four bonded atom types, in which atoms i-j, j-k, k-l are bonded.
- type of dihedral potential, e.g. `opls` for OPLS cosine series with
  four terms.
- dihedral potential parameters, with the coefficients in the form V_n/2
  (1 +/- cos(n phi)).

        CT  CT  CT  CT   opls    5.4392   -0.2092    0.8368    0.0000

`IMPROPER` records describe improper dihedral angles between
intramolecular atom types:
- four bonded atom types, in which atoms i-k, j-k, k-l are bonded.
- type of dihedral potential, e.g. `opls` for OPLS cosine series with
  four terms.
- dihedral potential parameters, with the coefficients in the form
  V_n/2 (1 +/- cos(n phi)).

        CA  CA  CA  HA   opls    0.0000    9.2048    0.0000    0.0000


Drude Polarizable Force Field
----------------------------
The `fftool` script can read a database of drude polarizable terms.
`examples/drude.inp` gives an example of the drude database.
The polarizable force field and input files can be exported for
LAMMPS, GROMACS and NAMD by `--drude` argument

        fftool 100 co2.zmat -b 30.0 --drude drude.inp -l
        fftool 100 co2.zmat -b 30.0 --drude drude.inp -n
        fftool 100 co2.zmat -b 30.0 --drude drude.inp -g
        
One should be careful that the implementation
of drude model are more or less different in these three software.
- GROMACS use SCF approach to solve the positions of drude particles
at each simulation step, whereas NAMD and LAMMPS use dual thermostats
(a hot one for normal atoms and a code one for drude particles)
to approximate the SCF results.
- Another difference is the handling of thole screening between
drude dipoles. LAMMPS applies thole screening to **EVERY** pairs of
drude dipoles, no matter intra- or inter-molecular, regardless of their
exclusion state defined by *special_bonds*.
NAMD by default only applies thole screening to 1-2 and 1-3 excluded
drude dipole pairs. But there's option to apply thole screening
to all pairs. GROMACS **CAN ONLY** apply thole screening to excluded
intra-molecular drude dipole pairs. Intermolecular thole screening is
not supported for GROMACS.
- The third difference is the handing of polarization catastrophe.
NAMD supports hardwall restraints and quadratic restraints.
GROMACS only supports quadratic restraints.
In LAMMPS, drude-core bonds are just normal bonds, so one
can apply any restraints supported by bond_style.

Although these differences, for a **REASONABLE** set of force field
parameters, all these three software should give *identical* results,
of course when applied with the same thermostat et. al.


Empirical Lennard-Jones parameters scaling
------------------------------
LJ parameters for non-polarizable force field implicitly include
mean contributions from polarizations. Introducing polarizability
directly into a well-developed non-polarizable force field will
cause the LJ interactions being overestimated. A method is proposed
to empirically scale down the LJ strength based on the charge,
dipole and polarizability of molecular fragments.
Read [the article](https://pubs.acs.org/doi/10.1021/acs.jctc.9b00689) for details.

The `fftool` script can read a database of properties of  monomers and dimers,
and calculate the scaling factor for LJ epsilon parameters between atom pairs.
`examples/ljscale.inp` gives an example of the LJ scaling database.

To apply the LJ scaling, use `--scaleeps` argument.

        fftool 100 DMSO.zmat -b 30.0 --drude drude.inp --scaleeps ljscale.inp -l
        fftool 100 DMSO.zmat -b 30.0 --drude drude.inp --scaleeps ljscale.inp -g
        fftool 100 DMSO.zmat -b 30.0 --drude drude.inp --scaleeps ljscale.inp -n

LJ sigma parameters can also be scaled by using `--scalesig` argument.
It can be used with or without `--scaleeps` argument.

        fftool 100 DMSO.zmat -b 30.0 --drude drude.inp --scaleeps ljscale.inp --scalesig 0.985 -l
        fftool 100 DMSO.zmat -b 30.0 --drude drude.inp --scalesig 1.1 -g
        fftool 100 DMSO.zmat -b 30.0 --scalesig 1.1 -n

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

* [GROMACS](http://www.gromacs.org/): H.J.C. Berendsen et. al.
  Comp Phys Commun, 91 (1995) 43,
  DOI: [10.1016/0010-4655(95)00042-E](https://doi.org/10.1016/0010-4655(95)00042-E)

* [NAMD](https://doi.org/10.1002/jcc.20289) James C Phillips et. al.
  J Comp Chem, 26 (2005) 1781–1802,
  DOI: [10.1002/jcc.20289](https://doi.org/10.1002/jcc.20289)

* [Empirical LJ scaling](https://pubs.acs.org/doi/10.1021/acs.jctc.9b00689)
  Kateryna Goloviznina et. al. J Chem Theory and Comput (2019),
  DOI: [10.1021/acs.jctc.9b00689](https://doi.org/10.1021/acs.jctc.9b00689)