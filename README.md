<img src="https://raw.githubusercontent.com/zorkzou/MolBO/master/molbo-logo.png" />

# MolBO
Generate NBO-47 file from the output file of MOLPRO.

## Latest Version
Version 2.1.5 (07/04/2014).

1. Some wrong comments in the source code have been corrected.
2. The format of geometry is a little adjusted for NBO 6.

## Features

1. It reads basis sets, overlap matrix, Fork matrix, density matrix, ... from [MOLPRO](http://www.molpro.net/)'s output file (use GPRINT and MATROP commands to obtain these data; see the examples), and prepares [NBO](http://nbo.chem.wisc.edu/)'s *.47 file.
2. It calculates Mayer's bond orders (MBO).
3. Electronic densities can be calculated at the following levels of theory by [MOLPRO](http://www.molpro.net/): SCF (RHF, ROHF, UHF, RDFT, RODFT, UDFT; with or without density fitting), Post-HF (MP2, MP3, CCSD, QCISD, QCISD(T), EOM-CCSD, Full-CI), Local Post-HF (LMP2), multi-configurational methods (MCSCF or CASSCF, CASVB, MRCI including SRCI, CASPT2 including SRMP2).
4. [MOLPRO](http://www.molpro.net/) 2006 - 2010 and NBO [3.1](http://www.ccl.net/cca/software/SOURCES/FORTRAN/nbo/index.shtml), 5.x, and [6.x](http://nbo.chem.wisc.edu/) have been tested.

## Limitations

1. Symmetry equivalent atoms are not allowed. It is necessary to decrease the symmetry until all the atoms are symmetry unique. For example, use C2v for CO2, Cs for H2O, and C1 for CF4. If you want to use high symmetry in your MOLPRO calculation, please save a Molden file and try the [Molden2AIM](https://github.com/zorkzou/Molden2AIM) program.
2. For a restart calculation of [MOLPRO](http://www.molpro.net/), do not use the permanent file 1.
3. For the all-electron relativistic calculations (e.g. DKHn), do not print kinetic energy and potential energy matrices because they are non-relativistic.

See [readme.html](https://zorkzou.github.io/MolBO/readme.html) for details.

See [NBO6](http://nbo.chem.wisc.edu/) webpage for better support of higher versions of [MOLPRO](http://www.molpro.net/): http://nbo6.chem.wisc.edu/INSTALL.molpro

Examples of applications can be found in W. Zou, D. Nori-Shargh, and J. E. Boggs, On the Covalent Character of Rare Gas Bonding Interactions: A New Kind of Weak Interaction, [J. Phys. Chem. A, 2013, 117, 207-212](http://pubs.acs.org/doi/abs/10.1021/jp3104535).
