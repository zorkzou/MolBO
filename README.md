<img src="https://github.com/zorkzou/MolBO/molbo-logo.png" />

# MolBO
Generate NBO-47 file from the output file of MOLPRO.

## Latest Version
Version 2.1.5 (07/04/2014).

1. Some wrong comments in the source code have been corrected.
2. The format of geometry is a little adjusted for NBO 6.

## Features

1. It reads basis sets, overlap matrix, Fork matrix, density matrix, ... from MOLPRO's output file (use GPRINT and MATROP commands to obtain these data; see the examples), and prepares NBO's *.47 file.
2. It calculates Mayer's bond orders (MBO).
3. Electronic densities can be calculated at the following levels of theory by MOLPRO: SCF (RHF, ROHF, UHF, RDFT, RODFT, UDFT; with or without density fitting), Post-HF (MP2, MP3, CCSD, QCISD, QCISD(T), EOM-CCSD, Full-CI), Local Post-HF (LMP2), multi-configurational methods (MCSCF or CASSCF, CASVB, MRCI including SRCI, CASPT2 including SRMP2).
4. MOLPRO 2006 - 2010 and NBO 3.1, 5.x, and 6.x have been tested.

See readme.html for details.

See NBO6 webpage for better support of higher versions of MOLPRO: http://nbo6.chem.wisc.edu/INSTALL.molpro
