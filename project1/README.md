<<<<<<< HEAD
Summary

A C-based program to calculate the Møller-Plesset second-order perturbation theory (MP2) correlation energy and Hartree-Fock (HF) energy for a closed-shell molecular system can be found in this directory. The code reads integral data, including orbital energies and one- and two-electron integrals, from a TREXIO file and outputs the total HF and MP2 energies for a specific molecule.

Directory Structure:

├── LICENSE                # It Holds the license text under which this software is released.
├── AUTHORS                # It contains lists oof the contributors to this project, with affiliations.
├── README.md              # This file provides an overview of the project, including objectives and directory structure.
├── INSTALL.md             # It provides step-by-step instructions for installing the required libraries and for compiling the code.
├── data                   # Stores the input .h5 files that contain one-electron and two-electron integrals, nuclear repulsion energy and orbital energies.
│   ├── c2h2.h5
│   ├── ch4.h5
│   ├── co2.h5
│   ├── h2o.h5
│   ├── h3coh.h5
│   └── hcn.h5             # Example TREXIO/HDF5 data files.
├── tests                  # The main C source code for this assignment is here. Although the folder’s name is “tests,” we can also store our core code here|                            The final compiled program is generated as compute_energy.
│   ├── hf_energy.h        # Header for functions handling HF integrals and energy calculations.
│   ├── hf_energy.c        # Source code for reading TREXIO data and computing Hartree-Fock energies.
│   ├── mp2_energy.h       # Header for functions computing MP2 correlation energy.
│   ├── mp2_energy.c       # Source code implementing the MP2 energy formula.
│   ├── main.c             # This is the main program orchestrating reading integrals, executing HF/MP2 calculations, and printing results.
│   └── compute_energy     # Compiled executable (produced after running 'make').
├── src
│   └── trexio-2.5.0       # TREXIO library source, placed here for local compilation and linking.
├── Makefile               # Automates compilation/linking with HDF5 and TREXIO libraries.

Project Description:

Hartree-Fock Computation:
The Hartree-Fock module reads the necessary integral data and orbital occupations from the TREXIO file, then applies the standard HF equations to produce a self-consistent-field energy.

MP2 Computation:
Leveraging the converged Hartree-Fock orbitals, the MP2 module computes the second-order perturbation correction to the HF energy, thus yielding the total correlated MP2 energy.

Usage:
Compile using the provided Makefile (e.g., make).
Execute the resulting compute_energy program by passing in a valid TREXIO .h5 file containing the integrals and orbital energies (e.g., ./compute_energy data/h2o.h5).
Data Requirements:
The .h5 file should include one-electron integrals (core Hamiltonian), two-electron integrals in sparse format, nuclear repulsion energy, and orbital energies. These data structures must conform to the TREXIO format.
=======
# HF-and-MP2-
>>>>>>> 2143981678f912e57b0dec288b8ff3b0295272d4
