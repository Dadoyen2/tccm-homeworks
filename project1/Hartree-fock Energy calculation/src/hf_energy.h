// File: src/hf_energy.h

#ifndef HF_ENERGY_H
#define HF_ENERGY_H

#include <trexio.h>

/**
 * @brief Reads the nuclear repulsion energy from a TREXIO file.
 *
 * @param trexio_file Pointer to an open TREXIO file.
 * @return Nuclear repulsion energy as a double.
 */
double read_nuclear_repulsion(trexio_t* trexio_file);

/**
 * @brief Reads the number of occupied orbitals (n_occ) from a TREXIO file.
 *
 * For a closed-shell system, n_occ is equal to the number of up-spin electrons.
 *
 * @param trexio_file Pointer to an open TREXIO file.
 * @return Number of occupied orbitals as an integer.
 */
int read_number_of_occupied_orbitals(trexio_t* trexio_file);

/**
 * @brief Reads one-electron integrals (core Hamiltonian) from a TREXIO file.
 *
 * The integrals are stored in a dynamically allocated array of size mo_num^2.
 *
 * @param trexio_file Pointer to an open TREXIO file.
 * @param mo_num Number of molecular orbitals.
 * @return Pointer to the array containing one-electron integrals.
 */
double* read_one_electron_integrals(trexio_t* trexio_file, int mo_num);

/**
 * @brief Reads two-electron integrals in sparse format from a TREXIO file.
 *
 * The integrals are stored in two arrays: index[] (contains indices) and value[] (contains integral values).
 *
 * @param trexio_file Pointer to an open TREXIO file.
 * @param n_integrals Pointer to store the number of non-zero integrals.
 * @param index Double pointer to store the indices array.
 * @param value Double pointer to store the values array.
 */
void read_two_electron_integrals(trexio_t* trexio_file,
                                 int64_t* n_integrals,
                                 int32_t** index,
                                 double** value);

/**
 * @brief Reads molecular orbital energies from a TREXIO file.
 *
 * The energies are stored in the mo_energy array of length mo_num.
 *
 * @param trexio_file Pointer to an open TREXIO file.
 * @param mo_num Number of molecular orbitals.
 * @param mo_energy Array to store molecular orbital energies.
 * @return TREXIO exit code indicating success or failure.
 */
trexio_exit_code read_mo_energies(trexio_t* trexio_file,
                                  int mo_num,
                                  double* mo_energy);

/**
 * @brief Computes the Hartree-Fock energy using integrals read from the TREXIO file.
 *
 * The formula used:
 * E(HF) = E_NN + 2 * sum_{i in occ} h_{ii} + sum_{i,j in occ} [2 <ij|ij> - <ij|ji>]
 *
 * @param E_NN Nuclear repulsion energy.
 * @param one_e_integrals Array of one-electron integrals.
 * @param n_integrals Number of non-zero two-electron integrals.
 * @param index Indices array for two-electron integrals.
 * @param value Values array for two-electron integrals.
 * @param mo_num Number of molecular orbitals.
 * @param n_occ Number of occupied orbitals.
 * @return Computed Hartree-Fock energy as a double.
 */
double compute_HF_energy(double E_NN,
                         double* one_e_integrals,
                         int64_t n_integrals,
                         int32_t* index,
                         double* value,
                         int mo_num,
                         int n_occ);

#endif // HF_ENERGY_H

