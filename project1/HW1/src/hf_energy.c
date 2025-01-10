// File: src/hf_energy.c

#include <stdio.h>
#include <stdlib.h>
#include <trexio.h>
#include "hf_energy.h"

// Macro for 4D indexing into a flattened 4D array
#define MO_TEI(i,j,k,l, arr, mo_num) \
    (arr[ (((size_t)(i)*(mo_num) + (j))*(mo_num) + (k))*(mo_num) + (l) ])

/**
 * @brief Reads the nuclear repulsion energy from a TREXIO file.
 */
double read_nuclear_repulsion(trexio_t* trexio_file) {
    trexio_exit_code rc;
    double E_NN = 0.0;

    rc = trexio_read_nucleus_repulsion(trexio_file, &E_NN);
    if (rc != TREXIO_SUCCESS) {
        fprintf(stderr, "TREXIO Error reading nuclear repulsion: %s\n",
                trexio_string_of_error(rc));
        exit(EXIT_FAILURE);
    }

    return E_NN;
}

/**
 * @brief Reads the number of occupied orbitals (n_occ) from a TREXIO file.
 */
int read_number_of_occupied_orbitals(trexio_t* trexio_file) {
    trexio_exit_code rc;
    int32_t n_up = 0;

    rc = trexio_read_electron_up_num(trexio_file, &n_up);
    if (rc != TREXIO_SUCCESS) {
        fprintf(stderr, "TREXIO Error reading number of up-spin electrons: %s\n",
                trexio_string_of_error(rc));
        exit(EXIT_FAILURE);
    }

    return (int)n_up; // For closed-shell systems, n_occ = n_up
}

/**
 * @brief Reads one-electron integrals (core Hamiltonian) from a TREXIO file.
 */
double* read_one_electron_integrals(trexio_t* trexio_file, int mo_num) {
    trexio_exit_code rc;
    // Allocate memory for mo_num x mo_num integrals
    double* integrals = (double*)malloc(mo_num * mo_num * sizeof(double));
    if (!integrals) {
        fprintf(stderr, "Memory allocation failed for one-electron integrals.\n");
        exit(EXIT_FAILURE);
    }

    rc = trexio_read_mo_1e_int_core_hamiltonian(trexio_file, integrals);
    if (rc != TREXIO_SUCCESS) {
        fprintf(stderr, "TREXIO Error reading one-electron integrals: %s\n",
                trexio_string_of_error(rc));
        free(integrals);
        exit(EXIT_FAILURE);
    }

    return integrals;
}

/**
 * @brief Reads two-electron integrals in sparse format from a TREXIO file.
 */
void read_two_electron_integrals(trexio_t* trexio_file,
                                 int64_t* n_integrals,
                                 int32_t** index,
                                 double** value) {
    trexio_exit_code rc;

    // Read the number of non-zero two-electron integrals
    rc = trexio_read_mo_2e_int_eri_size(trexio_file, n_integrals);
    if (rc != TREXIO_SUCCESS) {
        fprintf(stderr, "TREXIO Error reading number of two-electron integrals: %s\n",
                trexio_string_of_error(rc));
        exit(EXIT_FAILURE);
    }

    // Allocate memory for indices and values
    *index = (int32_t*)malloc(4 * (*n_integrals) * sizeof(int32_t));
    if (!(*index)) {
        fprintf(stderr, "Memory allocation failed for two-electron integrals indices.\n");
        exit(EXIT_FAILURE);
    }

    *value = (double*)malloc((*n_integrals) * sizeof(double));
    if (!(*value)) {
        fprintf(stderr, "Memory allocation failed for two-electron integrals values.\n");
        free(*index);
        exit(EXIT_FAILURE);
    }

    // Read the two-electron integrals
    int64_t buffer_size = *n_integrals;
    rc = trexio_read_mo_2e_int_eri(trexio_file, 0, &buffer_size, *index, *value);
    if (rc != TREXIO_SUCCESS) {
        fprintf(stderr, "TREXIO Error reading two-electron integrals: %s\n",
                trexio_string_of_error(rc));
        free(*index);
        free(*value);
        exit(EXIT_FAILURE);
    }

    // Verify that all integrals were read
    if (buffer_size != *n_integrals) {
        fprintf(stderr, "Mismatch in the number of two-electron integrals read.\n");
        free(*index);
        free(*value);
        exit(EXIT_FAILURE);
    }
}

/**
 * @brief Reads molecular orbital energies from a TREXIO file.
 */
trexio_exit_code read_mo_energies(trexio_t* trexio_file,
                                  int mo_num,
                                  double* mo_energy) {
    // Read molecular orbital energies
    trexio_exit_code rc = trexio_read_mo_energy(trexio_file, mo_energy);
    if (rc != TREXIO_SUCCESS) {
        fprintf(stderr, "TREXIO Error reading molecular orbital energies: %s\n",
                trexio_string_of_error(rc));
    }
    return rc;
}

/**
 * @brief Computes the Hartree-Fock energy using the provided integrals.
 */
double compute_HF_energy(double E_NN,
                         double* one_e_integrals,
                         int64_t n_integrals,
                         int32_t* index,
                         double* value,
                         int mo_num,
                         int n_occ) {
    double hf_energy = E_NN; // Start with nuclear repulsion energy

    // Add kinetic and electron-nucleus potential terms: 2 * sum_{i in occ} h_{ii}
    double sum_one_e = 0.0;
    for (int i = 0; i < n_occ; i++) {
        sum_one_e += 2.0 * one_e_integrals[i * mo_num + i];
    }
    hf_energy += sum_one_e;

    // Allocate memory for full two-electron integrals (MO TEI) using 8-fold symmetry
    size_t total_size = (size_t)mo_num * mo_num * mo_num * mo_num;
    double* mo_tei = (double*)calloc(total_size, sizeof(double));
    if (!mo_tei) {
        fprintf(stderr, "Memory allocation failed for MO two-electron integrals.\n");
        exit(EXIT_FAILURE);
    }

    // Populate the MO TEI array using the sparse integrals and symmetry
    for (int64_t n = 0; n < n_integrals; n++) {
        int i = index[4 * n + 0];
        int j = index[4 * n + 1];
        int k = index[4 * n + 2];
        int l = index[4 * n + 3];
        double val = value[n];

        // Applying 8-fold permutational symmetry
        MO_TEI(i, j, k, l, mo_tei, mo_num) = val;
        MO_TEI(i, l, k, j, mo_tei, mo_num) = val;
        MO_TEI(k, l, i, j, mo_tei, mo_num) = val;
        MO_TEI(k, j, i, l, mo_tei, mo_num) = val;
        MO_TEI(j, i, l, k, mo_tei, mo_num) = val;
        MO_TEI(l, i, j, k, mo_tei, mo_num) = val;
        MO_TEI(l, k, j, i, mo_tei, mo_num) = val;
        MO_TEI(j, k, l, i, mo_tei, mo_num) = val;
    }

    // Add two-electron integrals: sum_{i,j in occ} [2 <ij|ij> - <ij|ji>]
    double sum_two_e = 0.0;
    for (int i = 0; i < n_occ; i++) {
        for (int j = 0; j < n_occ; j++) {
            double coulomb = MO_TEI(i, j, i, j, mo_tei, mo_num);
            double exchange = MO_TEI(i, j, j, i, mo_tei, mo_num);
            sum_two_e += (2.0 * coulomb - exchange);
        }
    }
    hf_energy += sum_two_e;

    // Free allocated memory
    free(mo_tei);

    return hf_energy;
}

