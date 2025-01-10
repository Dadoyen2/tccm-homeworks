#include <stdio.h>
#include <stdlib.h>
#include <trexio.h>
#include "hf_energy.h"

//---------------------------------------------------------------------
// 1) Reading nuclear repulsion energy
//---------------------------------------------------------------------
double read_nuclear_repulsion(trexio_t* trexio_file) {
    trexio_exit_code rc;
    double energy;
    
    rc = trexio_read_nucleus_repulsion(trexio_file, &energy);
    if (rc != TREXIO_SUCCESS) {
        fprintf(stderr, "TREXIO Error reading nuclear repulsion energy: %s\n",
                trexio_string_of_error(rc));
        exit(1);
    }
    return energy;
}

//---------------------------------------------------------------------
// 2) Reading number of occupied orbitals
//---------------------------------------------------------------------
int read_number_of_occupied_orbitals(trexio_t* trexio_file) {
    trexio_exit_code rc;
    int32_t n_up;
    
    rc = trexio_read_electron_up_num(trexio_file, &n_up);
    if (rc != TREXIO_SUCCESS) {
        fprintf(stderr, "TREXIO Error reading number of occupied orbitals: %s\n",
                trexio_string_of_error(rc));
        exit(1);
    }
    return n_up; // For closed-shell systems, Nocc = n_up
}

//---------------------------------------------------------------------
// 3) Reading one-electron integrals
//---------------------------------------------------------------------
double* read_one_electron_integrals(trexio_t* trexio_file, int mo_num) {
    trexio_exit_code rc;

    // Allocate memory for a mo_num x mo_num array
    double* one_e_integrals = (double*)malloc(mo_num * mo_num * sizeof(double));
    if (one_e_integrals == NULL) {
        fprintf(stderr, "Error: Unable to allocate memory for one-electron integrals.\n");
        exit(1);
    }

    // Read the integrals from the file
    rc = trexio_read_mo_1e_int_core_hamiltonian(trexio_file, one_e_integrals);
    if (rc != TREXIO_SUCCESS) {
        fprintf(stderr, "TREXIO Error reading one-electron integrals: %s\n",
                trexio_string_of_error(rc));
        free(one_e_integrals);
        exit(1);
    }

    return one_e_integrals;
}

//---------------------------------------------------------------------
// 4) Reading two-electron integrals
//---------------------------------------------------------------------
void read_two_electron_integrals(trexio_t* trexio_file, 
                                 int64_t* n_integrals, 
                                 int32_t** index, 
                                 double** value) 
{
    trexio_exit_code rc;

    // First read the number of non-zero two-electron integrals
    rc = trexio_read_mo_2e_int_eri_size(trexio_file, n_integrals);
    if (rc != TREXIO_SUCCESS) {
        fprintf(stderr, "TREXIO Error reading size of 2e integrals: %s\n",
                trexio_string_of_error(rc));
        exit(1);
    }

    // Allocate memory for indices and values
    *index = (int32_t*)malloc(4 * (*n_integrals) * sizeof(int32_t));
    *value = (double*)malloc((*n_integrals) * sizeof(double));

    if (*index == NULL || *value == NULL) {
        fprintf(stderr, "Error: Unable to allocate memory for two-electron integrals.\n");
        exit(1);
    }

    // Read the integrals from the file
    int64_t buffer_size = *n_integrals;
    rc = trexio_read_mo_2e_int_eri(trexio_file,
                                   0,                // offset_file
                                   &buffer_size,     // pointer to buffer_size
                                   *index,
                                   *value);

    if (rc != TREXIO_SUCCESS) {
        fprintf(stderr, "TREXIO Error reading 2e integrals: %s\n",
                trexio_string_of_error(rc));
        free(*index);
        free(*value);
        exit(1);
    }
}

//---------------------------------------------------------------------
// 5) Computing the Hartree-Fock energy
//---------------------------------------------------------------------
//
// E(HF) = E_NN
//         + 2 * \sum_{i=1 to n_occ} h_{ii}
//         + \sum_{i=1 to n_occ}\sum_{j=1 to n_occ} [2 <ij|ij> - <ij|ji>]
//
// Indices i,j run over occupied MOs
//---------------------------------------------------------------------

// A helper macro to do 4D indexing in a 1D array
#define MO_TEI(i,j,k,l, arr, mo_num) \
    (arr[ (((size_t)(i) * (mo_num) + (j)) * (mo_num) + (k)) * (mo_num) + (l) ])

double compute_HF_energy(double E_NN, 
                         double* one_e_integrals, 
                         int64_t n_integrals, 
                         int32_t* index, 
                         double* value, 
                         int mo_num, 
                         int n_occ) 
{
    //-----------------------------------------------------------------
    // 5.1) Start with nuclear repulsion
    //-----------------------------------------------------------------
    double hf_energy = E_NN;

    //-----------------------------------------------------------------
    // 5.2) Add the one-electron part: 2 * sum_{i=1..n_occ} h_{ii}
    //-----------------------------------------------------------------
    double sum_one_e = 0.0;
    for (int i = 0; i < n_occ; i++) {
        sum_one_e += 2.0 * one_e_integrals[i * mo_num + i];
    }
    hf_energy += sum_one_e;

    //-----------------------------------------------------------------
    // 5.3) Build a full 4D array of MO two-electron integrals
    //      (only needed for the occupied block, but we’ll store all).
    //-----------------------------------------------------------------
    size_t total_size = (size_t)mo_num * mo_num * mo_num * mo_num;
    double* mo_tei = (double*) calloc(total_size, sizeof(double));
    if (mo_tei == NULL) {
        fprintf(stderr, "Error: Unable to allocate memory for mo_tei.\n");
        exit(1);
    }

    // Fill the mo_tei array, using all possible permutations
    // because integrals obey 8-fold symmetry.
    for (int64_t n = 0; n < n_integrals; n++) {
        int i = index[4*n + 0];
        int j = index[4*n + 1];
        int k = index[4*n + 2];
        int l = index[4*n + 3];

        double integral = value[n];

        // Store the integral in mo_tei for all symmetric permutations:
        MO_TEI(i,j,k,l, mo_tei, mo_num) = integral;
        MO_TEI(i,l,k,j, mo_tei, mo_num) = integral;
        MO_TEI(k,l,i,j, mo_tei, mo_num) = integral;
        MO_TEI(k,j,i,l, mo_tei, mo_num) = integral;

        MO_TEI(j,i,l,k, mo_tei, mo_num) = integral;
        MO_TEI(l,i,j,k, mo_tei, mo_num) = integral;
        MO_TEI(l,k,j,i, mo_tei, mo_num) = integral;
        MO_TEI(j,k,l,i, mo_tei, mo_num) = integral;
    }

    //-----------------------------------------------------------------
    // 5.4) Compute the two-electron part:
    //      sum_{i=1..n_occ} sum_{j=1..n_occ} [2 <ij|ij> - <ij|ji>]
    //-----------------------------------------------------------------
    double two_e_sum = 0.0;
    for (int i = 0; i < n_occ; i++) {
        for (int j = 0; j < n_occ; j++) {
            // coulomb  = <ij|ij>
            // exchange = <ij|ji>
            double coulomb  = MO_TEI(i, j, i, j, mo_tei, mo_num);
            double exchange = MO_TEI(i, j, j, i, mo_tei, mo_num);
            two_e_sum += (2.0 * coulomb - exchange);
        }
    }

    hf_energy += two_e_sum;

    //-----------------------------------------------------------------
    // 5.5) Free allocated memory
    //-----------------------------------------------------------------
    free(mo_tei);

    //-----------------------------------------------------------------
    // 5.6) Return final HF energy
    //-----------------------------------------------------------------
    return hf_energy;
}

