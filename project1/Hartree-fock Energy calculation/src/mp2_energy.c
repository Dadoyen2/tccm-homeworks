// File: src/mp2_energy.c

#include <stdio.h>
#include <stdlib.h>
#include "mp2_energy.h"

// Macro for 4D indexing into a flattened 4D array
#define MO_TEI(i,j,k,l, arr, mo_num) \
    (arr[ (((size_t)(i)*(mo_num) + (j))*(mo_num) + (k))*(mo_num) + (l) ])

/**
 * @brief Computes the MP2 correlation energy using the provided molecular orbital energies and two-electron integrals.
 */
double compute_MP2_energy(double* mo_energy,
                          int mo_num,
                          int n_occ,
                          int64_t n_integrals,
                          int32_t* index,
                          double* value) {
    // Allocate memory for full two-electron integrals (MO TEI) using 8-fold symmetry
    size_t total_size = (size_t)mo_num * mo_num * mo_num * mo_num;
    double* mo_tei = (double*)calloc(total_size, sizeof(double));
    if (!mo_tei) {
        fprintf(stderr, "Memory allocation failed for MO two-electron integrals in MP2.\n");
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

    double emp2 = 0.0; // Initialize MP2 energy

    // Loop over occupied orbitals
    for (int i = 0; i < n_occ; i++) {
        for (int j = 0; j < n_occ; j++) {
            // Loop over virtual (unoccupied) orbitals
            for (int a = n_occ; a < mo_num; a++) {
                for (int b = n_occ; b < mo_num; b++) {
                    // Energy denominator: e_i + e_j - e_a - e_b
                    double denom = (mo_energy[i] + mo_energy[j]) - (mo_energy[a] + mo_energy[b]);

                    // Two-electron integrals: <ij|ab> and <ij|ba>
                    double ijab = MO_TEI(i, j, a, b, mo_tei, mo_num);
                    double ijba = MO_TEI(i, j, b, a, mo_tei, mo_num);

                    // Contribution to MP2 energy
                    double numerator = ijab * (2.0 * ijab - ijba);
                    emp2 += (numerator / denom);
                }
            }
        }
    }

    // Free allocated memory
    free(mo_tei);

    return emp2;
}

