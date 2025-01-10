// File: src/main.c

#include <stdio.h>
#include <stdlib.h>
#include <trexio.h>
#include "hf_energy.h"
#include "mp2_energy.h"

int main(int argc, char* argv[]) {
    // Check for correct usage
    if (argc != 2) {
        fprintf(stderr, "Usage: %s <trexio_file>\n", argv[0]);
        return EXIT_FAILURE;
    }

    // 1. Open TREXIO file
    const char* filename = argv[1];
    trexio_exit_code rc;
    trexio_t* trexio_file = trexio_open(filename, 'r', TREXIO_AUTO, &rc);
    if (rc != TREXIO_SUCCESS) {
        fprintf(stderr, "TREXIO Error opening file '%s': %s\n", filename, trexio_string_of_error(rc));
        return EXIT_FAILURE;
    }

    // 2. Read nuclear repulsion energy
    double E_NN = read_nuclear_repulsion(trexio_file);
    printf("Nuclear repulsion energy (E_NN) = %.6f atomic units\n", E_NN);

    // 3. Read number of occupied orbitals
    int n_occ = read_number_of_occupied_orbitals(trexio_file);
    printf("Number of occupied orbitals (n_occ) = %d\n", n_occ);

    // 4. Read total number of molecular orbitals (mo_num)
    int32_t mo_num;
    rc = trexio_read_mo_num(trexio_file, &mo_num);
    if (rc != TREXIO_SUCCESS) {
        fprintf(stderr, "TREXIO Error reading number of MOs (mo_num): %s\n",
                trexio_string_of_error(rc));
        trexio_close(trexio_file);
        return EXIT_FAILURE;
    }
    printf("Number of molecular orbitals (mo_num) = %d\n", mo_num);

    // 5. Read one-electron integrals
    double* one_e_integrals = read_one_electron_integrals(trexio_file, mo_num);

    // 6. Read two-electron integrals
    int64_t n_integrals;
    int32_t* index;
    double* value;
    read_two_electron_integrals(trexio_file, &n_integrals, &index, &value);
    printf("Number of non-zero two-electron integrals = %ld\n", (long)n_integrals);

    // 7. Compute Hartree-Fock energy
    double hf_energy = compute_HF_energy(E_NN,
                                         one_e_integrals,
                                         n_integrals,
                                         index,
                                         value,
                                         mo_num,
                                         n_occ);
    printf("Computed Hartree-Fock energy (E_HF) = %.8f atomic units\n", hf_energy);

    // 8. Read molecular orbital energies
    double* mo_energy = (double*)malloc(mo_num * sizeof(double));
    if (!mo_energy) {
        fprintf(stderr, "Memory allocation failed for molecular orbital energies.\n");
        // Free previously allocated memory before exiting
        free(one_e_integrals);
        free(index);
        free(value);
        trexio_close(trexio_file);
        return EXIT_FAILURE;
    }
    rc = read_mo_energies(trexio_file, mo_num, mo_energy);
    if (rc != TREXIO_SUCCESS) {
        trexio_close(trexio_file);
        free(one_e_integrals);
        free(index);
        free(value);
        free(mo_energy);
        return EXIT_FAILURE;
    }

    // 9. Compute MP2 correlation energy
    double mp2_energy = compute_MP2_energy(mo_energy,
                                           mo_num,
                                           n_occ,
                                           n_integrals,
                                           index,
                                           value);
    printf("Computed MP2 correlation energy (EMP2) = %.8f atomic units\n", mp2_energy);

    // 10. Print total MP2 energy (E_HF + EMP2)
    printf("Total MP2 energy (E_HF + EMP2) = %.8f atomic units\n", hf_energy + mp2_energy);

    // Cleanup: Free allocated memory and close TREXIO file
    free(one_e_integrals);
    free(index);
    free(value);
    free(mo_energy);

    rc = trexio_close(trexio_file);
    if (rc != TREXIO_SUCCESS) {
        fprintf(stderr, "TREXIO Error closing file '%s': %s\n", filename, trexio_string_of_error(rc));
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}

