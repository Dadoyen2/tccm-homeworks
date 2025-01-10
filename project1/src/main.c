#include <stdio.h>
#include <stdlib.h>
#include <trexio.h>

#include "hf_energy.h"
#include "mp2_energy.h"  // if you implement MP2

int main(int argc, char* argv[]) {

    // 1) Check command-line arguments
    if (argc != 2) {
        fprintf(stderr, "Usage: %s <trexio_filename>\n", argv[0]);
        return 1;
    }

    const char* filename = argv[1];
    trexio_exit_code rc;

    // 2) Open the TREXIO file
    trexio_t* trexio_file = trexio_open(filename, 'r', TREXIO_AUTO, &rc);
    if (rc != TREXIO_SUCCESS) {
        fprintf(stderr, "TREXIO Error: %s\n", trexio_string_of_error(rc));
        return 1;
    }

    // 3) Read Nuclear Repulsion
    double E_NN = read_nuclear_repulsion(trexio_file);
    printf("Nuclear repulsion energy = %.6f\n", E_NN);

    // 4) Read Number of Occupied Orbitals
    int n_occ = read_number_of_occupied_orbitals(trexio_file);
    printf("Number of occupied orbitals (n_occ) = %d\n", n_occ);

    // 5) Read Total Number of MOs
    int32_t mo_num;
    rc = trexio_read_mo_num(trexio_file, &mo_num);
    if (rc != TREXIO_SUCCESS) {
        fprintf(stderr, "TREXIO Error reading mo_num: %s\n", trexio_string_of_error(rc));
        trexio_close(trexio_file);
        return 1;
    }
    printf("Number of molecular orbitals (mo_num) = %d\n", mo_num);

    // 6) Read One-Electron Integrals
    double* one_e_integrals = read_one_electron_integrals(trexio_file, mo_num);

    // 7) Read Two-Electron Integrals
    int64_t n_integrals;
    int32_t* index;
    double*  value;
    read_two_electron_integrals(trexio_file, &n_integrals, &index, &value);
    printf("Number of non-zero two-electron integrals = %ld\n", (long)n_integrals);

    // 8) Compute HF energy
    double hf_energy = compute_HF_energy(E_NN, 
                                         one_e_integrals, 
                                         n_integrals, 
                                         index, 
                                         value, 
                                         mo_num, 
                                         n_occ);
    printf("Computed HF energy = %.8f\n", hf_energy);

    // 9) (Optional) Compute MP2 energy if implemented
    // double mp2_energy = compute_MP2_energy(/* pass needed arguments */);
    // printf("Computed MP2 energy correction = %.8f\n", mp2_energy);
    // printf("Total MP2 energy = %.8f\n", hf_energy + mp2_energy);

    // 10) Clean up
    free(one_e_integrals);
    free(index);
    free(value);

    rc = trexio_close(trexio_file);
    if (rc != TREXIO_SUCCESS) {
        fprintf(stderr, "TREXIO Error: %s\n", trexio_string_of_error(rc));
        return 1;
    }

    return 0;
}

