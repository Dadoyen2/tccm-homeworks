#include <stdio.h>
#include <stdlib.h>
#include <trexio.h>   // Make sure you have installed TREXIO properly

// Helper function to look up 2e integrals from (index,value). 
// For now, a simple brute-force approach:
double get_2e_integral(int p, int q, int r, int s,
                       const int32_t* idx, const double* val,
                       int64_t n_int)
{
    // Because of 8-fold symmetry, you must find the matching integral
    // in the (idx,val) arrays. This is just a placeholder for demonstration:
    for (int64_t n = 0; n < n_int; n++) {
        int i = idx[4*n + 0];
        int j = idx[4*n + 1];
        int k = idx[4*n + 2];
        int l = idx[4*n + 3];

        // If (i,j,k,l) matches (p,q,r,s) under the stored convention
        // or any recognized symmetry, return val[n].
        // For simplicity, let's do a direct match:
        if (i == p && j == q && k == r && l == s) {
            return val[n];
        }
        // You could also check permutations if needed.
    }
    // If not found, integral = 0.
    return 0.0;
}

int main(int argc, char* argv[]) {

    if (argc < 2) {
        fprintf(stderr, "Usage: %s <trexio_file>\n", argv[0]);
        return 1;
    }

    // 1. Open file
    trexio_exit_code rc;
    trexio_t* file = trexio_open(argv[1], 'r', TREXIO_AUTO, &rc);
    if (rc != TREXIO_SUCCESS) {
        fprintf(stderr, "Cannot open TREXIO file: %s\n", trexio_string_of_error(rc));
        return 1;
    }

    // 2. Read nuclear repulsion energy
    double E_nuc;
    rc = trexio_read_nucleus_repulsion(file, &E_nuc);
    if (rc != TREXIO_SUCCESS) {
        fprintf(stderr, "Error reading nuclear repulsion: %s\n", trexio_string_of_error(rc));
        return 1;
    }

    // 3. Read number of up electrons
    int32_t n_up;
    rc = trexio_read_electron_up_num(file, &n_up);
    if (rc != TREXIO_SUCCESS) {
        fprintf(stderr, "Error reading electron_up_num: %s\n", trexio_string_of_error(rc));
        return 1;
    }
    int n_occ = n_up;  // for closed-shell

    // 4. Read mo_num
    int32_t mo_num;
    rc = trexio_read_mo_num(file, &mo_num);
    if (rc != TREXIO_SUCCESS) {
        fprintf(stderr, "Error reading mo_num: %s\n", trexio_string_of_error(rc));
        return 1;
    }

    // 5. Allocate & read the 1-electron integrals
    double* Hcore = (double*)malloc(mo_num * mo_num * sizeof(double));
    if (Hcore == NULL) {
        fprintf(stderr, "Malloc failed for Hcore\n");
        return 1;
    }
    rc = trexio_read_mo_1e_int_core_hamiltonian(file, Hcore);
    if (rc != TREXIO_SUCCESS) {
        fprintf(stderr, "Error reading 1e ints: %s\n", trexio_string_of_error(rc));
        return 1;
    }

    // 6. Read the 2-electron integrals (sparse)
    int64_t n_integrals;
    rc = trexio_read_mo_2e_int_eri_size(file, &n_integrals);
    if (rc != TREXIO_SUCCESS) {
        fprintf(stderr, "Error reading 2e size: %s\n", trexio_string_of_error(rc));
        return 1;
    }

    int32_t* idx = (int32_t*)malloc(4 * n_integrals * sizeof(int32_t));
    double*  val = (double*) malloc(n_integrals * sizeof(double));
    if (idx == NULL || val == NULL) {
        fprintf(stderr, "Malloc failed for 2e integrals\n");
        return 1;
    }
    int64_t buffer_size = n_integrals;
    rc = trexio_read_mo_2e_int_eri(file, 0, &buffer_size, idx, val);
    if (rc != TREXIO_SUCCESS) {
        fprintf(stderr, "Error reading 2e ints: %s\n", trexio_string_of_error(rc));
        return 1;
    }

    // 7. Read orbital energies
    double* mo_energy = (double*)malloc(mo_num * sizeof(double));
    if (mo_energy == NULL) {
        fprintf(stderr, "Malloc failed for mo_energy\n");
        return 1;
    }
    rc = trexio_read_mo_energy(file, mo_energy);
    if (rc != TREXIO_SUCCESS) {
        fprintf(stderr, "Error reading mo_energy: %s\n", trexio_string_of_error(rc));
        return 1;
    }

    // 8. Close the file
    rc = trexio_close(file);
    if (rc != TREXIO_SUCCESS) {
        fprintf(stderr, "Error closing file: %s\n", trexio_string_of_error(rc));
        return 1;
    }
    file = NULL; // not strictly necessary, but good practice

    // ---- Now compute HF energy (very simplified) ----
    // HF = E_nuc
    //    + 2 * sum_{i=1 to n_occ} <i|h|i>
    //    + sum_{i=1}^{n_occ} sum_{j=1}^{n_occ} [2 <ij|ij> - <ij|ji> ]
    double E_HF = E_nuc;

    // One-electron sum
    double sum_one = 0.0;
    for (int i = 0; i < n_occ; i++) {
        sum_one += Hcore[i*mo_num + i]; // <i|h|i>
    }
    E_HF += 2.0 * sum_one;

    // Two-electron sum
    double sum_two = 0.0;
    for (int i = 0; i < n_occ; i++) {
      for (int j = 0; j < n_occ; j++) {
        double ijij = get_2e_integral(i, j, i, j, idx, val, n_integrals);
        double ijji = get_2e_integral(i, j, j, i, idx, val, n_integrals);
        sum_two += (2.0 * ijij - ijji);
      }
    }
    E_HF += sum_two;

    // ---- MP2 correlation ----
    // E_MP2 = sum_{ij in occ} sum_{ab in virt} ( <ij|ab> * ( <ij|ab> - <ij|ba> ) ) / (ei + ej - ea - eb)
    double E_MP2 = 0.0;
    for (int i = 0; i < n_occ; i++) {
      for (int j = 0; j < n_occ; j++) {
        for (int a = n_occ; a < mo_num; a++) {
          for (int b = n_occ; b < mo_num; b++) {
            double ijab = get_2e_integral(i, j, a, b, idx, val, n_integrals);
            double ijba = get_2e_integral(i, j, b, a, idx, val, n_integrals);

            double numerator  = ijab * (ijab - ijba);
            double denominator = mo_energy[i] + mo_energy[j]
                               - mo_energy[a] - mo_energy[b];

            E_MP2 += numerator / denominator;
          }
        }
      }
    }

    // Print results
    printf("HF energy = %20.8f\n", E_HF);
    printf("MP2 corr. = %20.8f\n", E_MP2);
    printf("MP2 total = %20.8f\n", E_HF + E_MP2);

    // Free memory
    free(Hcore);
    free(idx);
    free(val);
    free(mo_energy);

    return 0;
}

