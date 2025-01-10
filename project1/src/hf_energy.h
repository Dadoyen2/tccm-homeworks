#ifndef HF_ENERGY_H
#define HF_ENERGY_H

#include <trexio.h>

// Functions to read data from TREXIO
double read_nuclear_repulsion(trexio_t* trexio_file);
int    read_number_of_occupied_orbitals(trexio_t* trexio_file);
double* read_one_electron_integrals(trexio_t* trexio_file, int mo_num);
void   read_two_electron_integrals(trexio_t* trexio_file, 
                                   int64_t* n_integrals, 
                                   int32_t** index, 
                                   double** value);

// Function to compute HF energy
double compute_HF_energy(double E_NN, 
                         double* one_e_integrals, 
                         int64_t n_integrals, 
                         int32_t* index, 
                         double* value, 
                         int mo_num, 
                         int n_occ);

#endif // HF_ENERGY_H

