#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#define MAX_ATOMS 100
#define EPSILON 0.0661 // Lennard-Jones epsilon (in kcal/mol)
#define SIGMA 0.3345   // Lennard-Jones sigma (in nm)
#define TOTAL_STEPS 1000 // Total simulation steps
#define OUTPUT_INTERVAL 10 // Interval for writing output
#define TIMESTEP 0.2 // Time step for simulation

// Function to allocate a 2D array
double** malloc_2d(size_t m, size_t n) {
    double** a = malloc(m * sizeof(double*));
    if (a == NULL) {
        return NULL;
    }
    a[0] = malloc(m * n * sizeof(double));
    if (a[0] == NULL) {
        free(a);
        return NULL;
    }
    for (size_t i = 1; i < m; i++) {
        a[i] = a[i - 1] + n;
    }
    return a;
}

// Function to free a 2D array
void free_2d(double** a) {
    if (a) {
        free(a[0]);
        free(a);
    }
}

// Struct to store atom data
typedef struct {
    char atom[3];
    double x, y, z;
    double mass;
} Atom;

// Function to calculate distances
void compute_distances(size_t Natoms, double** coord, double** distance) {
    for (size_t i = 0; i < Natoms; i++) {
        for (size_t j = 0; j < Natoms; j++) {
            if (i == j) {
                distance[i][j] = 0.0;
            } else {
                double dx = (coord[i][0] - coord[j][0]);
                double dy = (coord[i][1] - coord[j][1]) ;
                double dz = (coord[i][2] - coord[j][2]) ;
                distance[i][j] = sqrt(dx * dx + dy * dy + dz * dz);
            }
        }
    }
}

// Function to calculate Lennard-Jones potential energy
double compute_LJ_potential(double epsilon, double sigma, size_t Natoms, double** distance) {
    double total_potential = 0.0;

    for (size_t i = 0; i < Natoms; i++) {
        for (size_t j = i + 1; j < Natoms; j++) { // Only consider j > i
            double r = distance[i][j];
            if (r > 0) { // Avoid division by zero
                double r_inv = sigma / r;
                double r_inv6 = pow(r_inv, 6);
                double r_inv12 = r_inv6 * r_inv6;
                total_potential += 4 * epsilon * (r_inv12 - r_inv6);
            }
        }
    }

    return total_potential;
}

// Function to compute kinetic energy
double compute_kinetic_energy(size_t Natoms, double** velocity, double* mass) {
    double total_kinetic = 0.0;

    for (size_t i = 0; i < Natoms; i++) {
        double v2 = velocity[i][0] * velocity[i][0] +
                    velocity[i][1] * velocity[i][1] +
                    velocity[i][2] * velocity[i][2]; // v^2 = vx^2 + vy^2 + vz^2
        total_kinetic += 0.5 * (mass[i]) * v2; // T = 1/2 * m * v^2
    }

    return total_kinetic;
}

// Function to compute accelerations based on Lennard-Jones forces
void compute_acc(size_t Natoms, double** coord, double* mass, double** distance, double** acceleration) {
    // Constants for the Lennard-Jones potential
    double epsilon = 0.0661; // J/mol
    double sigma = 0.3345;  // nm

    // Initialize acceleration array to zero
    for (size_t i = 0; i < Natoms; i++) {
        acceleration[i][0] = 0.0; // ax
        acceleration[i][1] = 0.0; // ay
        acceleration[i][2] = 0.0; // az
    }

    // Compute accelerations for each atom
    for (size_t i = 0; i < Natoms; i++) {
        for (size_t j = 0; j < Natoms; j++) {
            if (i != j) {
                double dx = (coord[i][0] - coord[j][0]) ;
                double dy = (coord[i][1] - coord[j][1]) ;
                double dz = (coord[i][2] - coord[j][2]) ;
                double r = distance[i][j];

                if (r > 0) { // Avoid division by zero
                    double r_inv = sigma / r;
                    double r_inv6 = pow(r_inv, 6);
                    double r_inv12 = r_inv6 * r_inv6;
                    double force_mag = 24 * epsilon * (2 * r_inv12 - r_inv6) / (r * r); // Force magnitude

                    // Update acceleration components
                    acceleration[i][0] += force_mag * dx / (mass[i] );
                    acceleration[i][1] += force_mag * dy / (mass[i] );
                    acceleration[i][2] += force_mag * dz / (mass[i] );
                }
            }
        }
    }
}

// Verlet integration for updating positions and velocities
void verlet_update(size_t Natoms, double dt, double** coord, double** velocity, double** acceleration, double** distance, double* mass) {
    double** new_acceleration = malloc_2d(Natoms, 3);

    if (new_acceleration == NULL) {
        printf("Memory allocation failed for new_acceleration\n");
        exit(EXIT_FAILURE);
    }

    // Update positions
    for (size_t i = 0; i < Natoms; i++) {
        for (size_t dim = 0; dim < 3; dim++) {
            coord[i][dim] += velocity[i][dim] * dt + 0.5 * acceleration[i][dim] * dt * dt;
        }
    }

    // Update distances and accelerations
    compute_distances(Natoms, coord, distance);
    compute_acc(Natoms, coord, mass, distance, new_acceleration);

    // Update velocities
    for (size_t i = 0; i < Natoms; i++) {
        for (size_t dim = 0; dim < 3; dim++) {
            velocity[i][dim] += 0.5 * (acceleration[i][dim] + new_acceleration[i][dim]) * dt;
            acceleration[i][dim] = new_acceleration[i][dim];
        }
    }

    free_2d(new_acceleration);
}

// Write the trajectory to an XYZ file
void write_xyz(FILE* file, size_t Natoms, Atom* atoms, double** coord, double LJ_potential, double kinetic_energy) {
    fprintf(file, "%zu\n", Natoms); // Number of atoms
    fprintf(file, "LJ=%.6f, KE=%.6f, Total=%.6f\n", LJ_potential, kinetic_energy, LJ_potential + kinetic_energy); // Comment line

    for (size_t i = 0; i < Natoms; i++) {
        fprintf(file, "%s %.6f %.6f %.6f\n", atoms[i].atom, coord[i][0] , coord[i][1] , coord[i][2] );
    }
}

// Write energies to a CSV file for plotting
void write_energies(FILE* energy_file, int step, double LJ_potential, double kinetic_energy) {
    double total_energy = LJ_potential + kinetic_energy;
    fprintf(energy_file, "%d,%.6f,%.6f,%.6f\n", step, LJ_potential, kinetic_energy, total_energy);
}

int main() {
    // File path
    char input_file[] = "inp.txt";
    char trajectory_file[] = "trajectory.xyz";
    char energy_file[] = "energies.csv";

    // Open input file
    FILE* file = fopen(input_file, "r");
    if (file == NULL) {
        perror("Error opening input file");
        return EXIT_FAILURE;
    }

    // Read atom data
    Atom atoms[MAX_ATOMS];
    int atom_count = 0;
    char line[256];

    // Skip header line
    fgets(line, sizeof(line), file);

    // Parse atoms and coordinates
    while (fgets(line, sizeof(line), file)) {
        char atom[3];
        double x, y, z, mass;
        if (sscanf(line, "%2s %lf %lf %lf %lf", atom, &x, &y, &z, &mass) == 5) {
            strcpy(atoms[atom_count].atom, atom);
            atoms[atom_count].x = x;
            atoms[atom_count].y = y;
            atoms[atom_count].z = z;
            atoms[atom_count].mass = mass;
            atom_count++;
        }
    }
    fclose(file);

    // Allocate memory for coordinates, distances, velocities, accelerations, and masses
    double** coord = malloc_2d(atom_count, 3);
    double** distance = malloc_2d(atom_count, atom_count);
    double** velocity = malloc_2d(atom_count, 3);
    double** acceleration = malloc_2d(atom_count, 3);
    double* mass = malloc(atom_count * sizeof(double));

    if (coord == NULL || distance == NULL || velocity == NULL || acceleration == NULL || mass == NULL) {
        printf("Memory allocation failed\n");
        free_2d(coord);
        free_2d(distance);
        free_2d(velocity);
        free_2d(acceleration);
        free(mass);
        return EXIT_FAILURE;
    }

    // Initialize coordinates, velocities, and masses
    for (int i = 0; i < atom_count; i++) {
        coord[i][0] = atoms[i].x;
        coord[i][1] = atoms[i].y;
        coord[i][2] = atoms[i].z;

        velocity[i][0] = 0.0;
        velocity[i][1] = 0.0;
        velocity[i][2] = 0.0;

        mass[i] = atoms[i].mass;
    }

    compute_distances(atom_count, coord, distance);
    compute_acc(atom_count, coord, mass, distance, acceleration);

    // Compute initial energies
    double LJ_potential = compute_LJ_potential(0.0661, 0.3345, atom_count, distance);
    double kinetic_energy = compute_kinetic_energy(atom_count, velocity, mass);
    double total_energy = LJ_potential + kinetic_energy;

    // Print initial coordinates and energies
    printf("Initial Coordinates and Energies:\n");
    printf("Coordinates:\n");
    for (int i = 0; i < atom_count; i++) {
        printf("%s: (%.6f, %.6f, %.6f)\n", atoms[i].atom, coord[i][0], coord[i][1], coord[i][2]);
    }
    printf("Lennard-Jones Potential: %.6f\n", LJ_potential);
    printf("Kinetic Energy: %.6f\n", kinetic_energy);
    printf("Total Energy: %.6f\n\n", total_energy);

    // Open output files
    FILE* output = fopen(trajectory_file, "w");
    FILE* energy_output = fopen(energy_file, "w");
    if (output == NULL || energy_output == NULL) {
        perror("Error opening output files");
        return EXIT_FAILURE;
    }

    // Write CSV header for energies
    fprintf(energy_output, "Step,LJ_Potential,Kinetic_Energy,Total_Energy\n");

    // Simulation loop
    for (int step = 0; step < TOTAL_STEPS; step++) {
        if (step % OUTPUT_INTERVAL == 0) {
            write_xyz(output, atom_count, atoms, coord, LJ_potential, kinetic_energy);
        }

        write_energies(energy_output, step, LJ_potential, kinetic_energy);
        verlet_update(atom_count, TIMESTEP, coord, velocity, acceleration, distance, mass);

        // Recompute energies after update
        LJ_potential = compute_LJ_potential(0.0661, 0.3345, atom_count, distance);
        kinetic_energy = compute_kinetic_energy(atom_count, velocity, mass);
    }

    fclose(output);
    fclose(energy_output);

    // Free memory
    free_2d(coord);
    free_2d(distance);
    free_2d(velocity);
    free_2d(acceleration);
    free(mass);

    return EXIT_SUCCESS;
}
