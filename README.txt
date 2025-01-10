# Running the Simulation

This document provides a detailed guide on how to set up, compile, and execute the Lennard-Jones simulation program. Follow the instructions below to ensure smooth execution.

## **Prerequisites**
Before proceeding, ensure the following are available on your Linux system:
- A C compiler, preferably GCC.
- A suitable working environment to compile and run C programs.

## **Setup Instructions**

 **Create the Input File**
   - Create a text file named `input.txt` in the same directory as the program.
   - Format the file as follows:
     - The **first line** specifies the number of atoms.
     - Each subsequent line specifies the atom type, its Cartesian coordinates (x, y, z), and its mass.

     Example of `input.txt`:
     ```
     2
     Ar 0.0 0.0 0.0 39.95
     Ar 0.0 0.0 0.5 39.95
     ```

 **Compile the Program**
   - Ensure the source code file (e.g., `program.c`) is in the same directory as the `input.txt` file.
   - Compile the program using the GCC compiler with the math library linked:
     ```bash
     gcc program.c -o program -lm
     ```

**Run the Program**
   - Execute the compiled program:
     ```bash
     ./program
     ```

**Output Files**
   - After execution, two output files will be generated in the current directory:
     - `trajectory.xyz`: Contains the trajectory data of the simulation in XYZ format.
     - `energies.csv`: Contains the energy data (Lennard-Jones potential, kinetic energy, and total energy) for each step of the simulation.


To simplify compilation and execution, you can use a `Makefile`. Create a `Makefile` with the following content:

```makefile
all:
	gcc program.c -o program -lm

run:
	./program
```



## **Notes**
- Ensure the input file `input.txt` and the executable are in the same directory.
- The program assumes proper formatting of the input file.
- Both the `trajectory.xyz` and `energies.csv` files will be saved in the current working directory.

