# File: cc/Makefile

# ============================
#        Compiler Settings
# ============================

# Compiler to use
CC = gcc

# Compiler flags
# -I./src: Include headers from the src directory
# -O2: Optimization level 2
# -Wall: Enable all warnings
CFLAGS = -I./src -O2 -Wall $(shell pkg-config --cflags hdf5)

# Linker flags
# -L/usr/local/lib: Directory where TREXIO libraries are installed
# $(shell pkg-config --libs hdf5): Linker flags for HDF5
# -ltrexio: Link against TREXIO library
LDFLAGS = -L/usr/local/lib $(shell pkg-config --libs hdf5) -ltrexio

# ============================
#         Source Files
# ============================

# List of source files
SRC = src/main.c src/hf_energy.c src/mp2_energy.c

# Object files derived from source files
OBJ = $(SRC:.c=.o)

# Name of the final executable
EXEC = compute_energy

# ============================
#            Rules
# ============================

# Default target to build the executable
all: $(EXEC)

# Rule to link object files into the executable
$(EXEC): $(OBJ)
	$(CC) $(OBJ) -o $(EXEC) $(LDFLAGS)

# Pattern rule to compile .c files into .o files
%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

# Clean target to remove compiled object files and executable
clean:
	rm -f src/*.o $(EXEC)

