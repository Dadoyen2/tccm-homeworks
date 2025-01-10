Installation and Compilation Instructions

This guide provides comprehensive instructions for installing, compiling, and running your Hartree-Fock (HF) and MP2 energy calculation program.

1. Prerequisites

Ensure you have the following installed:

GNU Compiler Collection (gcc)
pkg-config**
HDF5 (libhdf5-dev on Ubuntu)
Make (for building the project)

On Ubuntu/Debian-based systems, you can install these prerequisites with:

```bash
sudo apt update && sudo apt install build-essential pkg-config libhdf5-dev


2. Install TREXIO

wget [https://github.com/TREX-CoE/trexio/releases/download/v2.5.0/trexio-2.5.0.tar.gz](https://github.com/TREX-CoE/trexio/releases/download/v2.5.0/trexio-2.5.0.tar.gz)
tar -zxvf trexio-2.5.0.tar.gz
cd trexio-2.5.0
./configure
make
sudo make install
sudo ldconfig /usr/local/lib

3. Set Environment Variables 
Update LD_LIBRARY_PATH so the runtime linker can find TREXIO and HDF5:

echo 'export LD_LIBRARY_PATH=/usr/local/lib:/usr/lib/x86_64-linux-gnu/hdf5/serial:$LD_LIBRARY_PATH' >> ~/.bashrc
source ~/.bashrc


4. Compile the Program
Navigate to your project directory (containing your Makefile and src folder):


cd ~/tccm-homeworks/HW1
Run make to build the executable:
make

./compute_energy data/h2o.h5

./compute_energy data/c2h2.h5
./compute_energy data/ch4.h5
./compute_energy data/co2.h5
./compute_energy data/h3coh.h5
./compute_energy data/hcn.h5