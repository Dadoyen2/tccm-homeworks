import psi4
import trexio

# Set memory and output options
psi4.set_memory('500 MB')
psi4.core.set_output_file('water_output.dat', False)

# Define the water molecule geometry
water = """
O
H 1 0.96
H 1 0.96 2 104.5
"""

# Set Psi4 options
psi4.set_options({
    'basis': 'sto-3g',
    'scf_type': 'pk',
    'reference': 'rhf',
    'e_convergence': 1e-8,
    'd_convergence': 1e-8
})

# Perform Hartree-Fock calculation
energy = psi4.energy('scf', molecule=water)

# Access wavefunction
wfn = psi4.core.Wavefunction.build(psi4.core.get_active_molecule(), psi4.core.get_options())

# Export to TREXIO
trexio_file = trexio.TrexioFile("water.trexio", mode="w")

trexio_file.write_nucleus_repulsion(wfn.nuclear_repulsion_energy())
trexio_file.write_electron_up_num(wfn.nalpha())
trexio_file.write_electron_down_num(wfn.nbeta())
trexio_file.write_mo_num(wfn.nmo())

# Write one-electron integrals
H_core = wfn.H().to_array()
trexio_file.write_mo_1e_int_core_hamiltonian(H_core.flatten())

# Write two-electron integrals (electron repulsion integrals)
eri = wfn.eri().to_array()
# Flatten the 4D ERI tensor into a list of non-zero integrals with their indices
# TREXIO expects integrals in the form <ij|kl>
# For simplicity, we'll write all integrals, but for larger systems, consider using symmetry
indices, values = trexio.convert_eri(eri)
trexio_file.write_mo_2e_int_eri(indices.flatten(), values)

# Write molecular orbital energies
mo_energies = wfn.epsilon_a().to_array()
trexio_file.write_mo_energy(mo_energies)

# Close TREXIO file
trexio_file.close()

