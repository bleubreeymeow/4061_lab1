LabA part3.1 Lennard Jones potential

For LJ Potential calculation, solid argon was being implemented. The values used for this calculation are stated below.

Reference: Introduction to solid state physics by Charles Kittel

SIGMA = 3.40 (Angstroms)
EPSILON = 0.0104233206 (eV)
distance cutoff =  2.5 * SIGMA (Angstroms)
lattice constant = 5.31 (Angstroms)
periodicty = (5,5,5)

LJ potential value from reference = -0.083 eV per atom
LJ potential from my code = -0.084 eV per atom


LabA part3.2 Coulomb buckingham potential

For Coulomb buckingham potential, NaCl was being implemented. The values used for this calculation are stated below. The subroutine of simple cubic lattice was being used, and hence there is a scaling of the periodicity and lattice constant by a factor of 2.

Reference: Investigation of H-center diffusion in sodium chloride by molecular dynamics simulation

lattice constant = 5.64 (Angstroms) (it was divided by a factor of two in the code because the lattice constant is the distance between each Na atom)
periodicty = (5,5,5)

Coulomb buckingham potential value from reference = -7.92 eV per ion pair
Coulomb buckingham potential from my code = -7.943796 eV per ion pair