# There are two cases for the conductivity problem: with internal measurements and with boundary measurements.
# measurement = "Internal"
# or
measurement = "Boundary"

## Type of problemns: Type = 1 (Equidist); Type = 2 (binary); Type = 3 (ternary);##
typeproblem = 1

## Types of initialization
# typeinit = 1 Just one initialization and one optimization (For additional initializations and optimizations, add 2, 3, etc. to ninit_list) 
# Warning: if typeinit = 2 or 3, then make ninit_list=(1).
# typeinit = 2 N perturbations of the ground truth (The value of N must be provided) 
# typeinit = 3 N random initializations (The value of N must be provided) 
typeinit = 1
N = 100

# Choose the width of the domain boundary A: Aeps = 0.0 for the case of internal measurements, and Aeps > 0.0 for the case of boundary measurements.
Aeps = 0.05

# Spectral projected gradient method parameters
maxit = 500 # Maximum number of iterations
eps = 1E-6 # Epsilon that determines the lack of progress in the movement of sites, determines whether the step is too small in line search and whether the search direction is small.
maxtime = 14400 # Maximum execution time for each instance

nsites_list = [9]
ninit_list = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
# Warning: for the case with boundary measurements, the number of sources (nsources) must be 1 or 3. For the case with internal measurements, choose between 1 and 4.
nsources_list = [1, 3]
nmesh_list = [128]
noise_coeff_list = [0.0, 0.0025, 0.005, 0.01]