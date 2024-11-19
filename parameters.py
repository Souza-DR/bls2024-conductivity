############################################
## Choose the width of the domain boundary A
############################################
## Aeps = 0.0 for the case of internal measurements, and Aeps > 0.0 for the case of boundary measurements.
Aeps = 0.05

##############################
## Choose the Number of Sites
##############################
nsites_list = [9]

########################################################
## Choose the Ground Truth Type: default, random, input
########################################################
## default: the default Ground Truths (3-19 sites) were constructed to satisfy the hypotheses required for the gradient calculation formulas of the objective function.
## random: draws (2*nsites) numbers from the interval (Aeps, 1.0 - Aeps)
## input: allows the user to manually input the coordinates of the Ground Truth sites, as long as they are within the interval (Aeps, 1.0 - Aeps).
type_Ground_Truth = "default"

###############################################################################################################
## There are two cases for the conductivity problem: with boundary measurements and with internal measurements.
###############################################################################################################
# measurement = "Internal"
## or
measurement = "Boundary"

##################################################################################
## Type of problemns: Type = 1 (Equidist); Type = 2 (binary); Type = 3 (ternary)
##################################################################################
typeproblem = 1

#######################################################################################################################################
## Types of initialization
    ## typeinit = 1 The initialization is a single perturbation of the ground truth (GT). Only one optimization is performed. (For additional initializations, such as new perturbations of the GT, and corresponding optimizations, add 2, 3, etc., to the ninit_list, i.e. ninit_list = [1, 2, 3, ... , n])
    ## typeinit = 2: N perturbations of the ground truth (the value of Num must be provided), and the perturbation that provides the least value of the function is chosen as the initial point.
    ## typeinit = 3 N random initializations (The value of Num must be provided) 
#######################################################################################################################################
typeinit = 1
ninit_list = [1]
Num = 100

################################################
## Spectral projected gradient method parameters
################################################
maxit = 500 ## Maximum number of iterations
eps = 1E-6 ## Epsilon that determines the lack of progress in the movement of sites, determines whether the step is too small in line search and whether the search direction is small.
maxtime = 14400 ## Maximum execution time for each instance

##############################################
## Choose the Number of Sources (Measurements)
##############################################
# Warning: for the case with boundary measurements, the number of sources (nsources) must be 1 or 3. For the case with internal measurements, choose between 1 and 4.
nsources_list = [1]

###############################################
## Choose the Mesh Size to Discretize the Domain
###############################################
nmesh_list = [128]

#################################################
## Choose the Noise Level to be Added to the Data 
#################################################
## the chosen value should be between 0 and 1, and the actual noise percentage will be close to the chosen value)
noise_coeff_list = [0.0]