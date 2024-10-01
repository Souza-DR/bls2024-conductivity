import numpy as np
import sys, os
import shutil as sh

## Type of problemns: Type = 1 (Equidist); Type = 2 (binary); Type = 3 (ternary);##
typeproblem = 1

## Types of initialization
# typeinit = 1 Just one initialization and one optimization (For additional initializations and optimizations, add 2, 3, etc. to ninit_list) 
# Warning: if typeinit = 2 or 3, then make ninit_list=(1).
# typeinit = 2 N perturbations of the ground truth (The value of N must be provided) 
# typeinit = 3 N random initializations (The value of N must be provided) 
typeinit = 1
N = 100

#Fixed boundary layer width
Aeps = 0.05

# Spectral projected gradient method parameters
maxit = 500 # Maximum number of iterations
eps = 1E-6 # Epsilon that determines the lack of progress in the movement of sites, determines whether the step is too small in line search and whether the search direction is small.
maxtime = 10800 # Maximum execution time for each instance

nsites_list = [5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]
ninit_list = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
nsources_list = [1, 3]
nmesh_list = [16, 32]
noise_coeff_list = [0.005, 0.01, 0.02]

#Accessing each tests folder to removing potential old results
for nsites in nsites_list:
    for nmesh in nmesh_list:
        for noise_coeff in noise_coeff_list:
            for nsources in nsources_list:                
                for ninit in ninit_list:

                    input_files = "./"+str(nsites)+'_'+str(nsources)+'_'+str(nmesh)+'_'+str(int(1000.0*noise_coeff))+'_'+str(ninit)

                    # Validating the existence of instance data.
                    if os.path.exists(input_files):
                        #Delete the directories and all your files.
                        sh.rmtree(input_files)

#Moving the necessary files to run the tests with the adjusted parameters
for nsites in nsites_list:
    for nmesh in nmesh_list:
        for noise_coeff in noise_coeff_list:
            for nsources in nsources_list:                
                for ninit in ninit_list:
                    
                    input_files = "./"+str(nsites)+'_'+str(nsources)+'_'+str(nmesh)+'_'+str(int(1000.0*noise_coeff))+'_'+str(ninit)

                    if not os.path.exists(input_files):
                        os.makedirs(input_files)
                    
                    sh.copyfile("./files/drand.f90", input_files+"/drand.f90")
                    sh.copyfile("./files/geometry.f90", input_files+"/geometry.f90")
                    sh.copyfile("./files/geompack2.f90", input_files+"/geompack2.f90")
                    sh.copyfile("./files/mergesort.f90", input_files+"/mergesort.f90")
                    sh.copyfile("./files/vorintpols.f90", input_files+"/vorintpols.f90")
                    sh.copyfile("./files/voro.f90", input_files+"/voro.f90")
                    sh.copyfile("./files/voroextras.f90", input_files+"/voroextras.f90")
                    sh.copyfile("./files/EIT.py", input_files+"/EIT.py")

                    with open(input_files+"/instance.sh", "w") as f:
                        f.write(r"""#!/bin/bash
FILE=./results/data.npz
FILEONE=./results/data1.npz
DIR=./gfortcache

if [ -e  $FILE ] ; then
  rm ./results/data.npz
fi
                                
if [ -e  $FILEONE ] ; then
  rm ./results/data1.npz
fi

if [ -e  $DIR ] ; then
  rm -r ./gfortcache
fi
                                
gfortran -fPIC -shared -O4 vorintpols.f90
gfortran -fPIC -shared -O4 -c voro.f90
gfortran -fPIC -shared -O4 -o libvoro.so mergesort.f90 geometry.f90 geompack2.f90 vorintpols.f90 voro.f90

python_script="EIT.py"

typeproblem="""+str(typeproblem)+r"""

typeinit="""+str(typeinit)+r"""

N="""+str(N)+r"""

Aeps="""+str(Aeps)+r"""

maxit="""+str(maxit)+r"""

eps="""+str(eps)+r"""

maxtime="""+str(maxtime)+r"""

nsites_list=("""+str(nsites)+r""")
ninit_list=("""+str(ninit)+r""")
nsources_list=("""+str(nsources)+r""")
nmesh_list=("""+str(nmesh)+r""")
noise_coeff_list=("""+str(noise_coeff)+r""")

for nsites in "${nsites_list[@]}"; do
for noise_coeff in "${noise_coeff_list[@]}"; do
for ninit in "${ninit_list[@]}"; do
for nsources in "${nsources_list[@]}"; do
for nmesh in "${nmesh_list[@]}"; do
    echo "Running $python_script with argument: $arg1"
    python3 EIT.py $nsites $ninit $nsources $nmesh $noise_coeff $typeproblem $typeinit $N $Aeps $maxit $eps $maxtime
done
done
done
done
done
    """)
                        
for nsites in nsites_list:
    for nmesh in nmesh_list:
        for noise_coeff in noise_coeff_list:
            for nsources in nsources_list:                
                # Creating a script to run all tests simultaneously
                with open("./compile.sh", "w") as f:
                    f.write(r"""#!/bin/bash
                """)
                    for ninit in ninit_list:
                        input_files = "./"+str(nsites)+'_'+str(nsources)+'_'+str(nmesh)+'_'+str(int(1000.0*noise_coeff))+'_'+str(ninit)
                        f.write(r"""
(cd """+str(input_files)+r""" && bash instance.sh) & """)

                    f.write(r"""
wait
                            """)
                    
                #Running ten tests simultaneously
                os.system("bash compile.sh")


# Function `add` that merges two arrays into a single one
def add(old, current):
        total = len(old) + len(current)
        new = np.zeros(total)
        new[0:len(old)] = old
        new[len(old):total] = current
        return new

#Creating the folder where the results will be stored
if not os.path.exists("./results"):
    # Create the directory if necessary
    os.makedirs('./results')

#Removing potential old results
if os.path.exists("./results/data.npz"):
    os.system('rm ./results/data.npz')
if os.path.exists("./results/data1.npz"):
    os.system('rm ./results/data1.npz')

#Accessing each test's folder to retrieve their respective data
for nsites in nsites_list:
    for nmesh in nmesh_list:
        for noise_coeff in noise_coeff_list:
            for nsources in nsources_list:
                for ninit in ninit_list:
                                    
                    input_files = "./"+str(nsites)+'_'+str(nsources)+'_'+str(nmesh)+'_'+str(int(1000.0*noise_coeff))+'_'+str(ninit)

                    # Validating the existence of instance data to confirm the success or failure of the instance creation.
                    if os.path.exists(input_files+"/results/data.npz"):                  
                        # Checking if any data has already been stored or if we are in the first loop
                        if not os.path.exists("./results/data.npz"):
                            # Saving the data of the first instance
                            sh.copyfile(input_files+"/results/data.npz", "./results/data.npz")
                        else:
                            # Loading data from previous instances
                            loaded_data = np.load("./results/data.npz")

                            nsites_array = loaded_data["nsites_array"]
                            nsources_array = loaded_data["nsources_array"]
                            nmesh_array = loaded_data["nmesh_array"]
                            noise_coeff_array = loaded_data["noise_coeff_array"]
                            noise_level_array = loaded_data["noise_level_array"]
                            Aeps_array = loaded_data["Aeps_array"]
                            finit_array = loaded_data["finit_array"]
                            normgpinit_array = loaded_data["normgpinit_array"]
                            ffinal_array = loaded_data["ffinal_array"]
                            normgpfinal_array  = loaded_data["normgpfinal_array"]
                            erroropt_array  = loaded_data["erroropt_array"]
                            errorinit_array  = loaded_data["errorinit_array"]
                            flagsol_array = loaded_data["flagsol_array"]
                            iter_array = loaded_data["iter_array"]
                            numevalf_array = loaded_data["numevalf_array"]
                            CPU_time_array = loaded_data["CPU_time_array"]
                            sigma_array = loaded_data["sigma_array"]
                            solx_array = loaded_data["solx_array"]
                            xini_array = loaded_data["xini_array"]
                            xfinal_array = loaded_data["xfinal_array"]

                            #Loading data from the current instance
                            loaded_temp = np.load(input_files+"/results/data.npz")

                            nsites_current = loaded_temp["nsites_array"]
                            nsources_current = loaded_temp["nsources_array"]
                            nmesh_current = loaded_temp["nmesh_array"]
                            noise_coeff_current = loaded_temp["noise_coeff_array"]
                            Aeps_current = loaded_temp["Aeps_array"]
                            finit_current = loaded_temp["finit_array"]
                            noise_level_current = loaded_temp["noise_level_array"]
                            normgpinit_current = loaded_temp["normgpinit_array"]
                            ffinal_current = loaded_temp["ffinal_array"]
                            normgpfinal_current = loaded_temp["normgpfinal_array"]
                            erroropt_current = loaded_temp["erroropt_array"]
                            errorinit_current = loaded_temp["errorinit_array"]
                            flagsol_current = loaded_temp["flagsol_array"]
                            iter_current = loaded_temp["iter_array"]
                            numevalf_current = loaded_temp["numevalf_array"]
                            CPU_time_current = loaded_temp["CPU_time_array"]
                            sigma_current = loaded_temp["sigma_array"]
                            solx_current = loaded_temp["solx_array"]
                            xini_current = loaded_temp["xini_array"]
                            xfinal_current = loaded_temp["xfinal_array"]

                            # Appending the current array to the existing array
                            nsites_array = add(nsites_array, nsites_current)
                            nsources_array = add(nsources_array, nsources_current)
                            nmesh_array = add(nmesh_array, nmesh_current)
                            noise_coeff_array= add(noise_coeff_array, noise_coeff_current)
                            Aeps_array = add(Aeps_array, Aeps_current)
                            finit_array = add(finit_array, finit_current)
                            noise_level_array = add(noise_level_array, noise_level_current)
                            normgpinit_array = add(normgpinit_array, normgpinit_current)
                            ffinal_array  = add(ffinal_array, ffinal_current)
                            normgpfinal_array  = add(normgpfinal_array, normgpfinal_current)
                            erroropt_array = add(erroropt_array, erroropt_current)
                            errorinit_array = add(errorinit_array, errorinit_current)
                            flagsol_array = add(flagsol_array, flagsol_current) 
                            iter_array = add(iter_array, iter_current)
                            numevalf_array = add(numevalf_array, numevalf_current)
                            CPU_time_array = add(CPU_time_array, CPU_time_current)
                            sigma_array = add(sigma_array, sigma_current)
                            solx_array = add(solx_array, solx_current)
                            xini_array = add(xini_array, xini_current)
                            xfinal_array = add(xfinal_array, xfinal_current)

                            np.savez("./results/data.npz", nsites_array=nsites_array, nsources_array=nsources_array, nmesh_array=nmesh_array, noise_coeff_array=noise_coeff_array, Aeps_array=Aeps_array, finit_array=finit_array, noise_level_array=noise_level_array, normgpinit_array=normgpinit_array, ffinal_array=ffinal_array, normgpfinal_array=normgpfinal_array, erroropt_array=erroropt_array, errorinit_array=errorinit_array, flagsol_array=flagsol_array, iter_array=iter_array, numevalf_array=numevalf_array, CPU_time_array=CPU_time_array, sigma_array=sigma_array, solx_array=solx_array, xini_array=xini_array, xfinal_array=xfinal_array)

#Accessing each test's folder to retrieve their respective data
for nsites in nsites_list:
    for nmesh in nmesh_list:
        for noise_coeff in noise_coeff_list:
            for nsources in nsources_list:
                for ninit in ninit_list:
                    
                    input_files = "./"+str(nsites)+'_'+str(nsources)+'_'+str(nmesh)+'_'+str(int(1000.0*noise_coeff))+'_'+str(ninit)

                    # Validating the existence of instance data to confirm the success or failure of the instance creation.
                    if os.path.exists(input_files+"/results/data1.npz"):
                        # Checking if any data has already been stored or if we are in the first loop
                        if not os.path.exists("./results/data1.npz"):
                            # Saving the data of the first instance
                            sh.copyfile(input_files+"/results/data1.npz", "./results/data1.npz")
                            
                        else:
                            # Loading data from previous instances
                            loaded_data = np.load("./results/data1.npz")

                            
                            finit_array = loaded_data["finit_array"]
                            normgpinit_array = loaded_data["normgpinit_array"]
                            ffinal_array = loaded_data["ffinal_array"]
                            normgpfinal_array  = loaded_data["normgpfinal_array"]
                            
                            #Loading data from the current instance
                            loaded_temp = np.load(input_files+"/results/data1.npz")

                            finit_current = loaded_temp["finit_array"]
                            normgpinit_current = loaded_temp["normgpinit_array"]
                            ffinal_current = loaded_temp["ffinal_array"]
                            normgpfinal_current = loaded_temp["normgpfinal_array"]
                            
                            # Appending the current array to the existing array
                            finit_array = add(finit_array, finit_current)
                            normgpinit_array = add(normgpinit_array, normgpinit_current)
                            ffinal_array  = add(ffinal_array, ffinal_current)
                            normgpfinal_array  = add(normgpfinal_array, normgpfinal_current)

                            np.savez("./results/data1.npz", finit_array=finit_array, normgpinit_array=normgpinit_array, ffinal_array=ffinal_array, normgpfinal_array=normgpfinal_array)
                        
                        #Delete the directories and all your files.
                        sh.rmtree(input_files)

            


