import numpy as np
import sys, os
import shutil as sh
import gfort2py as gf
import subprocess
from pdb import set_trace

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
Aeps = 0.0

# Spectral projected gradient method parameters
maxit = 500 # Maximum number of iterations
eps = 1E-6 # Epsilon that determines the lack of progress in the movement of sites, determines whether the step is too small in line search and whether the search direction is small.
maxtime = 10800 # Maximum execution time for each instance

nsites_list = [9]
ninit_list = [1]
# Warning: for the case with boundary measurements, the number of sources (nsources) must be 1 or 3. For the case with internal measurements, choose between 1 and 4.
nsources_list = [1]
nmesh_list = [128]
noise_coeff_list = [0.005]

#Moving the necessary files to run the tests with the adjusted parameters
for nsites in nsites_list:
    for nmesh in nmesh_list:
        for noise_coeff in noise_coeff_list:
            for nsources in nsources_list:                
                for ninit in ninit_list:
                    
                    input_files = "./instances/"+str(nsites)+'_'+str(nsources)+'_'+str(nmesh)+'_'+str(int(1000.0*noise_coeff))+'_'+str(ninit)

                    if not os.path.exists(input_files):
                        os.makedirs(input_files)
                    
                    sh.copyfile("./files/drand.f90", input_files+"/drand.f90")
                    sh.copyfile("./files/geometry.f90", input_files+"/geometry.f90")
                    sh.copyfile("./files/geompack2.f90", input_files+"/geompack2.f90")
                    sh.copyfile("./files/mergesort.f90", input_files+"/mergesort.f90")
                    sh.copyfile("./files/vorintpols.f90", input_files+"/vorintpols.f90")
                    sh.copyfile("./files/voro.f90", input_files+"/voro.f90")
                    sh.copyfile("./files/voroextras.f90", input_files+"/voroextras.f90")
                    sh.copyfile("./files/"+measurement+".py", input_files+"/"+measurement+".py")

                    with open(input_files+"/instance.sh", "w") as f:
                        f.write(r"""#!/bin/bash
DIR=./gfortcache

if [ -e  $DIR ] ; then
  rm -r ./gfortcache
fi
                                
gfortran -fPIC -shared -O4 vorintpols.f90
gfortran -fPIC -shared -O4 -c voro.f90
gfortran -fPIC -shared -O4 -o libvoro.so mergesort.f90 geometry.f90 geompack2.f90 vorintpols.f90 voro.f90

python_script=" """+measurement+r""".py"

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
    python3 """+measurement+r""".py $nsites $ninit $nsources $nmesh $noise_coeff $typeproblem $typeinit $N $Aeps $maxit $eps $maxtime
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
                        input_files = "./instances/"+str(nsites)+'_'+str(nsources)+'_'+str(nmesh)+'_'+str(int(1000.0*noise_coeff))+'_'+str(ninit)
                        f.write(r"""
(cd """+str(input_files)+r""" && bash instance.sh) & """)

                    f.write(r"""
wait
                            """)
                    
                #Running ten tests simultaneously
                os.system("bash compile.sh")

#Creating the folder where the results will be stored
if not os.path.exists("./results/data"):
    # Create the directory if necessary
    os.makedirs('./results/data')

#Accessing each test's folder to retrieve their respective data
for nsites in nsites_list:
    for nmesh in nmesh_list:
        for noise_coeff in noise_coeff_list:
            for nsources in nsources_list:
                for ninit in ninit_list:
                    # Validating the existence of instance data to confirm the success or failure of the instance creation
                    filename = str(nsites)+'_'+str(nsources)+'_'+str(nmesh)+'_'+str(int(1000.0*noise_coeff))+'_'+str(ninit)
                    if os.path.exists("./instances/"+filename):  
                        sh.move("./instances/"+filename+'/'+filename+'.npz', './results/data/'+filename+'.npz')
                   
#Delete the directorie and all your files.
sh.rmtree("./instances/")
            


