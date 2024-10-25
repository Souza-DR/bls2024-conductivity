import numpy as np
import sys, os
import shutil as sh
import gfort2py as gf
import subprocess
from pdb import set_trace
import parameters as prt

measurement = prt.measurement
typeproblem = prt.typeproblem
typeinit = prt.typeinit
N = prt.N
Aeps = prt.Aeps
maxit = prt.maxit
eps = prt.eps
maxtime = prt.maxtime
nsites_list = prt.nsites_list
ninit_list = prt.ninit_list
nsources_list = prt.nsources_list
nmesh_list = prt.nmesh_list
noise_coeff_list = prt.noise_coeff_list


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

#Creating a PDF with all tables and reconstructions from the tests
os.system("bash results.sh")