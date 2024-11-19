import numpy as np
import sys, os, random
import shutil as sh
import parameters as prt
import files.groundtruth as gt

measurement = prt.measurement
typeproblem = prt.typeproblem
typeinit = prt.typeinit
Num = prt.Num
Aeps = prt.Aeps
maxit = prt.maxit
eps = prt.eps
maxtime = prt.maxtime
nsites_list = prt.nsites_list
type_GT = prt.type_Ground_Truth
ninit_list = prt.ninit_list
nsources_list = prt.nsources_list
nmesh_list = prt.nmesh_list
noise_coeff_list = prt.noise_coeff_list


if measurement == 'Boundary' and Aeps == 0.0:
    print('Warning: For the case of boundary measurements, the width of the domain boundary should be Aeps > 0.0')
    sys.exit()

if measurement == 'Internal' and Aeps > 0.0:
    print('Warning: For the case of internal measurements, the width of the domain boundary should be Aeps = 0.0')
    sys.exit()

if measurement == 'Boundary' and not set(nsources_list).issubset({1, 3}):
    print('Warning: For the case with boundary measurements, the number of sources (nsources) must be 1 or 3')
    sys.exit()

if measurement == 'Internal' and not set(nsources_list).issubset({1, 2, 3, 4}):
    print('Warning: For cases with internal measurements, the number of sources (nsources) must be chosen between 1 and 4')
    sys.exit()

if typeinit not in [1, 2, 3]:
    print('Error: The variable typeinit must be equal to 1, 2, or 3.')
    sys.exit()

if (typeinit == 2 or typeinit == 3) and max(ninit_list) > 1:
    print('Warning: if typeinit = 2 or 3, then make ninit_list=(1)')
    sys.exit()


#Moving the necessary files to run the tests with the adjusted parameters
for nsites in nsites_list:
    for nmesh in nmesh_list:
        for noise_coeff in noise_coeff_list:
            for nsources in nsources_list:                
                for ninit in ninit_list:
                    
                    namefile = str(nsites)+'_'+str(nsources)+'_'+str(nmesh)+'_'+str(int(1000.0*noise_coeff))+'_'+str(ninit)

                    if not os.path.exists("./instances/"+namefile):
                        os.makedirs("./instances/"+namefile)
                    
                    if type_GT == "default":                            
                        solx = gt.get_data(nsites)
                    elif type_GT == "input":
                        solx = np.zeros(2*nsites)
                        for j in range(2*nsites):
                            if j%2 == 0:
                                try: 
                                    x = float(input(f"Provide the coordinate {Aeps} < x < {1.0 - Aeps} for site a_{j//2} [Example: 0.5]: "))

                                    if Aeps < x < 1.0 - Aeps:
                                        solx[j] = x
                                    else: 
                                        print(f"The provided coordinate is not within the interval ({Aeps}, {1.0 - Aeps})")    
                                        sys.exit()
                                except ValueError:
                                         print("Invalid input! Please provide a valid numeric value.")
                                         sys.exit()
                            else:
                                try:
                                    y = float(input(f"Provide the coordinate {Aeps} < y < {1.0 - Aeps} for site a_{j//2} [Example: 0.5]: "))
                                    if Aeps < y < 1.0 - Aeps:
                                        solx[j] = y
                                    else: 
                                        print(f"The provided coordinate is not within the interval ({Aeps}, {1.0 - Aeps})")
                                        sys.exit()
                                except ValueError:
                                        print("Invalid input! Please provide a valid numeric value.")
                                        sys.exit()
                    else:
                        solx = np.zeros(2*nsites)
                        for j in range(2*nsites):
                            lbda = random.random()
                            solx[j] = lbda*(Aeps + 0.01) + (1.0 - lbda)*(0.99 - Aeps)

                    np.savez("./instances/"+namefile+'/'+namefile+".npz", solx=solx)
                    
                    
                    sh.copyfile("./files/drand.f90", "./instances/"+namefile+"/drand.f90")
                    sh.copyfile("./files/geometry.f90", "./instances/"+namefile+"/geometry.f90")
                    sh.copyfile("./files/geompack2.f90", "./instances/"+namefile+"/geompack2.f90")
                    sh.copyfile("./files/mergesort.f90", "./instances/"+namefile+"/mergesort.f90")
                    sh.copyfile("./files/vorintpols.f90", "./instances/"+namefile+"/vorintpols.f90")
                    sh.copyfile("./files/voro.f90", "./instances/"+namefile+"/voro.f90")
                    sh.copyfile("./files/voroextras.f90", "./instances/"+namefile+"/voroextras.f90")
                    sh.copyfile("./files/"+measurement+".py", "./instances/"+namefile+"/"+measurement+".py")
                    sh.copyfile("./files/vorofunction.py", "./instances/"+namefile+"/vorofunction.py")
                                        


                    with open("./instances/"+namefile+"/instance.sh", "w") as f:
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

N="""+str(Num)+r"""

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

alphabet = []
for i in range(97, 123):
    alphabet.append(chr(i))

count = 0
if not os.path.exists('./results/figures'):
    # Create the directory if necessary
    os.makedirs('./results/figures')

#Accessing each test's folder to retrieve their respective data
for nsites in nsites_list:
    for nmesh in nmesh_list:
        for noise_coeff in noise_coeff_list:
            for nsources in nsources_list:
                for ninit in ninit_list:
                    # Validating the existence of instance data to confirm the success or failure of the instance creation
                    filename = str(nsites)+'_'+str(nsources)+'_'+str(nmesh)+'_'+str(int(1000.0*noise_coeff))+'_'+str(ninit)
                    if os.path.exists("./instances/"+filename):
                        
                        if count % 4 == 0:
                            l = 0
                            name = 'blsfig'+str((count//4))

                        sh.move("./instances/"+filename+'/'+filename+'.npz', './results/data/'+filename+'.npz')

                        sh.move("./instances/"+filename+'/solx.mp', './results/figures/'+name+str(alphabet[l])+'.mp')
                        sh.move("./instances/"+filename+'/solx.eps', './results/figures/'+name+str(alphabet[l])+'.eps')
                        l += 1

                        sh.move("./instances/"+filename+'/xini.mp', './results/figures/'+name+str(alphabet[l])+'.mp')
                        sh.move("./instances/"+filename+'/xini.eps', './results/figures/'+name+str(alphabet[l])+'.eps')
                        l += 1

                        sh.move("./instances/"+filename+'/xfinal.mp', './results/figures/'+name+str(alphabet[l])+'.mp')
                        sh.move("./instances/"+filename+'/xfinal.eps', './results/figures/'+name+str(alphabet[l])+'.eps')
                        l += 1
                        count += 1

                   
#Delete the directorie and all your files.
sh.rmtree("./instances/")

#Creating a PDF with all tables and reconstructions from the tests
os.system("bash results.sh")