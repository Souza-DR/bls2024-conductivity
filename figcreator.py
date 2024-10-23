import numpy as np
import gfort2py as gf
import os, sys
from pdb import set_trace
import shutil as sh
import subprocess

nsites_list = [9]
ninit_list = [1]
# Warning: for the case with boundary measurements, the number of sources (nsources) must be 1 or 3. For the case with internal measurements, choose between 1 and 4.
nsources_list = [1]
nmesh_list = [128]
noise_coeff_list = [0.005]

def voronoi(x, sigma, draw):

    An = 4
    Ax = np.array([Aeps, 1.0 - Aeps, 1.0 - Aeps, Aeps, Aeps])
    Ay = np.array([Aeps, Aeps, 1.0 - Aeps, 1.0 - Aeps, Aeps])
    Aflag = np.array([4, 2, 3, 1, 4])

    nsites = int(len(x)/2)
    sites = np.zeros((nsites,2))
    for i in range(nsites):
        sites[i,:] = x[2*i:2*i+2]

    SHARED_LIB_NAME=f'./files/libvoro.so'
    MOD_FILE_NAME='./files/voro.mod'
    user_cache_dir = './files/gfortcache'

    if not os.path.exists(user_cache_dir):
        # Create the directory if necessary
            os.makedirs(user_cache_dir)
    x=gf.fFort(SHARED_LIB_NAME,MOD_FILE_NAME, cache_folder=user_cache_dir)
    
    nvmax = 500
    nv = int(0)
    vx = np.zeros(nvmax,dtype=float)
    vy = np.zeros(nvmax,dtype=float)
    vflag = np.zeros(nvmax,dtype=int)
    sstart = np.zeros(nsites+1,dtype=int)
    istop = int(0)
    mydict = x.voronoi(nsites,sites.T,An,Ax,Ay,Aflag,nvmax,sstart,nv,vx,vy,vflag,istop)[1]
    istop = (mydict)['istop']

    vor = []
    if istop == 0:
        sstart = (mydict)['sstart']
        nv = (mydict)['nv']
        vx = (mydict)['vx']
        vy = (mydict)['vy']
        vectorvflag = (mydict)['vflag']

        #Ajustando os índices das células
        for i in range(len(vectorvflag)):
            if vectorvflag[i] > 0:
                vectorvflag[i] = vectorvflag[i] - 1

        for i in range(nsites):
            cell = []
            for k in range(sstart[i] - 1,sstart[i+1] - 1):
                cell.append([[vx[k],vy[k]],vectorvflag[k]])
            vor.append(cell)

        if draw:
            colors = sigma.astype(int)
            x.drawvor(nsites,sites.T,colors,An,Ax,Ay,sstart,nv,vx,vy)
        
    return vor, istop

subprocess.run(["bash", "libraryfortran.sh"], cwd="./files/")

alphabet = []
for i in range(97, 123):
    alphabet.append(chr(i))

k = 0
for nsites in nsites_list:
    for nmesh in nmesh_list:
        for noise_coeff in noise_coeff_list:
            for nsources in nsources_list:
                for ninit in ninit_list:
                    
                    input_files = str(nsites)+'_'+str(nsources)+'_'+str(nmesh)+'_'+str(int(1000.0*noise_coeff))+'_'+str(ninit)

                    loaded_data = np.load("./results/data/"+input_files+".npz")
                   
                    sigma = loaded_data['sigma']
                    solx = loaded_data['solx']
                    xini = loaded_data['xini']
                    xfinal = loaded_data['xfinal']
                    Aeps = loaded_data['Aeps']

                    if not os.path.exists('./results/figures'):
                    # Create the directory if necessary
                        os.makedirs('./results/figures')

                    if k % 4 == 0:
                        l = 0
                        name = 'blsfig'+str((k//4))

                        
                    voronoi(solx, sigma, True)
                    os.system('mpost voronoi.mp')
                    os.rename('./voronoi.mp', './results/figures/'+name+str(alphabet[l])+'.mp')
                    os.rename('./voronoi.mps','./results/figures/'+name+str(alphabet[l])+'.eps')
                    l += 1

                    voronoi(xini, sigma, True)
                    os.system('mpost voronoi.mp')
                    os.rename('./voronoi.mp', './results/figures/'+name+str(alphabet[l])+'.mp')
                    os.rename('./voronoi.mps','./results/figures/'+name+str(alphabet[l])+'.eps')
                    l += 1

                    voronoi(xfinal, sigma, True)
                    os.system('mpost voronoi.mp')
                    os.rename('./voronoi.mp', './results/figures/'+name+str(alphabet[l])+'.mp')
                    os.rename('./voronoi.mps','./results/figures/'+name+str(alphabet[l])+'.eps')
                    l += 1
                    k += 1

                os.system('rm voronoi.mpx voronoi.log')