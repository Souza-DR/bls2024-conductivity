import numpy as np
import gfort2py as gf
import os, sys
from pdb import set_trace
import shutil as sh
import subprocess

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

subprocess.run(["bash", "compile.sh"], cwd="./files/")

loaded_data = np.load("./results/data.npz")
loaded_data1 = np.load("./results/data1.npz")

nsites_array = loaded_data["nsites_array"]
nsources_array = loaded_data["nsources_array"]
nmesh_array = loaded_data["nmesh_array"]
noise_coeff_array = loaded_data["noise_coeff_array"]
sigma_array = loaded_data['sigma_array']
solx_array = loaded_data['solx_array']
xini_array = loaded_data['xini_array']
xfinal_array = loaded_data['xfinal_array']
Aeps_array = loaded_data['Aeps_array']

ffinal_array1 = loaded_data1["ffinal_array"]


ffinal1 = []
ffinal2 = []
ffinal3 = []
vencedor = []

for l in range(len(nsites_array)):
    if nsources_array[l] == 1:
        ffinal1.append(ffinal_array1[l])
    # if nsources_array[l] == 2:
    #     ffinal2.append(ffinal_array1[l])
    if nsources_array[l] == 3:
        ffinal3.append(ffinal_array1[l])

best_ffinal1 = np.min(np.array(ffinal1))
# best_ffinal2 = np.min(np.array(ffinal2))
best_ffinal3 = np.min(np.array(ffinal3))

# print(best_ffinal1)


for l in range(len(nsites_array)):
    if nsources_array[l] == 1 and ffinal_array1[l] <= best_ffinal1:
        vencedor.append(l)
    # if nsources_array[l] == 2 and ffinal_array1[l] <= best_ffinal2:
        # vencedor.append(l)
    if nsources_array[l] == 3 and ffinal_array1[l] <= best_ffinal3:
        vencedor.append(l)


position = np.zeros(len(nsites_array)+1)
sig_pst = np.zeros(len(nsites_array)+1)
for i in range(1, len(nsites_array)+1):
    for l in range(i):
        position[i] = position[i]+ nsites_array[l]
        sig_pst[i] = sig_pst[i] + nsites_array[l] + 1

position = np.array(position, dtype=int)
sig_pst = np.array(sig_pst, dtype=int)

alphabet = []
for i in range(97, 123):
    alphabet.append(chr(i))


l = 0
name = 'blsfig6'
for k in vencedor:
 
    nsites, nsources, nmesh, noise_coeff, Aeps, sigma, solx, xini, xfinal = nsites_array[k], nsources_array[k], nmesh_array[k], noise_coeff_array[k], Aeps_array[k], sigma_array[sig_pst[k]:sig_pst[k+1]], solx_array[2*position[k]:2*position[k+1]], xini_array[2*position[k]:2*position[k+1]], xfinal_array[2*position[k]:2*position[k+1]]

    

    if not os.path.exists('./results/figures'):
        # Create the directory if necessary
        os.makedirs('./results/figures')

    # if k % 4 == 0:
    #     l = 0
    #     name = 'blsfig'+str((k//4))

        
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

os.system('rm voronoi.mpx voronoi.log')