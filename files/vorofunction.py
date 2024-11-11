import numpy, sys, os
# from prettytable import PrettyTable
import gfort2py as gf
# import matplotlib.pyplot as plt
# from pdb import set_trace


def voronoi(x, sigma, Aeps, draw=False):

    An = 4
    Ax = numpy.array([Aeps, 1.0 - Aeps, 1.0 - Aeps, Aeps, Aeps])
    Ay = numpy.array([Aeps, Aeps, 1.0 - Aeps, 1.0 - Aeps, Aeps])
    Aflag = numpy.array([4, 2, 3, 1, 4])

    nsites = int(len(x)/2)
    sites = numpy.zeros((nsites,2))
    for i in range(nsites):
        sites[i,:] = x[2*i:2*i+2]

    vor = []    
    if nsites >= 3:
        SHARED_LIB_NAME=f'./libvoro.so'
        MOD_FILE_NAME='voro.mod'

        user_cache_dir = './gfortcache'
        if not os.path.exists(user_cache_dir):
        # Create the directory if necessary
            os.makedirs(user_cache_dir)

        x=gf.fFort(SHARED_LIB_NAME,MOD_FILE_NAME, cache_folder=user_cache_dir)
        
        nvmax = 500
        nv = int(0)
        vx = numpy.zeros(nvmax,dtype=float)
        vy = numpy.zeros(nvmax,dtype=float)
        vflag = numpy.zeros(nvmax,dtype=int)
        sstart = numpy.zeros(nsites+1,dtype=int)
        istop = int(0)
        mydict = x.voronoi(nsites,sites.T,An,Ax,Ay,Aflag,nvmax,sstart,nv,vx,vy,vflag,istop)[1]
        istop = (mydict)['istop']

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
    
    if nsites == 2:
        print('Warning: choose nsites >= 3')            
        sys.exit()
                    
    return vor, istop