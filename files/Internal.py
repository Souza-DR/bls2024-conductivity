from fenics import * #ARRUMAR
import numpy, time, random, sys, os
from prettytable import PrettyTable
import gfort2py as gf
import matplotlib.pyplot as plt
from pdb import set_trace
        
def entire_bd(x, on_boundary):
    return on_boundary        

def tag_internal_boundary():
  
  tdim = mesh.topology().dim()
  mesh.init(tdim-1, tdim)                        # Creates connectivities between facets and cells
  facet_to_cell = mesh.topology()(tdim-1, tdim)  # MeshConnectivity object (stores connections between facets and cells)
  domain_values = domains.array()                # For each cell, assumes 1 or 0 (depending on subdomain)
  facet_values = int_boundary.array()            # For each facet, it will assume 1 if is an internal boundary and 0 otherwise

  
  for facet in range(len(facet_values)):
      cells = facet_to_cell(facet)                 # Returns the index of the cells which are connected to this facet
    
      if len(cells) == 2:                           # If the facet is conected with two cells, then it is an internal facet
          values = domain_values[cells]
          if values[0] != values[1]:
            facet_values[facet] = numpy.max(values) + (nsites-1)*numpy.min(values)
            
def omega(x):
    for cell in cells(mesh):
        p1 = numpy.array([cell.get_vertex_coordinates()[0], cell.get_vertex_coordinates()[1]])
        p2 = numpy.array([cell.get_vertex_coordinates()[2], cell.get_vertex_coordinates()[3]])
        p3 = numpy.array([cell.get_vertex_coordinates()[4], cell.get_vertex_coordinates()[5]])
        ic = incenter(p1, p2, p3)

        kmin = 0
        distmin = numpy.linalg.norm(ic - x[0:2])
        for k in range(1, nsites):
            dist = numpy.linalg.norm(ic - x[2*k:2*k+2])
            if dist < distmin:
                kmin = k
                distmin = dist

        domains.array()[cell.index()] = kmin

def voronoi(x, draw):

    An = 4
    Ax = numpy.array([0.0, 1.0, 1.0, 0.0, 0.0])
    Ay = numpy.array([0.0, 0.0, 1.0, 1.0, 0.0])
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

def rotate(x):
    return numpy.array([-x[1], x[0]])

def incenter(A, B, C):
    a = numpy.linalg.norm(B-C)
    b = numpy.linalg.norm(C-A)
    c = numpy.linalg.norm(A-B)
    return numpy.array([(a*A[0]+b*B[0]+c*C[0])/(a+b+c), (a*A[1]+b*B[1]+c*C[1])/(a+b+c)])

def xinit(ninit):
    for j in range(ninit):
        # Perturbando a solução e projetando em D = [0,1]x[0,1]
        for i in range(2*nsites):
            xini[i] = max(0.0, min(1.0, solx[i] + 0.1*(2*random.random() - 1)))
    return xini

def noise():
    solx = get_data(nsites)
    noise_num = 0.0 
    noise_denom = 0.0
    #  Tag subdomains
    omega(solx)
    dx = Measure("dx", domain=mesh, subdomain_data=domains) 

    for alpha in range(nsources):
        zeta_alpha = 1.0 / weightsG[alpha]
        hsolve_clean = Function(V1)
        hsolve = Function(V1)
        htrial = TrialFunction(V1)
        vtest = TestFunction(V1)
        a = 0.0
        for i in range(nsites):
            a = a + (Constant(sigma[i]) * inner(grad(htrial), grad(vtest))) * dx(i)
        L = fsource[alpha] * vtest * dx        

        bcdir = DirichletBC(V1, Constant(0), entire_bd)
        # Solve
        solve(a == L, hsolve_clean, bcdir)

        # add noise to the data
        max_hsolve = numpy.abs(hsolve_clean.vector()[:]).max()
        h_perturb = numpy.random.default_rng(seed=123456).normal(loc= 0, scale=noise_coeff*max_hsolve, size=hsolve.vector()[:].size)
        hsolve.vector()[:] = hsolve_clean.vector()[:] + h_perturb
        noise_num = noise_num + zeta_alpha*assemble( (hsolve - hsolve_clean)**2 * dx ) 
        noise_denom = noise_denom +  zeta_alpha*assemble( hsolve_clean**2 * dx )
    
    noise_level = numpy.sqrt(noise_num / noise_denom)
    return noise_level

def evalfg(n, x, vor, greq=False):
    nsites = int(n/2)
    sites = numpy.zeros((nsites,2))
    for i in range(nsites):
        sites[i,:] = x[2*i:2*i+2]
    
    omega(x)
    dx = Measure("dx", domain=mesh, subdomain_data=domains)

    int_boundary.set_all(0)
    tag_internal_boundary()
    dS_int = Measure("dS", domain=mesh, subdomain_data=int_boundary)

    bcdir = DirichletBC(V1, Constant(0), entire_bd)

    G = 0.0
    gradGbound = numpy.zeros(2*nsites)

    for alpha in range(nsources):
 
        # Define rhs and lhs of PDE
        a = 0.0
        for i in range(nsites):
            a = a + (Constant(sigma[i]) * inner(grad(utf), grad(wtest))) * dx(i)
        L = fsource[alpha] * wtest * dx

        solve(a == L, usol, bcdir)

        # Compute cost function        
        Galpha = assemble( (usol - h[alpha])**2 * dx )
        if weightsG[alpha] == -1.0:
            weightsG[alpha] = 1.0/Galpha
    
        G = G + 0.5 * weightsG[alpha] * Galpha
        # print('Galpha =', Galpha)
        
        if greq:
            # Define rhs and lhs of PDE (?)
            a = 0.0
            for i in range(nsites):
                a = a + (Constant(sigma[i]) * inner(grad(ptf), grad(wtest))) * dx(i)
            L = - (usol - h[alpha]) * wtest * dx

            
            solve(a == L, psol, bcdir)

            for i in range(nsites):
                ### Boundary Expression ###
                grad_usol = grad(usol)
                grad_psol = grad(psol)

                # the terms that are continuous can be removed because they cancel out at the interface of two cells
                
                # tempplus = 0.5*(usol - h[alpha])**2 - fsource*psol  + Constant(sigma[i])*inner(grad_usol('+'), grad_psol('+'))

                tempplus = Constant(sigma[i])*inner(grad_usol('+'), grad_psol('+'))
                S1plus = Identity(2)*tempplus 
                S1plus = S1plus - outer(grad_psol('+'),Constant(sigma[i])*grad_usol('+')) - outer(Constant(sigma[i])*grad_usol('+'), grad_psol('+'))

                # the terms that are continuous can be removed because they cancel out at the interface of two cells

                #tempminus = 0.5*(usol - h[alpha])**2 - fsource*psol  + Constant(sigma[i])*inner(grad_usol('-'), grad_psol('-')) 

                tempminus = Constant(sigma[i])*inner(grad_usol('-'), grad_psol('-'))
                S1minus = Identity(2)*tempminus                
                S1minus = S1minus - outer(grad_psol('-'),Constant(sigma[i])*grad_usol('-')) - outer(Constant(sigma[i])*grad_usol('-'), grad_psol('-'))  

                for r in range(len(vor[i])):
                    v = ((vor[i])[r])[0]
                    vflag = ((vor[i])[r])[1]
                    vnext = ((vor[i])[(r+1)%len(vor[i])])[0]
                    if (vflag >= 0):
                        # As células vizinhas da aresta interna E são a_i e a_vflag
                        nu = rotate(numpy.array(v) - numpy.array(vnext))/numpy.linalg.norm(numpy.array(v) - numpy.array(vnext))
                        den = numpy.linalg.norm(sites[vflag,:] - sites[i, :])
                        xi = int(numpy.max([vflag, i]) + (nsites-1)*numpy.min([vflag, i]))

                        if i > vflag:
                            S1nunu = S1plus[0,0] * nu[0] * nu[0] + S1plus[1,1] * nu[1] * nu[1] + (S1plus[0,1] + S1plus[1,0]) * nu[0] * nu[1]
                            hkE0 = assemble(S1nunu*Expression("x[0] - ak", ak = sites[vflag,0], degree = 2)*dS_int(xi)+Constant(0)*dx)/den
                            hkE1 = assemble(S1nunu*Expression("x[1] - ak", ak = sites[vflag,1], degree = 2)*dS_int(xi)+Constant(0)*dx)/den
                            hiE0 = assemble(S1nunu*Expression("x[0] - ai", ai = sites[i, 0], degree = 2)*dS_int(xi)+Constant(0)*dx)/den
                            hiE1 = assemble(S1nunu*Expression("x[1] - ai", ai = sites[i, 1], degree = 2)*dS_int(xi)+Constant(0)*dx)/den
                        else:
                            S1nunu = S1minus[0,0] * nu[0] * nu[0] + S1minus[1,1] * nu[1] * nu[1] + (S1minus[0,1] + S1minus[1,0]) * nu[0] * nu[1]
                            hkE0 = assemble(S1nunu*Expression("x[0] - ak", ak = sites[vflag,0], degree = 2)*dS_int(xi)+Constant(0)*dx)/den
                            hkE1 = assemble(S1nunu*Expression("x[1] - ak", ak = sites[vflag,1], degree = 2)*dS_int(xi)+Constant(0)*dx)/den
                            hiE0 = assemble(S1nunu*Expression("x[0] - ai", ai = sites[i, 0], degree = 2)*dS_int(xi)+Constant(0)*dx)/den
                            hiE1 = assemble(S1nunu*Expression("x[1] - ai", ai = sites[i, 1], degree = 2)*dS_int(xi)+Constant(0)*dx)/den
                            
                        gradGbound[2*vflag:2*vflag+2] = gradGbound[2*vflag:2*vflag+2] - weightsG[alpha] * numpy.array([hkE0, hkE1])

                        gradGbound[2*i:2*i+2] = gradGbound[2*i:2*i+2] + weightsG[alpha] * numpy.array([hiE0, hiE1]) 
                                                        
    if greq:
        return G, gradGbound
    else:
        return G

def newfvalues(f, fvalues):
        M = len(fvalues)
        if max(fvalues) == 0.0:
            fvalues[0] = f
        else:
            oldfvalues = numpy.copy(fvalues[0:M-1])
            fvalues[0] = f
            fvalues[1:M] = oldfvalues
        return fvalues

def projectintoA(p):
    p[0] = max(0.01, min(p[0], 0.99))
    p[1] = max(0.01, min(p[1], 0.99))

def project(n, x):
    for i in range(int(n/2)):
        projectintoA(x[2*i:2*i+2])

def projectedgradient(n, x, eps, maxit, maxtime):
    starttime = time.time()
    project(n, x)

    vor, vorflag = voronoi(x, draw = True)
    if vorflag != 0:
        print('In projectedgradient, voronoi diagram is not well defined at (projection of) the initial guess')
        sys.exit()

    global weightsG
    weightsG = -numpy.ones(nsources)

    f, g = evalfg(n, x, vor, greq = True)

    M = 5
    fvalues = numpy.zeros(M)
    fvalues = newfvalues(f, fvalues)

    numevalf = 1

    gp = x - g
    project(n, gp)
    gp = gp - x
    normgp = numpy.linalg.norm(gp)
    
    lamspg = 1.0
    d = x - lamspg*g  
    project(n, d)
    d = d - x

    # Saving the values of G(x^0) and |grad G(x^0)|
    finit, normgpinit = f, normgp

    iter = 0

    myTable = PrettyTable(["iter", "fcnt", "G", "||gP||", 'normxdiff', "x"])
    myTable.add_row([iter, numevalf, f, normgp, 0.0, x])
    data = myTable.get_string()
    with open('./saida.txt', 'w') as txt:
        txt.write(data)
    print(myTable)
    smallstep = False
    smallxdiff = False
    TIME = False
    while normgp > eps and iter < maxit and not smallstep and not smallxdiff and not TIME:
        iter = iter + 1

        alpha = 1.0
        xtrial = x + alpha * d
        vor, vorflag = voronoi(xtrial, draw = False)
        if vorflag != 0:
            ftrial = float('inf')
        else:
            ftrial = evalfg(n, xtrial, vor)
            numevalf = numevalf + 1

        print('-->', alpha, ftrial, numevalf)
        
        gtd = numpy.inner(g, d)

        fmax = max(fvalues)

        while not (ftrial <= fmax + 1E-4 * alpha * gtd) and not smallstep:
            atrial = (- gtd * alpha**2) / (2.0*(ftrial - f - alpha * gtd))
            if not (0.1*alpha <= atrial and atrial <= 0.9*alpha):
                atrial = alpha / 2
            alpha = atrial
            xtrial = x + alpha * d
            vor, vorflag = voronoi(xtrial, draw = False)
            if vorflag != 0:
                ftrial = float('inf')
            else:
                ftrial = evalfg(n, xtrial, vor)
                numevalf = numevalf + 1

            print('-->', alpha, ftrial, numevalf)
        
            if alpha < eps:
                smallstep = True

        xdiff = xtrial - x

        x = xtrial
        vor, vorflag = voronoi(x, draw = True)
        if vorflag != 0:
            print('In projectedgradient, vorflag must be zero here and it is not.')

        gdiff = g

        f, g = evalfg(n, x, vor, greq = True)

        numevalf = numevalf + 1

        fvalues = newfvalues(f, fvalues)

        gdiff = g - gdiff

        gp = x - g
        project(n, gp)
        gp = gp - x
        normgp = numpy.linalg.norm(gp)

        lamspg = max(1.0E-3, min(numpy.inner(xdiff, xdiff) / numpy.inner(xdiff, gdiff), 1.0E+3))
        d = x - lamspg*g
        project(n, d)
        d = d - x

        myTable.add_row([iter, numevalf, f, normgp, numpy.linalg.norm(xdiff),  x])
        data = myTable.get_string()
        with open('./saida.txt', 'w') as txt:
            txt.write(data)
        print(myTable)

        if numpy.linalg.norm(xdiff) < eps:
            smallxdiff = True

        finaltime = time.time()
        CPU_time = finaltime - starttime

        if CPU_time >= maxtime*nsources:
            TIME = True
            

    if iter > maxit-1:
        print('Maximum number of iterations reached')
        flagsol = 0
    elif normgp <= eps:
        print('Small search direction')
        flagsol = 1
    elif smallstep:
        print('Too small step in line search')
        flagsol = 2
    elif smallxdiff:
        print("Lack of progress in the movement of sites.")
        flagsol = 3
    elif TIME:
        print('Timeout exceeded')
        flagsol = 4
    else:
        print('In projectedgradient, main loop ended by an unknown criterion.')
        flagsol = -1
    
    return flagsol, x, finit, normgpinit, f, normgp, iter, numevalf
    
    
def randomInit(ntrials, nsites):
    Gtrial = numpy.zeros(ntrials)
    xtrial_mat = numpy.random.rand(2*nsites, ntrials)  
    
    for k in range(ntrials):
        xtrial = xtrial_mat[:,k]
        vor, istop = voronoi(xtrial, draw = False)
        if istop != 0:
            print('The Voronoi method encountered an error when constructing a diagram from a random initialization.')
            sys.exit()
        Gtrial[k] = evalfg(2*nsites, xtrial, vor)
    
    print('---------------------')        
    print('Random initialization')    
    print('Gtrial = ', Gtrial)
    print('---------------------')        

    kmin = numpy.argmin(Gtrial)        
    xiniRand = xtrial_mat[:,kmin]
    
    return xiniRand


def vorDiag(sigma, m, mesh, V, sites):
    # Create a MeshFunction to store values on cells
    cell_values = MeshFunction('double', mesh, dim=2)

    domains = MeshFunction("size_t", mesh, mesh.topology().dim(), 0)
    domain_values = domains.array() 

    # Assign values to each cell
    for cell in cells(mesh):

        vertex_coordinates = cell.get_vertex_coordinates()
            
        midpoint_x = (vertex_coordinates[0] + vertex_coordinates[2] + vertex_coordinates[4])/3.0
        midpoint_y = (vertex_coordinates[1] + vertex_coordinates[3] + vertex_coordinates[5])/3.0
       
        dist = numpy.zeros(m) 
        for k in range(m):
            dist[k] = (midpoint_x - sites[2*k])**2 + (midpoint_y - sites[2*k+1])**2  
    
        cell_values[cell] = sigma[numpy.argmin(dist)]
        domain_values[cell.index()] = numpy.argmin(dist) # mark the subdomains with the index of the Voronoi cell
    
    f = Function(V)

    dm = V.dofmap()
    for cell in cells(mesh):
        f.vector()[dm.cell_dofs(cell.index())] = cell_values[cell] 
        
    return f


def evalerror(sigma, nsites, mesh, V, solx, xini, xfinal):

    ground = vorDiag(sigma, nsites, mesh, V, solx)  
  
    fini = vorDiag(sigma, nsites, mesh, V, xini)  

    ffinal = vorDiag(sigma, nsites, mesh, V, xfinal)  

    groundnorm = assemble( ground * dx )
    error_init = assemble( abs(fini - ground) * dx )/groundnorm
    error_opt = assemble( abs(ffinal - ground) * dx )/groundnorm

    print('error_init = ', error_init)
    print('error_opt = ', error_opt)

    return error_opt, error_init 

def get_data(nsites):
    #----------------
    # Ground truth data
    #----------------
    if nsites == 2:
        print('Warning: choose nsites >= 3')            
        sys.exit()
    if nsites == 3:
        solx = numpy.array([0.23, 0.65,
                            0.87, 0.78,
                            0.72, 0.32])
    if nsites == 4:
        solx = numpy.array([0.17, 0.89,
                            0.25, 0.20,
                            0.64, 0.41,
                            0.73, 0.93])
    if nsites == 5:
        solx = numpy.array([0.23, 0.91,
                            0.24, 0.33,
                            0.93, 0.38,
                            0.90, 0.85,
                            0.45, 0.54])
    if nsites == 6:
        solx = numpy.array([0.56, 0.67,
                            0.42, 0.37,
                            0.91, 0.10,
                            0.15, 0.50,
                            0.90, 0.86,
                            0.18, 0.90])
    if nsites == 7:
        solx = numpy.array([0.15, 0.94,
                            0.21, 0.47,
                            0.83, 0.39,
                            0.14, 0.10,
                            0.89, 0.78,
                            0.49, 0.76,
                            0.52, 0.30])
    if nsites == 8:
        solx = numpy.array([0.41, 0.29,
                            0.16, 0.54,
                            0.80, 0.48,
                            0.14, 0.10,
                            0.90, 0.85,
                            0.51, 0.77,
                            0.91, 0.21,
                            0.13, 0.96])
    if nsites == 9:
        solx = numpy.array([0.50, 0.21,
                            0.60, 0.25,
                            0.37, 0.11,
                            0.37, 0.87,
                            0.88, 0.69,
                            0.10, 0.85,
                            0.79, 0.71,
                            0.24, 0.17,
                            0.43, 0.80])
    if nsites == 10:
        solx = numpy.array([0.38, 0.78,
                            0.55, 0.33,
                            0.31, 0.84,
                            0.50, 0.49,
                            0.16, 0.36,
                            0.44, 0.15,
                            0.28, 0.69,
                            0.57, 0.13,
                            0.45, 0.44,
                            0.51, 0.89])
    if nsites == 11:    
        solx = numpy.array([0.15, 0.30,
                            0.51, 0.31,
                            0.71, 0.80,
                            0.32, 0.81,
                            0.70, 0.17,
                            0.89, 0.73,
                            0.52, 0.57,
                            0.74, 0.48,
                            0.76, 0.12,
                            0.61, 0.36,
                            0.20, 0.75])
    if nsites == 12:
        solx = numpy.array([0.89, 0.36,
                            0.70, 0.70,
                            0.31, 0.82,
                            0.16, 0.63,
                            0.50, 0.37,
                            0.71, 0.89,
                            0.42, 0.80,
                            0.87, 0.19,
                            0.25, 0.16,
                            0.60, 0.21,
                            0.27, 0.89,
                            0.35, 0.14])
    if nsites == 13:
        solx = numpy.array([0.71, 0.80, 
                            0.32, 0.81,
                            0.70, 0.17,
                            0.74, 0.69,
                            0.52, 0.57,
                            0.74, 0.48,
                            0.76, 0.12,
                            0.61, 0.36,
                            0.18, 0.73,
                            0.75, 0.37,
                            0.19, 0.34,
                            0.18, 0.19,
                            0.38, 0.86])
    if nsites == 14:
        solx = numpy.array([0.20, 0.39,
                            0.12, 0.35,
                            0.10, 0.56,
                            0.42, 0.71,
                            0.65, 0.13,
                            0.29, 0.77,
                            0.21, 0.15,
                            0.31, 0.55,
                            0.82, 0.45,
                            0.58, 0.36,
                            0.32, 0.88,
                            0.09, 0.73,
                            0.75, 0.83,
                            0.54, 0.58])
    if nsites == 15:
        solx = numpy.array([0.41, 0.29,
                            0.16, 0.54,
                            0.80, 0.48,
                            0.14, 0.10,
                            0.90, 0.85,
                            0.51, 0.80,
                            0.91, 0.21,
                            0.15, 0.90,
                            0.53, 0.06,
                            0.50, 0.50,
                            0.20, 0.75,
                            0.68, 0.97,
                            0.13, 0.34,
                            0.67, 0.27,
                            0.70, 0.07])

    if nsites == 16:
        solx = numpy.array([0.90, 0.90,
                            0.15, 0.90,
                            0.20, 0.29,
                            0.10, 0.79,
                            0.58, 0.11,
                            0.22, 0.77,
                            0.32, 0.64,
                            0.55, 0.61,
                            0.43, 0.27,
                            0.11, 0.56,
                            0.75, 0.40,
                            0.90, 0.11,
                            0.37, 0.26,
                            0.46, 0.75,
                            0.33, 0.45,
                            0.65, 0.65])

    if nsites == 17:
        solx = numpy.array([0.61, 0.85,
                            0.50, 0.10,
                            0.88, 0.50,
                            0.22, 0.86,
                            0.82, 0.62,
                            0.57, 0.50,
                            0.76, 0.53,
                            0.10, 0.35,
                            0.15, 0.16,
                            0.77, 0.40,
                            0.85, 0.20,
                            0.38, 0.45,
                            0.44, 0.84,
                            0.85, 0.80,
                            0.39, 0.70,
                            0.28, 0.59,
                            0.63, 0.73])

    if nsites == 18:
        solx = numpy.array([0.59, 0.28,
                            0.70, 0.82,
                            0.19, 0.60,
                            0.60, 0.53,
                            0.67, 0.32,
                            0.39, 0.66,
                            0.16, 0.37,
                            0.54, 0.65,
                            0.17, 0.86,
                            0.38, 0.20,
                            0.44, 0.86,
                            0.53, 0.56,
                            0.75, 0.36,
                            0.77, 0.59,
                            0.23, 0.53,
                            0.18, 0.56,
                            0.39, 0.87,
                            0.88, 0.15])

    if nsites == 19:
        solx = numpy.array([0.03, 0.42,
                            0.37, 0.03,
                            0.32, 0.90,
                            0.22, 0.58,
                            0.88, 0.14,
                            0.65, 0.83,
                            0.90, 0.54,
                            0.50, 0.70,
                            0.38, 0.23,
                            0.10, 0.75,
                            0.40, 0.50,
                            0.79, 0.96,
                            0.77, 0.56,
                            0.85, 0.23,
                            0.25, 0.45,
                            0.86, 0.46,
                            0.55, 0.23,
                            0.66, 0.46,
                            0.90, 0.74])

    return solx


    
################## 
################## 
# MAIN ALGORITHM
##################
##################

# Start Time
starttime = time.time()

nsites = int(sys.argv[1])
ninit = int(sys.argv[2])
nsources = int(sys.argv[3])
nmesh = int(sys.argv[4])
noise_coeff = float(sys.argv[5])
typeProblem = int(sys.argv[6])
typeinit = int(sys.argv[7])
Num = int(sys.argv[8])
Aeps = float(sys.argv[9])
maxit = int(sys.argv[10])
eps = float(sys.argv[11])
maxtime = int(sys.argv[12])

equidist = False
binary = False
ternary = False

if typeProblem == 1:
    equidist = True
elif typeProblem == 2:
    binary = True
elif typeProblem == 3:
    ternary = True
else:
    print("Invalid value for the problem type (choose 1, 2, or 3).")
    sys.exit()

# Create unit square mesh
mesh = UnitSquareMesh(nmesh,nmesh, 'crossed')
P1 = FiniteElement("Lagrange", mesh.ufl_cell(), 1)
V1 = FunctionSpace(mesh, P1)

# Define domains
domains = MeshFunction("size_t", mesh, mesh.topology().dim(), 0)

# Define internal interface domain
int_boundary = MeshFunction("size_t", mesh, mesh.topology().dim()-1, 0)

#----------------
# Data
#----------------
if equidist:
    sigma = numpy.arange(1.0, nsites+1)
if binary:
    sigma = [5.0]*numpy.ones(nsites) 
    sigma[-int(nsites/2):] = [10.0]
if ternary:
    sigma = numpy.zeros(nsites) 
    ternary_values = numpy.array([1.0, 5.0, 10.0])
    for l in range(nsites):
        for i in range(3):
            if l % 3 == i:
                sigma[l] = ternary_values[i]


fsource_list = [Expression('1.0', degree = 2), Expression("cos(pi*x[0])*cos(pi*x[1])", degree=2), Expression("sin(pi*x[0])*sin(pi*x[1])", degree=2), Expression("cos(2*pi*x[0])*cos(2*pi*x[1])", degree=2)]

if nsources == 1:
    fsource = fsource_list[0:1]
if nsources == 2:
    fsource = fsource_list[0:2]
if nsources == 3:
    fsource = fsource_list[0:3]
if nsources == 4:
    fsource = fsource_list[0:4]
if nsources >= 5:
    print('The maximum number of sources is 4')
    sys.exit()


solx = get_data(nsites)
   
#  Tag subdomains
omega(solx)

dx = Measure("dx", domain=mesh, subdomain_data=domains)

h = []
for alpha in range(nsources):

    hsolve_clean = Function(V1)
    hsolve = Function(V1)
    htrial = TrialFunction(V1)
    vtest = TestFunction(V1)
    a = 0.0
    for i in range(nsites):
        a = a + (Constant(sigma[i]) * inner(grad(htrial), grad(vtest))) * dx(i)
    L = fsource[alpha] * vtest * dx        

    bcdir = DirichletBC(V1, Constant(0), entire_bd)

    # Solve
    solve(a == L, hsolve_clean, bcdir)

    max_hsolve = numpy.abs(hsolve_clean.vector()[:]).max()
    h_perturb = numpy.random.default_rng(seed=123456).normal(loc= 0, scale=noise_coeff*max_hsolve, size=hsolve.vector()[:].size)
    hsolve.vector()[:] = hsolve_clean.vector()[:] + h_perturb

    # save data
    h.append(hsolve)


# Define trial functions, test functions and solutions
utf = TrialFunction(V1)
ptf = TrialFunction(V1)
wtest = TestFunction(V1)
usol = Function(V1)
psol = Function(V1)

#Initial parameter for projected gradient
random.seed(123456)
xini = numpy.zeros(2*nsites)

if typeinit == 1:
    xini = xinit(ninit)
elif typeinit == 2:
    global weightsG
    weightsG = numpy.ones(nsources)
    xini = xinit(ninit)
    vor, istop = voronoi(xini, draw = False)
    if istop != 0:
        print('The Voronoi method encountered an error while constructing the manufactured solution')
        sys.exit()
    feval = evalfg(2*nsites, xini, vor)
    for j in range(Num - 1):
        xcurrent = xinit(ninit)
        vor, istop = voronoi(xcurrent, draw = False)
        if istop != 0:
            print('The Voronoi method encountered an error while constructing the manufactured solution')
            sys.exit()
        fevaltrial = evalfg(2*nsites, xcurrent, vor)
        if fevaltrial < feval:
            xini = numpy.copy(xcurrent)
            feval = fevaltrial
elif typeinit == 3: 
    ntrials = Num
    xini = randomInit(ntrials, nsites)

flagsol, xfinal, finit, normgpinit, ffinal, normgpfinal, iter, numevalf = projectedgradient(2*nsites, xini, eps, maxit, maxtime)

## Compute noise level##
noise_level = noise()

# Final time
finaltime = time.time()
CPU_time = finaltime - starttime

#Computing the error
erroropt, errorinit = evalerror(sigma, nsites, mesh, V1, solx, xini, xfinal)

# Calculating the initial value of the cost function and the gradient, without weights
weightsG = numpy.ones(nsources)
vor, istop = voronoi(xini, draw = False)
if istop != 0:
    print('The Voronoi method encountered an error while constructing the manufactured solution')
    sys.exit()
finit1, ginit1 = evalfg(2*nsites, xini, vor, greq = True)
gp = xini - ginit1
project(2*nsites, gp)
gp = gp - xini
normgpinit1 = numpy.linalg.norm(gp)


vor, istop = voronoi(xfinal, draw = False)
if istop != 0:
    print('The Voronoi method encountered an error while constructing the manufactured solution')
    sys.exit()
ffinal1, gfinal1 = evalfg(2*nsites, xfinal, vor, greq = True)
gp = xfinal - gfinal1
project(2*nsites, gp)
gp = gp - xfinal
normgpfinal1 = numpy.linalg.norm(gp)

namefile = str(nsites)+'_'+str(nsources)+'_'+str(nmesh)+'_'+str(int(1000.0*noise_coeff))+'_'+str(ninit)

numpy.savez("./"+namefile+".npz", nsites=nsites, nsources=nsources, nmesh=nmesh, noise_coeff=noise_coeff, Aeps=Aeps, 
finit=finit, noise_level=noise_level, normgpinit=normgpinit, ffinal=ffinal, normgpfinal=normgpfinal, erroropt=erroropt, errorinit=errorinit, flagsol=flagsol, iter=iter, numevalf=numevalf, CPU_time=CPU_time, sigma=sigma, solx=solx, xini=xini, xfinal=xfinal, finit1=finit1, ffinal1=ffinal1, normgpinit1=normgpinit1, normgpfinal1=normgpfinal1)

    
print('End of main loop !!')