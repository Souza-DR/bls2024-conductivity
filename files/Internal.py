from fenics import *
import numpy as np
import time, random, sys, os
from prettytable import PrettyTable
import vorofunction as vf
        
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
            facet_values[facet] = np.max(values) + (nsites-1)*np.min(values)
            
def omega(x):

    for cell in cells(mesh):
        p1 = np.array([cell.get_vertex_coordinates()[0], cell.get_vertex_coordinates()[1]])
        p2 = np.array([cell.get_vertex_coordinates()[2], cell.get_vertex_coordinates()[3]])
        p3 = np.array([cell.get_vertex_coordinates()[4], cell.get_vertex_coordinates()[5]])
        ic = incenter(p1, p2, p3)

        kmin = 0
        distmin = np.linalg.norm(ic - x[0:2])
        for k in range(1, nsites):
            dist = np.linalg.norm(ic - x[2*k:2*k+2])
            if dist < distmin:
                kmin = k
                distmin = dist

        domains.array()[cell.index()] = kmin

def rotate(x):
    return np.array([-x[1], x[0]])

def incenter(A, B, C):
    a = np.linalg.norm(B-C)
    b = np.linalg.norm(C-A)
    c = np.linalg.norm(A-B)
    return np.array([(a*A[0]+b*B[0]+c*C[0])/(a+b+c), (a*A[1]+b*B[1]+c*C[1])/(a+b+c)])

def xinit(ninit):
    for j in range(ninit):
        # Perturbing the solution and projecting onto D = [Aeps, 1 - Aeps]x[Aeps, 1 - Aeps]
        for i in range(2*nsites):
            xini[i] = max(Aeps, min(1.0 - Aeps, solx[i] + 0.1*(2*random.random() - 1)))
    return xini

def noise(solx):
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
        max_hsolve = np.abs(hsolve_clean.vector()[:]).max()
        h_perturb = np.random.default_rng(seed=123456).normal(loc= 0, scale=noise_coeff*max_hsolve, size=hsolve.vector()[:].size)
        hsolve.vector()[:] = hsolve_clean.vector()[:] + h_perturb
        noise_num = noise_num + zeta_alpha*assemble( (hsolve - hsolve_clean)**2 * dx ) 
        noise_denom = noise_denom +  zeta_alpha*assemble( hsolve_clean**2 * dx )
    
    noise_level = np.sqrt(noise_num / noise_denom)
    return noise_level

def randomInit(ntrials, nsites):
    Gtrial = np.zeros(ntrials)
    xtrial_mat = np.random.rand(2*nsites, ntrials)  
    
    for k in range(ntrials):
        xtrial = xtrial_mat[:,k]
        vor, istop = vf.voronoi(xtrial, sigma, Aeps, draw = False)
        if istop != 0:
            print('The Voronoi method encountered an error when constructing a diagram from a random initialization.')
            sys.exit()
        Gtrial[k] = evalfg(2*nsites, xtrial, vor)
    
    print('---------------------')        
    print('Random initialization')    
    print('Gtrial = ', Gtrial)
    print('---------------------')        

    kmin = np.argmin(Gtrial)        
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
       
        dist = np.zeros(m) 
        for k in range(m):
            dist[k] = (midpoint_x - sites[2*k])**2 + (midpoint_y - sites[2*k+1])**2  
    
        cell_values[cell] = sigma[np.argmin(dist)]
        domain_values[cell.index()] = np.argmin(dist) # mark the subdomains with the index of the Voronoi cell
    
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

def evalfg(n, x, vor, greq=False):
    nsites = int(n/2)
    sites = np.zeros((nsites,2))
    for i in range(nsites):
        sites[i,:] = x[2*i:2*i+2]
    
    omega(x)
    dx = Measure("dx", domain=mesh, subdomain_data=domains)

    int_boundary.set_all(0)
    tag_internal_boundary()
    dS_int = Measure("dS", domain=mesh, subdomain_data=int_boundary)

    bcdir = DirichletBC(V1, Constant(0), entire_bd)

    G = 0.0
    gradGbound = np.zeros(2*nsites)

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
        
        if greq:
            # Define rhs and lhs of PDE
            a = 0.0
            for i in range(nsites):
                a = a + (Constant(sigma[i]) * inner(grad(ptf), grad(wtest))) * dx(i)
            L = - (usol - h[alpha]) * wtest * dx

            
            solve(a == L, psol, bcdir)

            for i in range(nsites):

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
                        # The neighboring cells of the inner edge E are a_i and a_vflag
                        nu = rotate(np.array(v) - np.array(vnext))/np.linalg.norm(np.array(v) - np.array(vnext))
                        den = np.linalg.norm(sites[vflag,:] - sites[i, :])
                        xi = int(np.max([vflag, i]) + (nsites-1)*np.min([vflag, i]))

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
                            
                        gradGbound[2*vflag:2*vflag+2] = gradGbound[2*vflag:2*vflag+2] - weightsG[alpha] * np.array([hkE0, hkE1])

                        gradGbound[2*i:2*i+2] = gradGbound[2*i:2*i+2] + weightsG[alpha] * np.array([hiE0, hiE1]) 
                                                        
    if greq:
        return G, gradGbound
    else:
        return G

def newfvalues(f, fvalues):
        M = len(fvalues)
        if max(fvalues) == 0.0:
            fvalues[0] = f
        else:
            oldfvalues = np.copy(fvalues[0:M-1])
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

    vor, vorflag = vf.voronoi(x, sigma, Aeps, draw = True)
    if vorflag != 0:
        print('In projectedgradient, voronoi diagram is not well defined at (projection of) the initial guess')
        sys.exit()

    global weightsG
    weightsG = -np.ones(nsources)

    f, g = evalfg(n, x, vor, greq = True)

    M = 5
    fvalues = np.zeros(M)
    fvalues = newfvalues(f, fvalues)

    numevalf = 1

    gp = x - g
    project(n, gp)
    gp = gp - x
    normgp = np.linalg.norm(gp)
    
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
    with open('./output.txt', 'w') as txt:
        txt.write(data)
    print(myTable)
    smallstep = False
    smallxdiff = False
    TIME = False
    while normgp > eps and iter < maxit and not smallstep and not smallxdiff and not TIME:
        iter = iter + 1

        alpha = 1.0
        xtrial = x + alpha * d
        vor, vorflag = vf.voronoi(xtrial, sigma, Aeps, draw = False)
        if vorflag != 0:
            ftrial = float('inf')
        else:
            ftrial = evalfg(n, xtrial, vor)
            numevalf = numevalf + 1

        print('-->', alpha, ftrial, numevalf)
        
        gtd = np.inner(g, d)

        fmax = max(fvalues)

        while not (ftrial <= fmax + 1E-4 * alpha * gtd) and not smallstep:
            atrial = (- gtd * alpha**2) / (2.0*(ftrial - f - alpha * gtd))
            if not (0.1*alpha <= atrial and atrial <= 0.9*alpha):
                atrial = alpha / 2
            alpha = atrial
            xtrial = x + alpha * d
            vor, vorflag = vf.voronoi(xtrial, sigma, Aeps, draw = False)
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
        vor, vorflag = vf.voronoi(x, sigma, Aeps, draw = True)
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
        normgp = np.linalg.norm(gp)

        lamspg = max(1.0E-3, min(np.inner(xdiff, xdiff) / np.inner(xdiff, gdiff), 1.0E+3))
        d = x - lamspg*g
        project(n, d)
        d = d - x

        myTable.add_row([iter, numevalf, f, normgp, np.linalg.norm(xdiff),  x])
        data = myTable.get_string()
        with open('./output.txt', 'w') as txt:
            txt.write(data)
        print(myTable)

        if np.linalg.norm(xdiff) < eps:
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
typeproblem = int(sys.argv[6])
typeinit = int(sys.argv[7])
Num = int(sys.argv[8])
Aeps = float(sys.argv[9])
maxit = int(sys.argv[10])
eps = float(sys.argv[11])
maxtime = int(sys.argv[12])

equidist = False
binary = False
ternary = False

if typeproblem == 1:
    equidist = True
elif typeproblem == 2:
    binary = True
elif typeproblem == 3:
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
    sigma = np.arange(1.0, nsites+1)
if binary:
    sigma = [5.0]*np.ones(nsites) 
    sigma[-int(nsites/2):] = [10.0]
if ternary:
    sigma = np.zeros(nsites) 
    ternary_values = np.array([3.0, 6.0, 9.0])
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

#----------------
# Ground truth data
#---------------- 
namefile = str(nsites)+'_'+str(nsources)+'_'+str(nmesh)+'_'+str(int(1000.0*noise_coeff))+'_'+str(ninit)
solx = (np.load("./"+namefile+".npz"))['solx']
# solx = gt.get_data(nsites)

# Checking if it was possible to construct the Voronoi diagram from the provided Ground Truth.
vor, istop = vf.voronoi(solx, sigma, Aeps, draw = False)
if istop != 0:
    print('The Voronoi method encountered an error while constructing the manufactured solution')
    sys.exit()
   
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

    max_hsolve = np.abs(hsolve_clean.vector()[:]).max()
    h_perturb = np.random.default_rng(seed=123456).normal(loc= 0, scale=noise_coeff*max_hsolve, size=hsolve.vector()[:].size)
    hsolve.vector()[:] = hsolve_clean.vector()[:] + h_perturb

    # save data
    h.append(hsolve)


# Define trial functions, test functions and solutions
utf = TrialFunction(V1)
ptf = TrialFunction(V1)
wtest = TestFunction(V1)
usol = Function(V1)
psol = Function(V1)

# Defining the seed to be used in the random function for generating perturbations and random initializations.
random.seed(123456)
xini = np.zeros(2*nsites)

if typeinit == 1:
    xini = xinit(ninit)
elif typeinit == 2:
    global weightsG
    weightsG = np.ones(nsources)
    xini = xinit(ninit)
    vor, istop = vf.voronoi(xini, sigma, Aeps, draw = False)
    if istop != 0:
        print('The Voronoi method encountered an error while constructing the initial point')
        sys.exit()
    feval = evalfg(2*nsites, xini, vor)
    for j in range(Num - 1):
        xcurrent = xinit(ninit)
        vor, istop = vf.voronoi(xcurrent, sigma, Aeps, draw = False)
        if istop != 0:
            print('The Voronoi method encountered an error while constructing the initial point')
            sys.exit()
        fevaltrial = evalfg(2*nsites, xcurrent, vor)
        if fevaltrial < feval:
            xini = np.copy(xcurrent)
            feval = fevaltrial
elif typeinit == 3: 
    ntrials = Num
    xini = randomInit(ntrials, nsites)

flagsol, xfinal, finit, normgpinit, ffinal, normgpfinal, iter, numevalf = projectedgradient(2*nsites, xini, eps, maxit, maxtime)

## Compute noise level##
noise_level = noise(solx)

# Final time
finaltime = time.time()
CPU_time = finaltime - starttime

#Computing the error
erroropt, errorinit = evalerror(sigma, nsites, mesh, V1, solx, xini, xfinal)

# Drawing the Voronoi diagrams for the ground truth, initialization, and reconstruction
vorsol, istop = vf.voronoi(solx, sigma, Aeps, draw = True)
os.system('mpost voronoi.mp')
os.rename('./voronoi.mp', './solx.mp')
os.rename('./voronoi.mps', './solx.eps')

vorini, istop = vf.voronoi(xini, sigma, Aeps, draw = True)
os.system('mpost voronoi.mp')
os.rename('./voronoi.mp', './xini.mp')
os.rename('./voronoi.mps', './xini.eps')

vorfinal, istop = vf.voronoi(xfinal, sigma, Aeps, draw = True)
os.system('mpost voronoi.mp')
os.rename('./voronoi.mp', './xfinal.mp')
os.rename('./voronoi.mps', './xfinal.eps')


# Calculating the initial value of the cost function and the gradient, without weights
weightsG = np.ones(nsources)
finit1, ginit1 = evalfg(2*nsites, xini, vorini, greq = True)
gp = xini - ginit1
project(2*nsites, gp)
gp = gp - xini
normgpinit1 = np.linalg.norm(gp)


ffinal1, gfinal1 = evalfg(2*nsites, xfinal, vorfinal, greq = True)
gp = xfinal - gfinal1
project(2*nsites, gp)
gp = gp - xfinal
normgpfinal1 = np.linalg.norm(gp)

namefile = str(nsites)+'_'+str(nsources)+'_'+str(nmesh)+'_'+str(int(1000.0*noise_coeff))+'_'+str(ninit)

np.savez("./"+namefile+".npz", nsites=nsites, nsources=nsources, nmesh=nmesh, noise_coeff=noise_coeff, Aeps=Aeps, 
finit=finit, noise_level=noise_level, normgpinit=normgpinit, ffinal=ffinal, normgpfinal=normgpfinal, erroropt=erroropt, errorinit=errorinit, flagsol=flagsol, iter=iter, numevalf=numevalf, CPU_time=CPU_time, sigma=sigma, solx=solx, xini=xini, xfinal=xfinal, finit1=finit1, ffinal1=ffinal1, normgpinit1=normgpinit1, normgpfinal1=normgpfinal1)

    
print('End of main loop !!')