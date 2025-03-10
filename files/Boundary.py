from fenics import *
import numpy as np
import time, random, sys, os
from prettytable import PrettyTable
import vorofunction as vf
       
class Left(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and near(x[0], 0.0, DOLFIN_EPS)

class Right(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and near(x[0], 1.0, DOLFIN_EPS)

class Bottom(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and near(x[1], 0.0, DOLFIN_EPS)

class Top(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and near(x[1], 1.0, DOLFIN_EPS)

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

        if ic[0] <= Aeps or ic[0] >= 1.0 - Aeps or ic[1] <= Aeps or ic[1] >= 1.0 - Aeps:
            kmin = nsites
        else:
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

def psiv(a, v, w):
    return Expression("-((x[0] - w0)*aw0 + (x[1] - w1)*aw1) / denpsiv", degree=2, w0 = w[0], w1 = w[1], aw0 = (rotate(a - w))[0], aw1 = (rotate(a - w))[1], denpsiv = np.inner(a - v, rotate(a - w)))

def gradpsiv(a, v, w):
    den = np.inner(a - v, rotate(a - w))
    if near(den, 0.0, DOLFIN_EPS):
        print('In gradpsiv, den is close to zero.')
    return - rotate(a - w) / den

def Mvl(v, gradphil, aj, ai):
    detval = np.linalg.det(np.array([aj - ai, gradphil]))
    if near(detval, 0.0, DOLFIN_EPS):
        print('In Mvl, detval is close to zero')
    return - np.outer(rotate(gradphil), v - ai) / detval

def Mv(v,ai, aj, ak):
    detval = np.linalg.det(np.array([aj - ai, ak - ai]))
    if near(detval, 0.0, DOLFIN_EPS):
        print('In Mv, detval is close to zero')
    return np.outer(rotate(ai - aj), (v - ak)) / detval

def tcoeff(p0, p1, p2):
    tcoeffval = np.zeros((3,3))

    T = np.array([p0, p1, p2])
    
    for i in range(3):
        ip1 = (i + 1) % 3
        ip2 = (i + 2) % 3
    
        a = T[ip1,1] - T[i,1]
        b = T[i,0] - T[ip1,0]
        c = (T[i,1] - T[ip1,1]) * T[i,0] + (T[ip1,0] - T[i,0]) * T[i,1]

        tcoeffval[i,:] = [a,b,c]

        if a * T[ip2,0] + b * T[ip2,1] + c <= 0.0:
            tcoeffval[i,:] = - tcoeffval[i,:]

    return tcoeffval

def triangintdom(site, v, vprev, vnext):
    # # Coefficients of the straight lines that determine a triangle.
    tcoeffvalprev = tcoeff(site,v,vprev)
    a0 = tcoeffvalprev[0,0]; b0 = tcoeffvalprev[0,1]; c0 = tcoeffvalprev[0,2]
    a1 = tcoeffvalprev[1,0]; b1 = tcoeffvalprev[1,1]; c1 = tcoeffvalprev[1,2]
    a2 = tcoeffvalprev[2,0]; b2 = tcoeffvalprev[2,1]; c2 = tcoeffvalprev[2,2]

    tcoeffvalnext = tcoeff(site,v,vnext)
    A0 = tcoeffvalnext[0,0]; B0 = tcoeffvalnext[0,1]; C0 = tcoeffvalnext[0,2]
    A1 = tcoeffvalnext[1,0]; B1 = tcoeffvalnext[1,1]; C1 = tcoeffvalnext[1,2]
    A2 = tcoeffvalnext[2,0]; B2 = tcoeffvalnext[2,1]; C2 = tcoeffvalnext[2,2]

    for cell in cells(mesh):
        p1 = np.array([cell.get_vertex_coordinates()[0], cell.get_vertex_coordinates()[1]])
        p2 = np.array([cell.get_vertex_coordinates()[2], cell.get_vertex_coordinates()[3]])
        p3 = np.array([cell.get_vertex_coordinates()[4], cell.get_vertex_coordinates()[5]])
        ic = incenter(p1, p2, p3)

        if a0 * ic[0] + b0 * ic[1] + c0 >= 0.0 and a1 * ic[0] + b1 * ic[1] + c1 >= 0.0 and a2 * ic[0] + b2 * ic[1] + c2 >= 0.0:
            triangdomains.array()[cell.index()] = 1
        elif A0 * ic[0] + B0 * ic[1] + C0 >= 0.0 and A1 * ic[0] + B1 * ic[1] + C1 >= 0.0 and A2 * ic[0] + B2 * ic[1] + C2 >= 0.0:
            triangdomains.array()[cell.index()] = 2
        else: 
            triangdomains.array()[cell.index()] = 0

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
        hsolve_clean = Function(W)
        hsolve = Function(V1)
        (htrial, zeta_1) = TrialFunctions(W)
        (vtest, zeta_2) = TestFunctions(W) 

        # Define PDE
        a = 0.0
        for i in range(nsites+1):
            a = a + (inner(Constant(sigma[i]) * grad(htrial), grad(vtest)) + zeta_1 * vtest + htrial * zeta_2) * dx(i) 
        L = (gvec_list[alpha])[0] * vtest * ds(1) + (gvec_list[alpha])[1] * vtest * ds(2) + (gvec_list[alpha])[2] * vtest * ds(3) + (gvec_list[alpha])[3] * vtest * ds(4)

        # Compute solution of (2)-(3)
        solve(a == L, hsolve_clean)
        (hsolve_temp, csolve) = hsolve_clean.split(True)

        # add noise to the data
        max_hsolve = np.abs(hsolve_temp.vector()[:]).max()
        h_perturb = np.random.default_rng(seed=123456).normal(loc= 0, scale=noise_coeff*max_hsolve, size=hsolve.vector()[:].size)
        hsolve.vector()[:] = hsolve_temp.vector()[:] + h_perturb

        noise_num = noise_num + zeta_alpha*assemble( (hsolve - hsolve_temp)**2 * ds) 
        noise_denom = noise_denom +  zeta_alpha*assemble( hsolve_temp**2 * ds )
    
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

def evalfg(n, x, vor, greq=False, gbound=True, gdistr=False):
    nsites = int(n/2)
    sites = np.zeros((nsites,2))
    for i in range(nsites):
        sites[i,:] = x[2*i:2*i+2]
    
    omega(x)
    dx = Measure("dx", domain=mesh, subdomain_data=domains)

    int_boundary.set_all(0)
    tag_internal_boundary()
    dS_int = Measure("dS", domain=mesh, subdomain_data=int_boundary)

    G = 0.0
    gradGbound = np.zeros(2*nsites)
    gradGdistr = np.zeros(2*nsites)

    for alpha in range(nsources):
 
        # Define rhs and lhs of PDE
        a = 0.0
        for i in range(nsites+1):
            a = a + (Constant(sigma[i]) * inner(grad(utf), grad(wtest))) * dx(i)
        L = (gvec_list[alpha])[0] * wtest * ds(1) + (gvec_list[alpha])[1] * wtest * ds(2)

        # Assemble matrices
        aa = assemble(a)
        LL = assemble(L)

        # Dirichlet boundary conditions, applied on top and bottom of square (Gamma_a)
        DirichletBC(V1,h[alpha],bottom).apply(aa)
        DirichletBC(V1,h[alpha],bottom).apply(LL)
        DirichletBC(V1,h[alpha],top).apply(aa)
        DirichletBC(V1,h[alpha],top).apply(LL)
  
        # Solve
        solve(aa, usol.vector(), LL)

        # Define rhs and lhs of PDE
        a = 0.0
        for i in range(nsites+1):
            a = a + (Constant(sigma[i]) * inner(grad(vtf), grad(wtest))) * dx(i)
        L = (gvec_list[alpha])[2] * wtest * ds(3) + (gvec_list[alpha])[3] * wtest * ds(4)

        # Assemble matrices
        aa = assemble(a)
        LL = assemble(L)

        # Dirichlet boundary conditions, applied on left and right of square (Gamma_b)
        DirichletBC(V1,h[alpha],left).apply(aa)
        DirichletBC(V1,h[alpha],left).apply(LL)
        DirichletBC(V1,h[alpha],right).apply(aa)
        DirichletBC(V1,h[alpha],right).apply(LL)

        # Solve
        solve(aa, vsol.vector(), LL)

        # Compute cost function        
        Galpha = assemble( (usol - vsol)**2 * dx )
        if weightsG[alpha] == -1.0:
            weightsG[alpha] = 1.0/Galpha
    
        G = G + 0.5 * weightsG[alpha] * Galpha
        
        if greq:
            # Define rhs and lhs of PDE
            a = 0.0
            for i in range(nsites+1):
                a = a + (Constant(sigma[i]) * inner(grad(ptf), grad(wtest))) * dx(i)
            L = - (usol - vsol) * wtest * dx

            # Assemble matrices
            aa = assemble(a)
            LL = assemble(L)

            # Dirichlet boundary conditions, applied on top and bottom of square (Gamma_a)
            DirichletBC(V1,0.0,bottom).apply(aa)
            DirichletBC(V1,0.0,bottom).apply(LL)
            DirichletBC(V1,0.0,top).apply(aa)
            DirichletBC(V1,0.0,top).apply(LL)

            #Solve
            solve(aa, psol.vector(), LL)

            # Define rhs and lhs of PDE
            a = 0.0
            for i in range(nsites+1):
                a = a + (Constant(sigma[i]) * inner(grad(qtf), grad(wtest))) * dx(i)
            L = (usol - vsol) * wtest * dx

            # Assemble matrices
            aa = assemble(a)
            LL = assemble(L)

            # Dirichlet boundary conditions, applied on left and right of square (Gamma_b)
            DirichletBC(V1,0.0,left).apply(aa)
            DirichletBC(V1,0.0,left).apply(LL)
            DirichletBC(V1,0.0,right).apply(aa)
            DirichletBC(V1,0.0,right).apply(LL)

            # Solve
            solve(aa, qsol.vector(), LL)

            # # Define rhs and lhs of PDE
            a = 0.0
            for i in range(nsites+1):
                a = a + inner(grad(Htrial), grad(vhtest)) * dx(i)
            L = (vhtest-vhtest)* dx

            # Assemble matrices
            aa = assemble(a)
            LL = assemble(L)

            # Dirichlet boundary conditions, applied on top and bottom of square (Gamma_a)
            DirichletBC(V2,h[alpha],bottom).apply(aa)
            DirichletBC(V2,h[alpha],bottom).apply(LL)
            DirichletBC(V2,h[alpha],top).apply(aa)
            DirichletBC(V2,h[alpha],top).apply(LL)
            
            # Dirichlet boundary conditions, applied on left and right of square (Gamma_b)
            DirichletBC(V2,h[alpha],left).apply(aa)
            DirichletBC(V2,h[alpha],left).apply(LL)
            DirichletBC(V2,h[alpha],right).apply(aa)
            DirichletBC(V2,h[alpha],right).apply(LL)
            
            # Solve
            solve(aa, Hsol.vector(), LL)

            for i in range(nsites):

                if gdistr:
                    
                    ### Volumetric Approach ###
                    temp = 0.5*(usol - h[alpha])**2 - fsource[alpha]*psol + inner(grad(usol), grad(psol)) + Constant(sigma[i])*usol*psol

                    # Warning: In case of non-constant fsource and Volumetric Approach, the derivative of fsource must be calculated explicitly.
                    # S0 = (h[alpha] - usol)*grad(h[alpha]) - grad(fsource[alpha])
                    S0 = (h[alpha] - usol)*grad(h[alpha])

                    S1 = Identity(2)*temp
                    S1 = S1 - outer(grad(psol),grad(usol)) - outer(grad(usol),grad(psol))
                    
                if gbound:
                    
                    grad_usol = grad(usol)
                    grad_vsol = grad(vsol)
                    grad_psol = grad(psol)
                    grad_qsol = grad(qsol)
                    grad_Hsol = grad(Hsol)

                    # the terms that are continuous can be removed because they cancel out at the interface of two cells

                    # tempplus = 0.5*(usol - vsol)**2  + Constant(sigma[i])*(inner(grad_usol('+'), grad_psol('+')) + inner(grad_vsol('+'), grad_qsol('+')))

                    tempplus = Constant(sigma[i])*(inner(grad_usol('+'), grad_psol('+')) + inner(grad_vsol('+'), grad_qsol('+')))

                    S1plus = Identity(2)*tempplus 

                    S1plus = S1plus - outer(grad_psol('+'),Constant(sigma[i])*grad_usol('+')) - outer(Constant(sigma[i])*(grad_usol('+') - grad_Hsol('+')), grad_psol('+')) - outer(grad_qsol('+'),Constant(sigma[i])*grad_vsol('+')) - outer(Constant(sigma[i])*(grad_vsol('+') - grad_Hsol('+')),grad_qsol('+'))
                
                    # the terms that are continuous can be removed because they cancel out at the interface of two cells
                    
                    # tempminus = 0.5*(usol - vsol)**2  + Constant(sigma[i])*(inner(grad_usol('-'), grad_psol('-')) + inner(grad_vsol('-'), grad_qsol('-')))

                    tempminus = Constant(sigma[i])*(inner(grad_usol('-'), grad_psol('-')) + inner(grad_vsol('-'), grad_qsol('-')))

                    S1minus = Identity(2)*tempminus
                    
                    S1minus = S1minus - outer(grad_psol('-'),Constant(sigma[i])*grad_usol('-')) - outer(Constant(sigma[i])*(grad_usol('-') - grad_Hsol('-')), grad_psol('-')) - outer(grad_qsol('-'),Constant(sigma[i])*grad_vsol('-')) - outer(Constant(sigma[i])*(grad_vsol('-') - grad_Hsol('-')),grad_qsol('-'))

                for r in range(len(vor[i])):
                    v = ((vor[i])[r])[0]
                    vflag = ((vor[i])[r])[1]
                    vnext = ((vor[i])[(r+1)%len(vor[i])])[0]

                    ### Boundary Expression ###
                    if gbound and (vflag >= 0):
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

                    ### Volumetric Approach ###
                    if gdistr and (vprevflag >= 0 or vflag >= 0):    
                        # Tag the integration domain
                        triangintdom(sites[i,:],v,vprev, vnext)
                        dt = Measure('dx')(subdomain_data=triangdomains)                  
                        
                        psivvalvprev = psiv(sites[i,:], v, vprev)
                        gradpsivvalTemp = gradpsiv(sites[i,:],v,vprev)
                        gradpsivval = Constant((gradpsivvalTemp[0], gradpsivvalTemp[1]))
                  
                        Tivw0 = assemble(((S1*gradpsivval + S0*psivvalvprev)[0])*dt(1))
                        Tivw1 = assemble(((S1*gradpsivval + S0*psivvalvprev)[1])*dt(1))
                        
                        psivvalvnext = psiv(sites[i,:], v, vnext)
                        gradpsivvalTemp = gradpsiv(sites[i,:],v,vnext)
                        gradpsivval = Constant((gradpsivvalTemp[0], gradpsivvalTemp[1]))                  
                    
                        Tivw0 = Tivw0 + assemble(((S1*gradpsivval + S0*psivvalvnext)[0])*dt(2))
                        Tivw1 = Tivw1 + assemble(((S1*gradpsivval + S0*psivvalvnext)[1])*dt(2))

                        if (vprevflag < 0 or vflag < 0):
                            j = max(vprevflag, vflag)
                            l = int(- min(vprevflag, vflag))

                            # Definindo qual vizinho está no intervalo de integração
                            if vflag >= 0:
                                w = vprev
                            else:
                                w = vnext

                            gradpsivvalTemp = gradpsiv(sites[i,:],v,w)
                            gradpsivval = Constant((gradpsivvalTemp[0], gradpsivvalTemp[1])) 

                            mvl = Mvl(v,gradphi[l],sites[j,:],sites[i,:])
                            gradGdistr[2*i:2*i+2] = gradGdistr[2*i:2*i+2] + weightsG[alpha] * np.matmul(mvl.T,[Tivw0, Tivw1])
                            
                            mvl = Mvl(v,gradphi[l],sites[i,:],sites[j,:])
                            gradGdistr[2*j:2*j+2] = gradGdistr[2*j:2*j+2] + weightsG[alpha] * np.matmul(mvl.T,[Tivw0, Tivw1])
                        else:
                            j = vprevflag; k = vflag #it could also be j = vflag; k = vprevflag 

                            mvjki = Mv(v,sites[j,:],sites[k,:], sites[i,:])
                            gradGdistr[2*i:2*i+2] = gradGdistr[2*i:2*i+2] + weightsG[alpha] * np.matmul(mvjki.T,[Tivw0, Tivw1])
                            
                            mvkij = Mv(v,sites[k,:],sites[i,:], sites[j,:])
                            gradGdistr[2*j:2*j+2] = gradGdistr[2*j:2*j+2] + weightsG[alpha] * np.matmul(mvkij.T,[Tivw0, Tivw1])

                            mvijk = Mv(v,sites[i,:],sites[j,:], sites[k,:])
                            gradGdistr[2*k:2*k+2] = gradGdistr[2*k:2*k+2] + weightsG[alpha] * np.matmul(mvijk.T,[Tivw0, Tivw1])
                            
                    vprev = v
                    vprevflag = vflag 
                             
    if greq:
        return G, gradGbound
    else:
        return G
    
    # if greq:
    #     if gdistr == True and gbound == True:
    #         return G, gradGdistr, gradGbound
    #     if gdistr == True and gbound == False:
    #         return G, gradGdistr
    #     if gdistr == False and gbound == True:
    #         return G, gradGbound
    # else:
    #     return G

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
    p[0] = max(Aeps, min(p[0], 1.0 - Aeps))
    p[1] = max(Aeps, min(p[1], 1.0 - Aeps))

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
# Build function space with Lagrange multiplier
P1 = FiniteElement("Lagrange", mesh.ufl_cell(), 1)
V1 = FunctionSpace(mesh, P1)
R = FiniteElement("Real", mesh.ufl_cell(), 0)
W = FunctionSpace(mesh, P1 * R)
P2 = FiniteElement("Lagrange", mesh.ufl_cell(), 2)
V2 = FunctionSpace(mesh, P2)

# Define domains
domains = MeshFunction("size_t", mesh, mesh.topology().dim(), 0)

# Define internal interface domain
int_boundary = MeshFunction("size_t", mesh, mesh.topology().dim()-1, 0)

# Definindo o domínio triangular de integração
triangdomains = MeshFunction("size_t", mesh, mesh.topology().dim(), 0)

#Funções phi_l que definem o domínio D: left-> phi_1 = -x[0], right-> phi_2 = x[0] - 1, top-> phi_3 = -x[1], bottom-> phi_4 = x[1] - 1 
#Gradientes:

gradphi_1 = [-1.0, 0.0]
gradphi_2 = [1.0, 0.0]
gradphi_3 = [0.0, -1.0]
gradphi_4 = [0.0, 1.0]
gradphi = [[0.0, 0.0], gradphi_1, gradphi_2, gradphi_3, gradphi_4]

#----------------
# Data
#----------------
sigma = np.ones(nsites+1)
if equidist:
    for l in range(nsites):
        sigma[l] = l + 1
if binary:
    sigma[0:int(np.ceil(nsites/2))] = [5.0]
    sigma[int(np.ceil(nsites/2)):nsites] = [10.0]
if ternary:
    ternary_values = np.array([3.0, 6.0, 9.0])
    for l in range(nsites):
        for i in range(3):
            if l % 3 == i:
                sigma[l] = ternary_values[i]


gvec_list = [as_vector([Constant("1.0"), Constant("1.0"), Constant("-1.0"), Constant("-1.0")]), as_vector([Constant("1.0"), Constant("-1.0"), Constant("1.0"), Constant("-1.0")]), as_vector([Constant("1.0"), Constant("-1.0"), Constant("-1.0"), Constant("1.0")])]

#----------------
# Ground truth data
#---------------- 
namefile = str(nsites)+'_'+str(nsources)+'_'+str(nmesh)+'_'+str(int(1000.0*noise_coeff))+'_'+str(ninit)
solx = (np.load("./"+namefile+".npz"))['solx']

# Checking if it was possible to construct the Voronoi diagram from the provided Ground Truth.
vor, istop = vf.voronoi(solx, sigma, Aeps, draw = False)
if istop != 0:
    print('The Voronoi method encountered an error while constructing the manufactured solution')
    sys.exit()
    
#  Tag subdomains
omega(solx)

dx = Measure("dx", domain=mesh, subdomain_data=domains)

# Tag boundaries of A
left   = Left()
top    = Top()
right  = Right()
bottom = Bottom()

aboundaries = MeshFunction("size_t", mesh, mesh.topology().dim() - 1)

left.mark(  aboundaries, 1)
right.mark( aboundaries, 2)
top.mark(   aboundaries, 3)
bottom.mark(aboundaries, 4)

ds = Measure("ds")(subdomain_data=aboundaries)

h = []
g = []
for alpha in range(nsources):

    hsolve_clean = Function(W)
    hsolve = Function(V1)
    (htrial, zeta_1) = TrialFunctions(W)
    (vtest, zeta_2) = TestFunctions(W) 

    # Define PDE
    a = 0.0
    for i in range(nsites+1):
        a = a + (inner(Constant(sigma[i]) * grad(htrial), grad(vtest)) + zeta_1 * vtest + htrial * zeta_2) * dx(i) 
    L = (gvec_list[alpha])[0] * vtest * ds(1) + (gvec_list[alpha])[1] * vtest * ds(2) + (gvec_list[alpha])[2] * vtest * ds(3) + (gvec_list[alpha])[3] * vtest * ds(4)

    # Compute solution
    solve(a == L, hsolve_clean)
    (hsolve_temp, csolve) = hsolve_clean.split(True)

    # add noise to the data
    max_hsolve = np.abs(hsolve_temp.vector()[:]).max()
    h_perturb = np.random.default_rng(seed=123456).normal(loc= 0, scale=noise_coeff*max_hsolve, size=hsolve.vector()[:].size)
    hsolve.vector()[:] = hsolve_temp.vector()[:] + h_perturb

    # save data
    h.append(hsolve)


# Define trial functions, test functions and solutions
utf = TrialFunction(V1)
vtf = TrialFunction(V1)
ptf = TrialFunction(V1)
qtf = TrialFunction(V1)
wtest = TestFunction(V1)
usol = Function(V1)
vsol = Function(V1)
psol = Function(V1)
qsol = Function(V1)
Hsol = Function(V2)
Htrial = TrialFunction(V2)
vhtest = TestFunction(V2)

#Initial parameter for projected gradient
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