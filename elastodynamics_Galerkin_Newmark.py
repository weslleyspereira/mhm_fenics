from __future__ import print_function
from fenics import *
import numpy as np
#from mshr import *
#from math import ceil

# Load model parameters
from demos.data.physics_elastodynamics2D import *

# Stress tensor
def sigma(r):
	Eps = sym(grad(r))
	return 2.0*mu*Eps + lmbda*tr(Eps)*Identity(len(r))

###
### Load mesh
###

mesh = Mesh("demos/meshes/tet01.xml.gz")
nDim = mesh.mesh.topology().dim()

# Here we consider 0 as a internal face
faceSubdomains = MeshFunction("size_t", mesh, "mesh_func.xml")
#dirichlet_sub_domains = MeshFunction("bool", mesh, nDim-1)

## Plot mesh
#plot(mesh, interactive=True)

## Save Dolfin XML format of the subdomains
#File("mesh_func.xml") << faceSubdomains

# Save sub domains to VTK files
file = File("subdomains.pvd")
file << faceSubdomains
#file = File("dirichlet.pvd")
#file << dirichlet_sub_domains

# Get the set of subdomains
setFaceSubdomains = set(faceSubdomains.array())
nBoundarySubdomains = len(setFaceSubdomains)-1 #exclude 0 labels
	
# Define measures
ds = Measure('ds', domain=mesh, subdomain_data=faceSubdomains)

###
### Variational squeme
###

# Load Newmark parameters
import demos.params.newmark as newmark
# Load Finite Element Space params
import demos.params.galerkin as fem

# Defining function space
V = VectorFunctionSpace(mesh, fem.family, fem.order)

# Trial and test functions
u = TrialFunction(V)
v = TestFunction(V)

# weak operators
massOp = rho*inner(u, v)*dx
stiffOp = inner(sigma(u), grad(v))*dx
basisOp = [ inner(psi, v)*ds(i) for psi in psiBasis for i in range(nFacets) ]

## Matrices and vectors
#M = PETScMatrix()
#R = PETScMatrix()
#P = [ PETScVector() for i in range(len(basisOp)) ]

## Assembling
#assemble(massOp, tensor=M)
#assemble(stiffOp, tensor=R)
#[ assemble(basisOp[i], tensor=P[i]) for i in range(len(P)) ]

# Assembling
M = assemble(massOp)
R = assemble(stiffOp)
FPsi = [ assemble(op) for op in basisOp ]
nPsi = len(FPsi)

###
### Output configuration
###

outputFolder = "./results/"

pvdFileU = [ File(outputFolder+"/u_psi_"+str(i)+".pvd") for i in range(nPsi)]
pvdFileV = [ File(outputFolder+"/v_psi_"+str(i)+".pvd") for i in range(nPsi)]

###
### Solver
###

nT	= 100
T = 1.0

dt 	= (T-T0)/nT

# Preparing solvers
solverM = LUSolver(M)
solverMdtR = LUSolver(M + newmarkBeta*dt*dt*R)
solverM.parameters['reuse_factorization'] = True
solverMdtR.parameters['reuse_factorization'] = True

# Class representing the intial conditions
u_0ini = project(u0, V)
v_0ini = project(v0, V)

# Psi loop
for i in range(nPsi):

	print("Psi: ", i)
	
	u = u_0ini
	v = v_0ini

	# Save solutions in VTK format
	print("Time: ", T0)
	pvdFileU[i] << u
	pvdFileV[i] << v

	tk = T0
	for k in range(1, nT+1):

		t0 = tk
		tk = T0 + k*dt
		print("Time: ", tk)
	
		# Some additional terms
		u_0 	= u.vector()
		v_0 	= v.vector()
		Mv0 	= M*v_0
	
		# Solve first system
		solverMdtR.solve( u.vector(),
			M*u_0 + dt*Mv0 + 
			dt*dt*(newmarkBeta-0.5)*R*u_0 + 
			dt*dt*(newmarkBeta*(-(tk-T0))*FPsi[i] + (0.5-newmarkBeta)*(-(t0-T0))*FPsi[i])
		)
	
		# Solve second system
		solverM.solve( v.vector(),
			Mv0 - 
			dt*R*( newmarkGamma*u.vector() + (1-newmarkGamma)*u_0 ) + 
			dt*( newmarkGamma*(-(tk-T0))*FPsi[i] + (1-newmarkGamma)*(-(t0-T0))*FPsi[i])
		)

		# Save solutions in VTK format
		pvdFileU[i] << u
		pvdFileV[i] << v

## Plot solution
#plot(u, interactive=True)
