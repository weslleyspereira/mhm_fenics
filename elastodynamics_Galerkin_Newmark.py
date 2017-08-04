'''
Created on 3 de ago de 2017

@author: weslley
'''

###
### Python modules
###

# python future
from __future__ import print_function

# os path
from os.path import isfile

# numpy
from numpy import array

# fenics
from fenics import dx

# ufl
from ufl.operators import inner, grad
from ufl.measure import Measure

# Dolfin functions
from dolfin.functions.functionspace import VectorFunctionSpace,\
	TensorFunctionSpace
from dolfin.functions.function import TrialFunction, TestFunction
from dolfin.functions.expression import Expression

# Dolfin fem
from dolfin.fem.projection import project
from dolfin.fem.bcs import DirichletBC
from dolfin.fem.assembling import assemble

# Dolfin cpp
from dolfin.cpp.mesh import MeshFunction, SubDomain, Mesh
from dolfin.cpp.io import File
from dolfin.cpp.la import LUSolver
from demos.meshes.geometry import markBoundarySubdomains

###
### Configuration
###

verbose=True
reutilizeBoundaryInfo=False

meshPath = "demos/meshes/triangle.7.xml.gz"
outputFolder = "./results_galerkin/"

### Load elasticity module
import demos.data.elasticity as elasticity
### Load model parameters
import demos.data.physics_elastodynamics2D_dirichlet as physics
### Load geometry
from demos.meshes.square1x1 import geo
#import demos.meshes.triangle1x1 as geo
### Load Finite Element Space params
from demos.params.galerkin import finiteElement
# Load Newmark parameters
import demos.params.newmark as newmark
# Time interval
from demos.params.timeInterval import timeInterval

###
### Load mesh
###

# Mesh paths
#dirichletPath = "dirichletBnd.xml.gz"
bndPath = "boundaries.xml.gz"
bndPvd = "bndSubdomains.pvd"

# Load mesh
if isfile(meshPath):
	mesh = Mesh(meshPath)
else:
	print("Mesh bndFile does not exist!")
	exit(1)

# Here we consider 0 as a internal face
#dirichletSubdomain = MeshFunction("bool", mesh, geo.nDim-1)
bndSubdomains = MeshFunction("size_t", mesh, geo._nDim-1)
	
# If there exist marked boundary
if 'boundaries' in dir(physics):
	if (reutilizeBoundaryInfo and isfile(bndPath)):
		bndSubdomains = MeshFunction("size_t", mesh, bndPath)
	else:
		bndSubdomains = markBoundarySubdomains(mesh, physics.boundaries, outputXmlGz=bndPath, outputPvd=bndPvd)
else: 
	bndSubdomains.set_all(0)

## Get the set of subdomains
#setneumannSubdomains = set(neumannSubdomains.array())
#nBoundarySubdomains = len(setneumannSubdomains)-1 #exclude 0 labels
	
# Define measures
ds = Measure('ds', domain=mesh, subdomain_data=bndSubdomains) # boundary measure

# Plot mesh
#plot(mesh, interactive=True)
#plot(bndSubdomains, interactive=True)
#plot(dirichletSubdomain, interactive=True)

###
### Variational scheme
###

# Defining displacement space
V = VectorFunctionSpace(mesh, finiteElement.family, finiteElement.order)
# Stresss space
S = TensorFunctionSpace(mesh, finiteElement.family, finiteElement.order)

# Trial and test functions
uTrial = TrialFunction(V)
w = TestFunction(V)

# weak operators
massOp = physics.rho*inner(uTrial, w)*dx
stiffOp = inner(elasticity.stress(physics.mu, physics.lmbda, uTrial), grad(w))*dx

# Assembling
M = assemble(massOp)
R = assemble(stiffOp)

## Create zero boundary condition
#bc = DirichletBC(V, Constant(geo.nDim*[0.0]), dirichletSubdomain)

## Construct dirichlet contribution vectors
#dirichletM = PETScVector(MPI_COMM_WORLD, M.size(0))
#bc.zero(M)
#bc.zero_columns(M, dirichletM, 1)
#dirichletR = PETScVector(MPI_COMM_WORLD, R.size(0))
#bc.zero(R)
#bc.zero_columns(R, dirichletR, 1)

###
### Output configuration
###

pvdFileU = File(outputFolder+"/u.pvd", "compressed")
pvdFileV = File(outputFolder+"/v.pvd", "compressed")
pvdFileStress = File(outputFolder+"/stress.pvd", "compressed")

## Save sub domains to VTK files
#bndFile = File(outputFolder+"/subdomains.pvd")
#bndFile << neumannSubdomains

###
### Solver
###

# Class representing the initial conditions
u0 = project(physics.u0, V)
v0 = project(physics.v0, V)

# Initial condition
tk = physics.t0
u = u0
v = v0

# Source
physics.f.t = tk
timef = False
fk = assemble(inner(physics.f, w)*dx)
if isinstance(physics.f, Expression):
	if 't' in physics.f.user_parameters: timef = True

# Neumann boundary condition
timeNeumann = False
neumannk = 0
if 'neumannBC' in dir(physics):
	neumannOp = 0
	for i in physics.neumannBC:
		neumannOp += inner(physics.neumannBC[i], w)*ds(i)
		if isinstance(physics.neumannBC[i], Expression):
			if 't' in physics.neumannBC[i].user_parameters:
				timeNeumann = True
	neumannk = assemble(neumannOp)

# Variational RHS
pk = fk+neumannk

# Dirichlet boundary condition
dirichletBC = {}
if 'dirichletBC' in dir(physics):
	for i in physics.dirichletBC:
		physics.dirichletBC[i].t = tk
		ug = project(physics.dirichletBC[i], V)
		dirichletBC[i] = DirichletBC(V, ug, bndSubdomains, i)
		dirichletBC[i].apply(M)
		dirichletBC[i].apply(R)

# Timestep
dt = (timeInterval.tf-physics.t0)/timeInterval.nt

# Preparing linear solvers
solverM = LUSolver(M)
solverMdtR = LUSolver(M + newmark.beta*dt*dt*R)
solverM.parameters['reuse_factorization'] = True
solverMdtR.parameters['reuse_factorization'] = True

# Save solutions in VTK format
print("Time: ", tk)
u.rename("u", "u")
pvdFileU << (u, tk)
v.rename("v", "v")
pvdFileV << (v, tk)
s = project(elasticity.stress(physics.mu, physics.lmbda, u), S)
s.rename("s","s")
pvdFileStress << (s, tk)

for k in range(1, timeInterval.nt+1):

	t0 = tk
	tk = physics.t0 + k*dt
	print("Time: ", tk)

	# Some additional terms
	u_0  = u.vector()
	v_0  = v.vector()
	p0   = pk
	
	# Updating  dirichlet boundary condition
	for i in dirichletBC:
		if isinstance(physics.dirichletBC[i], Expression):
			if 't' in physics.dirichletBC[i].user_parameters:
				if verbose: print("Updating Dirichlet BC...")
				physics.dirichletBC[i].t = tk
				ug = project(physics.dirichletBC[i], V)
				dirichletBC[i] = DirichletBC(V, ug, bndSubdomains, i)
	
	# Updating source
	if timef:
		if verbose: print("Updating f...")
		physics.f.t = tk
		fk = assemble(inner(physics.f, w)*dx)
	
	# Updating neumann boundary condition
	if timeNeumann:
		if verbose: print("Updating Neumann BC...")
		neumannOp = 0
		for i in physics.neumannBC:
			neumannOp += inner(physics.neumannBC[i], w)*ds(i)
		neumannk = assemble(neumannOp)
	
	# Updating variational RHS
	pk = fk+neumannk
	
	# RHS for the first system
	b1 = M*(u_0 + dt*v_0) \
		+ R*(dt*dt*(newmark.beta-0.5)*u_0) \
		+ dt*dt*( newmark.beta*pk + (0.5-newmark.beta)*p0 )
	
	# Apply boundary condition
	if len(dirichletBC) > 0:
		if verbose: print("Applying boundary condition...")
		for i in dirichletBC: dirichletBC[i].apply(b1)

	# Solve first system
	solverMdtR.solve( u.vector(), b1 )
	
	# RHS for the second system
	b2 = M*v_0 + R*(-dt*( newmark.gamma*u.vector() + (1-newmark.gamma)*u_0 )) \
		+ dt*( newmark.gamma*pk + (1-newmark.gamma)*p0 )
	
	# Apply boundary condition
	if len(dirichletBC) > 0:
		if verbose: print("Applying boundary condition...")
		for i in dirichletBC: dirichletBC[i].apply(b2)

	# Solve second system
	solverM.solve( v.vector(), b2 )

	# Save solutions in VTK format
	pvdFileU << (u, tk)
	pvdFileV << (v, tk)
	s = project(elasticity.stress(physics.mu, physics.lmbda, u), S)
	s.rename("s","s")
	pvdFileStress << (s, tk)
