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
from dolfin.cpp.mesh import MeshFunction, Mesh
from dolfin.cpp.io import File
from dolfin.cpp.la import LUSolver
from demos.meshes.geometry import markBoundarySubdomains, markBoundariesOfMesh

#from mshr import *
#from math import ceil

###
### Controled by MHM
###

verbose=True
reutilizeBoundaryInfo=False

meshPath = "demos/meshes/triangle.7.xml.gz"
bndPath = "boundaries.xml.gz"
bndPvd = "bndSubdomains.pvd"
outputFolder = "./new_results/"

### Load model parameters
import demos.data.physics_elastodynamics2D_dirichlet as physics
### Load geometry
from demos.meshes.square1x1 import geo
### Load Finite Element Space params
from demos.params.galerkin import finiteElement
# Load Newmark parameters
import demos.params.newmark as newmark
# Time interval
from demos.params.timeInterval import timeInterval

###
### Configuration
###

### Load elasticity module
import demos.data.elasticity as elasticity

# Load mesh
if isfile(meshPath):
	mesh = Mesh(meshPath)
else:
	print("Mesh bndFile does not exist!")
	exit(1)
#plot(mesh, interactive=True)

###
### Mesh
###

# Here we consider 0 as a internal face
bndSubdomains = markBoundariesOfMesh(mesh, geo, outputPvd="localBoundary.pvd")
globalBndSubdomains = MeshFunction("size_t", mesh, geo._nDim-1)
	
# If there exist marked boundary
if 'boundaries' in dir(physics):
	if (reutilizeBoundaryInfo and isfile(bndPath)):
		globalBndSubdomains = MeshFunction("size_t", mesh, bndPath)
	else:
		globalBndSubdomains = markBoundarySubdomains(mesh, physics.boundaries, outputXmlGz=bndPath, outputPvd=bndPvd)
else: 
	globalBndSubdomains.set_all(0)
	
# Boundary measures
ds = Measure('ds', domain=mesh, subdomain_data=bndSubdomains)
dsGlobal = Measure('ds', domain=mesh, subdomain_data=globalBndSubdomains)

###
### Numerical configuration
###

# MHM parameters
psiL = 1

# Psi scalar strings
# psiString = [ "1.0" ]
# psiString = []
# if psiL == 1:
# 	for i in range(geo._nDim):
# 		psiString.append("x["+str(i)+"]")
psiString = [ "x[0]", "1-x[0]" ]

# Psi basis
psiBasis = []
for n in range(geo._nDim):
	for s in psiString:
		vecS = []
		
		for i in range(n):
			vecS.append("0.0")
		#print(vecS)
		
		vecS.append(s)
		#print(vecS)
		
		for i in range(geo._nDim-n-1):
			vecS.append("0.0")
		
		print(vecS)
		psiBasis.append( Expression(vecS, degree=2) )
		#Expression(vecS, degree=2)
		
# psiBasis = (
# 	Expression(("1.","0.","0."), degree=2),
# 	Expression(("0.","1.","0."), degree=2),
# 	Expression(("0.","0.","1."), degree=2),
# 	Expression(("-x[1]","x[0]","0."), degree=2),
# 	Expression(("-x[2]","0.","x[0]"), degree=2),
# 	Expression(("0.","-x[2]","x[1]"), degree=2)
# )

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
# basisOp = [ inner(psi, w)*ds(i) for psi in psiBasis for i in range(geo.getNoFaces()) ]
basisOp = [ inner(psi, w)*ds(1) for psi in psiBasis ]

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

outputFolder = "./new_results/"

pvdFileU = [ File(outputFolder+"/"+str(i)+"/u_psi.pvd", "compressed") for i in range(nPsi)]
pvdFileV = [ File(outputFolder+"/"+str(i)+"/v_psi.pvd", "compressed") for i in range(nPsi)]

###
### Solver
###

# Timestep
dt = (timeInterval.tf-physics.t0)/timeInterval.nt

# Dirichlet boundary condition
dirichletBC = {}
if 'dirichletBC' in dir(physics):
	for i in physics.dirichletBC:
		physics.dirichletBC[i].t = physics.t0
		ug = project(physics.dirichletBC[i], V)
		dirichletBC[i] = DirichletBC(V, ug, globalBndSubdomains, i)
		dirichletBC[i].apply(M)
		dirichletBC[i].apply(R)

# Preparing solvers
solverM = LUSolver(M)
solverMdtR = LUSolver(M + newmark.beta*dt*dt*R)
solverM.parameters['reuse_factorization'] = True
solverMdtR.parameters['reuse_factorization'] = True

# Psi loop
for iPsi in range(nPsi):

	if verbose: print("Psi: ", iPsi)
	
	# Initial condition
	tk = physics.t0
	u = project(physics.u0, V) # This can be improved by not projecting for each psi, but be aware! Do not use: u = u0_ini!!!
	v = project(physics.v0, V)

	# Neumann boundary condition
	timeNeumann = False
	neumannk = 0
	if 'neumannBC' in dir(physics):
		neumannOp = 0
		for i in physics.neumannBC:
			neumannOp += inner(physics.neumannBC[i], w)*dsGlobal[i]
			if isinstance(physics.neumannBC[i], Expression):
				if 't' in physics.neumannBC[i].user_parameters:
					timeNeumann = True
		neumannk = assemble(neumannOp)

	# Variational RHS
	pk = -FPsi[iPsi] + neumannk

	# Save solutions in VTK format
	print("Time: ", tk)
	pvdFileU[iPsi] << (u, tk)
	pvdFileV[iPsi] << (v, tk)

	for k in range(1, timeInterval.nt+1):

		t0 = tk
		tk = physics.t0 + k*dt
		print("Time: ", tk)
	
		# Some additional terms
		u_0 = u.vector()
		v_0 = v.vector()
		p0  = pk
	
		# Updating  dirichlet boundary condition
		for i in dirichletBC:
			if isinstance(physics.dirichletBC[i], Expression):
				if 't' in physics.dirichletBC[i].user_parameters:
					if verbose: print("Updating Dirichlet BC...")
					physics.dirichletBC[i].t = tk
					ug = project(physics.dirichletBC[i], V)
					dirichletBC[i] = DirichletBC(V, ug, globalBndSubdomains, i)
	
		# Updating neumann boundary condition
		if timeNeumann:
			if verbose: print("Updating Neumann BC...")
			neumannOp = 0
			for i in physics.neumannBC:
				neumannOp += inner(physics.neumannBC[i], w)*ds(i)
			neumannk = assemble(neumannOp)
	
		# Updating variational RHS
		pk = -FPsi[iPsi] + neumannk
	
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
		pvdFileU[iPsi] << (u, tk)
		pvdFileV[iPsi] << (v, tk)

## Plot solution
#plot(u, interactive=True)
