from __future__ import print_function
from fenics import *
#from mshr import *
#from math import ceil

# Load model parameters
from demos.data.physics_elastodynamics3D

# Stress tensor
def sigma(r):
	Eps = sym(grad(r))
	return 2.0*mu*Eps + lmbda*tr(Eps)*Identity(len(r))

###
### Mesh
###

# Topological dimension
nDim = 2

# Points used for face elements
geoPoints = (
	Point(0.0,0.0,0.0),
	Point(1.0,0.0,0.0),
	Point(0.0,1.0,0.0),
	Point(0.0,0.0,1.0)
)

# Number of faces elements
nFaces = 4

# Face elements
geoFaces = (
	(1,2,3),
	(0,2,3),
	(0,1,3),
	(0,1,2)
)

# Reference points in faces
geoBndPoints = [ face[0] for face in geoFaces ]

# Face normals
geoNormals = [
	(geoPoints[face[1]]-geoPoints[face[0]]).cross(geoPoints[face[2]]-geoPoints[face[0]]) for face in geoFaces
]

# Mesh
mesh = Mesh("demos/meshes/tet01.xml.gz")
#plot(mesh, interactive=True)

# Create mesh functions over the cell facets
boundary_sub_domains = MeshFunction("size_t", mesh, nDim-1)
dirichlet_sub_domains = MeshFunction("bool", mesh, nDim-1)

# Mark all facets of the domain as nFaces
boundary_sub_domains.set_all(nFaces)

# Mark all facets of the domain as False for Dirichlet
dirichlet_sub_domains.set_all(False)

# Mark boundary faces for face elements
for i in range(nFaces):
	class SubDomain(SubDomain):
		def inside(self, x, on_boundary):
			return near( geoNormals[i].dot(Point(x)-geoPoints[geoBndPoints[i]]), 0.0 ) and on_boundary
	aSubDomain = SubDomain()
	aSubDomain.mark(boundary_sub_domains, i)
	
# Boundary measures
ds = Measure('ds', domain=mesh, subdomain_data=boundary_sub_domains)

# Save Dolfin XML format of the subdomains
File("mesh_func.xml") << boundary_sub_domains

# Save sub domains to VTK files
file = File("subdomains.pvd")
file << boundary_sub_domains
file = File("dirichlet.pvd")
file << dirichlet_sub_domains

###
### Numerical configuration
###

## MHM parameters
#psiL 	= 1

## Psi scalar strings
#psiString = [ "1.0" ]
#if psiL == 1:
#	for i in range(nDim):
#		psiString.append("x["+str(i)+"]")

# Psi basis
#psiBasis = []
#for n in range(nDim):
#	for s in psiString:
#		vecS = []
#		for i in range(n):
#			vecS.append("0.0")
#		print(vecS)
#		vecS.append(s)
#		print(vecS)
#		for i in range(nDim-n-1):
#			vecS.append("0.0")
#		print(vecS)
#		#psiBasis.append( Expression(s, degree=2) )
#		Expression(vecS, degree=2)
psiBasis = (
	Expression(("1.","0.","0."), degree=2),
	Expression(("0.","1.","0."), degree=2),
	Expression(("0.","0.","1."), degree=2),
	Expression(("-x[1]","x[0]","0."), degree=2),
	Expression(("-x[2]","0.","x[0]"), degree=2),
	Expression(("0.","-x[2]","x[1]"), degree=2)
)

## Newmark parameters
newmarkGamma 	= 0.5
newmarkDelta 	= newmarkGamma-0.5
newmarkBeta 	= 0.75
nT				= 100

# Standard FEM parameters
femFamily 	= "CG"
femOrder 	= 3

# Defining function space
V = VectorFunctionSpace(mesh, femFamily, femOrder)

###
### Variational squeme
###

# Trial and test functions
u = TrialFunction(V)
v = TestFunction(V)

# weak operators
massOp = rho*inner(u, v)*dx
stiffOp = inner(sigma(u), grad(v))*dx
basisOp = [ inner(psi, v)*ds(i) for psi in psiBasis for i in range(nFaces) ]

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
