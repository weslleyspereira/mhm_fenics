##	Test meshes
#	
#	

from dolfin import *

# Set PETSc MUMPS paramter (this is required to prevent a memory error
# in some cases when using MUMPS LU solver).
if has_petsc():
	PETScOptions.set("mat_mumps_icntl_14", 40.0)

# Define boundary condition values
g = Constant(0.0)

# Create mesh and define function space
mesh = Mesh("../demos/meshes/refTetrahedra.xml.gz")
V = FunctionSpace(mesh, "Lagrange", 1)

for i in range(2):

	# Define variational problem
	u = TrialFunction(V)
	v = TestFunction(V)
	f = Constant(1.0)
	a = dot(grad(u), grad(v))*dx
	L = f*v*dx

	# Define boundary conditions
	bc1 = DirichletBC(V, g, 1)
	bc2 = DirichletBC(V, g, 2)
	bc3 = DirichletBC(V, g, 3)
	bc4 = DirichletBC(V, g, 4)

	# Compute solutions
	vtk_file = File("test_mesh"+str(i)+".pvd")
	u = Function(V)

	solve(a == L, u, bc1)
	vtk_file << u
	solve(a == L, u, bc2)
	vtk_file << u
	solve(a == L, u, bc3)
	vtk_file << u
	solve(a == L, u, bc4)
	vtk_file << u

	# Plot solution
	plot(u, mode="displacement", interactive=True)

	mesh = refine(mesh)
	V = FunctionSpace(mesh, "Lagrange", 1)
