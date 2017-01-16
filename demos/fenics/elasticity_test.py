"""
Example to solve an elasticity equation
"""

from dolfin import *

# Create mesh and define function space
mesh = UnitCubeMesh(16, 16, 16)
V = VectorFunctionSpace(mesh, "Lagrange", 2, 3)

# Define Dirichlet boundary (x = 0 or x = 1)
def boundary(x):
    return x[0] < DOLFIN_EPS or x[0] > 1.0 - DOLFIN_EPS

# Define boundary condition
u0 = Constant((0.0,0.0,0.0))
bc = DirichletBC(V, u0, boundary)

# Define variational problem
u = TrialFunction(V)
v = TestFunction(V)
f = Expression((
"pi*pi*(-4*1.*sin(pi*(-2*x[0] + 2*x[1] + 2*x[2])) + 2*1.*sin(pi*(2*x[0] - 2*x[1] + 2*x[2])) - 6*1.*sin(pi*(2*x[0] + 2*x[1] + 2*x[2])) - sin(pi*(-2*x[0] + 2*x[1] + 2*x[2])) + 5*sin(pi*(2*x[0] - 2*x[1] + 2*x[2])) + 3*sin(pi*(2*x[0] + 2*x[1] - 2*x[2])) - 9*sin(pi*(2*x[0] + 2*x[1] + 2*x[2])))",
"2*pi*pi*(2*1.*sin(pi*(-2*x[0] + 2*x[1] + 2*x[2])) - 1.*sin(pi*(2*x[0] - 2*x[1] + 2*x[2])) - 3*1.*sin(pi*(2*x[0] + 2*x[1] + 2*x[2])) + 5*sin(pi*(-2*x[0] + 2*x[1] + 2*x[2])) + 2*sin(pi*(2*x[0] - 2*x[1] + 2*x[2])) + 3*sin(pi*(2*x[0] + 2*x[1] - 2*x[2])) - 6*sin(pi*(2*x[0] + 2*x[1] + 2*x[2])))",
"pi*pi*(4*1.*sin(pi*(-2*x[0] + 2*x[1] + 2*x[2])) + 2*1.*sin(pi*(2*x[0] - 2*x[1] + 2*x[2])) - 6*1.*sin(pi*(2*x[0] + 2*x[1] + 2*x[2])) + 13*sin(pi*(-2*x[0] + 2*x[1] + 2*x[2])) + 11*sin(pi*(2*x[0] - 2*x[1] + 2*x[2])) + 9*sin(pi*(2*x[0] + 2*x[1] - 2*x[2])) - 15*sin(pi*(2*x[0] + 2*x[1] + 2*x[2])))"
), degree=2)
#g = Expression("0", degree=2)
h = Constant((0.0,0.0,0.0))

cteMu = Constant(1.0)
cteLambda = Constant(.25)

eps_u = sym(grad(u))
eps_v = sym(grad(v))

a = inner(2.*cteMu*eps_u, eps_v)*dx + inner(cteLambda*tr(eps_u),tr(eps_v))*dx
L = inner(f,v)*dx

# Compute solution
u_h = Function(V)
solve(a == L, u_h.vector(), bc, "lu")

# Save solution in VTK format
file = File("elasticity.pvd")
file << u

# Plot solution
plot(u, interactive=True)
