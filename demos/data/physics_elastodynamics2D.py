from fenics import Expression, Constant
from numpy import sqrt, pi

# Physical parameters
rho 	= 1.0
E 		= 1.0
nu 		= 0.25
mu		= E / (2.0*(1.0 + nu))
lmbda 	= E*nu / ((1.0 + nu)*(1.0 - 2.0*nu))
omega 	= sqrt(2.0)*pi*sqrt(mu/rho)

# Source
f 		= Expression(("2.4*pi*pi*sin(omega*t)*sin(2*pi*x[1])*cos(2*pi*x[0])", "-2.4*pi*pi*sin(omega*t)*sin(2*pi*x[0])*cos(2*pi*x[1])"), degree = 2, t = 0.0, omega = omega)
#f 		= Expression((0.0,0.0))
f 		= Constant((0.0, 0.0))

# Boundaries
#dirichletBnd = Expression("x[0] < -1", degree = 2)
#dirichletBnd = Expression("x[0] == 0", degree = 2)
boundaries = {
	1: Expression("x[0] == 0", degree = 2),
	2: Expression("x[0]+x[1] > 0.9", degree = 2)
}

# Boundary conditions
dirichletBC = {
	1 : Expression(("0","0"), degree = 2)
}
neumannBC  = {
	2 : Constant((1.0, 1.0))
}

# Initial conditions
t0 		= 0.0
u0 		= Constant((0.0,0.0))
v0 		= Constant((0.0,0.0))
#v0 		= Expression(("omega*sin(2*pi*x[1])*cos(2*pi*x[0])",
#"-omega*sin(2*pi*x[0])*cos(2*pi*x[1])"), degree=2, omega = omega)

# Exact solution

u 		= Expression((
"sin(omega*t)*cos(2*pi*x[0])*sin(2*pi*x[1])",
"-sin(omega*t)*sin(2*pi*x[0])*cos(2*pi*x[1])"),
degree = 2, t = 0.0, omega = omega)
v		= Expression((
"omega*sin(2*pi*x[1])*cos(omega*t)*cos(2*pi*x[0])",
"-omega*sin(2*pi*x[0])*cos(omega*t)*cos(2*pi*x[1])"),
degree = 2, t = 0.0, omega = omega)
