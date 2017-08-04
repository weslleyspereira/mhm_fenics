'''
Created on 3 de ago de 2017

@author: weslley
'''

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
f 		= Expression(("1.0*omega*omega*sin(2*pi*x[0])*sin(2*pi*x[1])*cos(omega*t) - pi*pi*(3.2*cos(pi*(2*x[0] - 2*x[1])) - 1.6*cos(pi*(2*x[0] + 2*x[1])))*(cos(omega*t) - 1) - 1.6*pi*pi*(cos(omega*t) - 1)*cos(pi*(2*x[0] - 2*x[1]))", "-1.0*omega*omega*sin(2*pi*x[0])*sin(2*pi*x[1])*cos(omega*t) + pi*pi*(3.2*cos(pi*(2*x[0] - 2*x[1])) - 1.6*cos(pi*(2*x[0] + 2*x[1])))*(cos(omega*t) - 1) + 1.6*pi*pi*(cos(omega*t) - 1)*cos(pi*(2*x[0] - 2*x[1]))"), degree = 2, t = 0.0, omega = omega)
#f 		= Expression((0.0,0.0))
#f 		= Constant((0.0, 0.0))

# Boundaries
#dirichletBnd = Expression("x[0] < -1", degree = 2)
#dirichletBnd = Expression("x[0] == 0", degree = 2)
boundaries = {
	1: Expression("x[1] == 0.0", degree = 2),
	2: Expression("x[0] == 1.0", degree = 2),
	3: Expression("x[1] == 1.0", degree = 2),
	4: Expression("x[0] == 0.0", degree = 2)
}

# Boundary conditions
dirichletBC = {
	1 : Constant((0.0, 0.0))#,
	#2 : Expression(("1.0","0"), degree = 2)
}
#neumannBC  = {
#	2 : Constant((0.0, 0.0))
#}
absorbingBC  = {
	4 : Constant((0.0, 0.0))
}

# Initial conditions
t0 		= 0.0
u0 		= Constant((0.0,0.0))
v0 		= Constant((0.0,0.0))
#v0 		= Expression(("omega*sin(2*pi*x[1])*cos(2*pi*x[0])",
#"-omega*sin(2*pi*x[0])*cos(2*pi*x[1])"), degree=2, omega = omega)

# Exact solution

u 		= Expression((
"(-cos(omega*t) + 1)*sin(2*pi*x[0])*sin(2*pi*x[1])",
"(cos(omega*t) - 1)*sin(2*pi*x[0])*sin(2*pi*x[1])"),
degree = 2, t = 0.0, omega = omega)
v		= Expression((
"omega*sin(2*pi*x[0])*sin(2*pi*x[1])*sin(omega*t)",
"-omega*sin(2*pi*x[0])*sin(2*pi*x[1])*sin(omega*t)"),
degree = 2, t = 0.0, omega = omega)
