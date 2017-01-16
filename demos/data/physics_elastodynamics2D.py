# Physical parameters
rho 	= 1.0
E 		= 1.0
nu 		= 0.25

# Source
f 		= ("2.4*pi*pi*sin(0.894427190999916*pi*t)*sin(2*pi*x[1])*cos(2*pi*x[0])", "-2.4*pi*pi*sin(0.894427190999916*pi*t)*sin(2*pi*x[0])*cos(2*pi*x[1])")

# Boundary conditions
g 		= {
			1 : (0.0, 0.0)
}
h 		= {
			1 : (0.0, 0.0)
			2 : (0.0, 0.0)
			3 : (0.0, 0.0)
			4 : (0.0, 0.0)
}

# Initial conditions
T0 		= 0.0
u0 		= (0.0,0.0)
v0 		= (0.0,0.0)
