from fenics import Point

# Topological dimension
nDim = 3

# Points used for face elements
points = (
	Point(0.0,0.0,0.0),
	Point(1.0,0.0,0.0),
	Point(0.0,1.0,0.0),
	Point(0.0,0.0,1.0)
)

# Number of faces elements
nFaces = 4

# Face elements
faces = (
	(1,2,3),
	(0,2,3),
	(0,1,3),
	(0,1,2)
)
