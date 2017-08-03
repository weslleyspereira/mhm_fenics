from fenics import Point

# Topological dimension
nDim = 2

# Points used for face elements
points = (
	Point(0.0,0.0),
	Point(1.0,0.0),
	Point(1.0,1.0),
	Point(0.0,1.0)
)

# Number of faces elements
nFaces = 4

# Face elements
faces = (
	(0,1),
	(1,2),
	(2,3),
	(3,0)
)

# Face markers
facesMarkers = (
	1,
	2,
	3,
	4
)
