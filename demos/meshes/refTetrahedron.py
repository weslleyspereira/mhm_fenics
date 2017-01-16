from fenics import Point

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
nFacets = 4

# Face elements
geoFacets = (
	(1,2,3),
	(0,2,3),
	(0,1,3),
	(0,1,2)
)
