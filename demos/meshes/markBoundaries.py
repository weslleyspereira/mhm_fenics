'''
Created on 3 de ago de 2017

@author: weslley
'''

from dolfin.cpp.mesh import Point, MeshFunction, SubDomain
from dolfin.cpp.function import near
from dolfin.cpp.io import File

def markBoundariesOfMesh(mesh,geo,**kwargs):
	'''Mark the boundaries of a mesh given a geometry
	
	All interior faces are marked with zero
	
	If the geometry has no marks, it will be mared using the convention:
	face[0] = 1
	face[1] = 2
	...
	face[N] = N+1
	'''
	
	# Reference point indexes in each face
	geoBndPoints = [ face[0] for face in geo._faces ]

	# Face normals
	geoNormals = geo.getNormals()

	# Create mesh functions over the cell facets
	faceSubdomains = MeshFunction("size_t", mesh, geo._nDim-1)

	# Mark all facets of the domain as 0
	faceSubdomains.set_all(0)

	# Mark boundary faces for face elements
	for i in range(geo.getNoFaces()):
		class SubDomain(SubDomain):
			def inside(self, x, on_boundary):
				return near( geoNormals[i].dot(Point(x)-geo._points[geoBndPoints[i]]), 0.0 ) and on_boundary
		aSubDomain = SubDomain()
		label = i+1
		if hasattr(geo, 'facesMarkers'):
			label = geo.facesMarkers[i]
		aSubDomain.mark(faceSubdomains, label)
		
	# Save Dolfin XML format of the subdomains
	if 'outputXmlGz' in kwargs:
		File(kwargs['outputXmlGz']) << faceSubdomains
		
	# Save sub domains to VTK files
	if 'outputPvd' in kwargs:
		File(kwargs['outputPvd']) << faceSubdomains