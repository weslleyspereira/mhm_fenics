'''
Created on 3 de ago de 2017

@author: weslley
'''

from dolfin.cpp.mesh import Point, MeshFunction, SubDomain
from dolfin.cpp.function import near
from dolfin.cpp.io import File
from numpy import array

class Geometry():
    """Geometry class"""
    
    def __init__(self, points, faces):
        self._points = points # Points used for face elements
        self._faces = faces # Face elements
        if len(self._faces) > 0:
            if len(self._faces[0]) == 2:
                self._nDim = 2 # Topological dimension
            else:
                self._nDim = 3 # Topological dimension
        
    def getNoFaces(self):
        return len(self._faces)
    
    def getNormals(self):
        geoNormals = [Point() for i in range(self.getNoFaces())]
        if self._nDim == 2:
            for i in range(self.getNoFaces()):
                v = self._points[self._faces[i][1]]-self._points[self._faces[i][0]]
                geoNormals[i] = Point(v[1],-v[0])
        elif self._nDim == 3:
            geoNormals = [
                (self._points[face[1]]-self._points[face[0]]).cross(self._points[face[2]]-self._points[face[0]]) for face in self._faces
            ]
        else:
            print("Invalid dimension!")
            exit(1)
        return geoNormals

def markBoundariesOfMesh(mesh, geo, **kwargs):
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
        class BndSubDomain(SubDomain):
            def inside(self, x, on_boundary):
                return near( geoNormals[i].dot(Point(x)-geo._points[geoBndPoints[i]]), 0.0 ) and on_boundary
        aSubDomain = BndSubDomain()
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
        
    return faceSubdomains

def markBoundarySubdomains(mesh, boundaries, **kwargs):
    '''Mark the boundaries of a mesh given a geometry
    
    All interior faces are marked with zero
    
    boundaries is a dictionary {label: definition of boundary}
    '''
    
    bndSubdomains = MeshFunction("size_t", mesh, mesh.geometry().dim()-1)
    bndSubdomains.set_all(0)
    
    for i in boundaries:
        class BndSubDomain(SubDomain):
            def inside(self, x, on_boundary):
                value = array([0.0])
                boundaries[i].eval(value, x)
                return (value==1.0) and on_boundary
        aux = BndSubDomain()
        aux.mark(bndSubdomains, i)
        
    # Save Dolfin XML format of the subdomains
    if 'outputXmlGz' in kwargs:
        File(kwargs['outputXmlGz']) << bndSubdomains
        
    # Save sub domains to VTK files
    if 'outputPvd' in kwargs:
        File(kwargs['outputPvd']) << bndSubdomains
        
    return bndSubdomains
    
## If there exist dirichlet boundary
#if hasDirichlet:
#    if os.path.isfile(dirichletPath):
#        dirichletSubdomain = MeshFunction("bool", mesh, dirichletPath)
#    else:
#        # Mark faces
#        dirichletSubdomain.set_all(False)
#        class Aux(SubDomain):
#            def inside(self, x, on_boundary):
#                value = array([0.0])
#                physics.dirichletBnd.eval(value, x)
#                return (value==1.0) and on_boundary
#        aux = Aux()
#        aux.mark(dirichletSubdomain, True)
#        aux.mark(bndSubdomains, 0)