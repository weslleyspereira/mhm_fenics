'''
Created on 3 de ago de 2017

@author: weslley
'''

from dolfin.cpp.mesh import Point

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
        geoNormals = [Point() for i in range(self._faces)]
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