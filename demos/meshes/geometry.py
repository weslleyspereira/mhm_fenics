'''
Created on 3 de ago de 2017

@author: weslley
'''

# from dolfin.cpp.mesh import Point

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