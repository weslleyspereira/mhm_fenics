'''
Created on 3 de ago de 2017

@author: weslley
'''
import unittest
from demos.meshes.geometry import Geometry, markBoundariesOfMesh
from dolfin.cpp.mesh import Point, Mesh, MeshFunction


class Test(unittest.TestCase):


    def testName(self):
        pass
    
    
    def testMarkSquare(self):
        
        geo = Geometry(
    
            # Points used for face elements
            (
                Point(0.0,0.0),
                Point(1.0,0.0),
                Point(1.0,1.0),
                Point(0.0,1.0)
            ),
                    
            # Face elements
            (
                (0,1),
                (1,2),
                (2,3),
                (3,0)
            )
                    )
        
        mesh = Mesh("triangle.7.nml.gz")
        
        faceSubdomains = markBoundariesOfMesh(mesh, geo)
        
        self.assertTrue(isinstance(faceSubdomains, MeshFunction), 
                        "faceSubdomains is not a MeshFunction")

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()