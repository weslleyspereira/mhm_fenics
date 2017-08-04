'''
Created on 3 de ago de 2017

@author: weslley
'''

from dolfin.cpp.mesh import Point
from demos.meshes.geometry import Geometry

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

## Face markers
#facesMarkers = (
#	1,
#	2,
#	3,
#	4
#)
