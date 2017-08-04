'''
Created on 3 de ago de 2017

@author: weslley
'''

from __future__ import print_function
from ufl.operators import grad, sym, tr
from ufl.constantvalue import Identity

###
### Elasticity definitions
###

# Stress tensor
def stress(mu,lmbda,u):
	Eps = sym(grad(u))
	return 2.0*mu*Eps + lmbda*tr(Eps)*Identity(len(u))
