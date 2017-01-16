from cython.parallel cimport parallel
cimport openmp

def sumsquare(n):
	with nogil:
		"This works even if the code does access the GIL. It doesn't release it."
	return sum( i*i for i in range(n))
