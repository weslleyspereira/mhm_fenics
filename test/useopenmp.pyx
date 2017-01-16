from cython.parallel cimport parallel
cimport openmp

def getnumthreads():
    cdef int num_threads
    openmp.omp_set_dynamic(0)
    openmp.omp_set_num_threads(8)
    with nogil, parallel():
        num_threads = openmp.omp_get_num_threads()
        with gil:
            print openmp.omp_get_thread_num()
            return num_threads
