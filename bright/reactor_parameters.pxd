################################################
#                 WARNING!                     #
# This file has been auto-generated by Bright. #
# Do not modify!!!                             #
#                                              #
#                                              #
#                    Come on, guys. I mean it! #
################################################
cimport numpy as np
from bright cimport cpp_reactor_parameters
from libcpp.map cimport map as cpp_map
from libcpp.string cimport string as std_string
from libcpp.vector cimport vector as cpp_vector
from pyne cimport stlconverters as conv

cdef class ReactorParameters:
    cdef void * _inst
    cdef public bint _free_inst
    cdef public np.ndarray _burn_times
    cdef public conv._MapStrDouble _cladding_form
    cdef public conv._MapStrDouble _coolant_form
    cdef public conv._MapStrDouble _fuel_form


