################################################
#                 WARNING!                     #
# This file has been auto-generated by xdress. #
# Do not modify!!!                             #
#                                              #
#                                              #
#                    Come on, guys. I mean it! #
################################################


cimport numpy as np
cimport pyne.stlcontainers
from bright cimport cpp_reactor_parameters
from libcpp.map cimport map as cpp_map
from libcpp.string cimport string as std_string
from libcpp.vector cimport vector as cpp_vector

cdef class ReactorParameters:
    cdef void * _inst
    cdef public bint _free_inst
    cdef public np.ndarray _burn_times
    cdef public pyne.stlcontainers._MapStrDouble _cladding_form
    cdef public pyne.stlcontainers._MapStrDouble _coolant_form
    cdef public pyne.stlcontainers._MapStrDouble _fuel_form




