################################################
#                 WARNING!                     #
# This file has been auto-generated by xdress. #
# Do not modify!!!                             #
#                                              #
#                                              #
#                    Come on, guys. I mean it! #
################################################


cimport fccomp
cimport pyne.stlcontainers
cimport reactor1g
from bright cimport cpp_fuel_fabrication
from bright cimport cpp_reactor1g
from libcpp.map cimport map as cpp_map
from libcpp.string cimport string as std_string
from pyne cimport cpp_material
from pyne cimport material



cdef class FuelFabrication(fccomp.FCComp):
    cdef public pyne.stlcontainers._MapStrDouble _deltaRs
    cdef public pyne.stlcontainers._MapStrDouble _mass_weights_in
    cdef public pyne.stlcontainers._MapStrDouble _mass_weights_out
    cdef public material._MapStrMaterial _materials
    cdef public reactor1g.Reactor1G _reactor
    pass    




