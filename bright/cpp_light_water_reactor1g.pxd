################################################
#                 WARNING!                     #
# This file has been auto-generated by xdress. #
# Do not modify!!!                             #
#                                              #
#                                              #
#                    Come on, guys. I mean it! #
################################################


from bright cimport cpp_reactor1g
from bright cimport cpp_reactor_parameters
from libcpp.string cimport string as std_string

cdef extern from "light_water_reactor1g.h" namespace "bright":

    cdef cppclass LightWaterReactor1G(cpp_reactor1g.Reactor1G):
        # constructors
        LightWaterReactor1G() except +
        LightWaterReactor1G(std_string) except +
        LightWaterReactor1G(std_string, std_string) except +
        LightWaterReactor1G(std_string, cpp_reactor_parameters.ReactorParameters) except +
        LightWaterReactor1G(std_string, cpp_reactor_parameters.ReactorParameters, std_string) except +
        LightWaterReactor1G(cpp_reactor_parameters.ReactorParameters) except +
        LightWaterReactor1G(cpp_reactor_parameters.ReactorParameters, std_string) except +

        # attributes


        # methods
        void calc_params() except +




