################################################
#                 WARNING!                     #
# This file has been auto-generated by xdress. #
# Do not modify!!!                             #
#                                              #
#                                              #
#                    Come on, guys. I mean it! #
################################################




cdef extern from "fluence_point.h" namespace "bright":

    cdef cppclass FluencePoint:
        # constructors
        FluencePoint() except +
        FluencePoint(int) except +
        FluencePoint(int, double) except +
        FluencePoint(int, double, double) except +

        # attributes
        double F
        int f
        double m

        # methods

        pass




