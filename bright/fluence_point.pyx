################################################
#                 WARNING!                     #
# This file has been auto-generated by Bright. #
# Do not modify!!!                             #
#                                              #
#                                              #
#                    Come on, guys. I mean it! #
################################################
"""Python wrapper for the fluence point.
"""




cdef class FluencePoint:
    """This class holds three simple data points that represent a fluence point."""

    # constuctors
    def __cinit__(self, *args, **kwargs):
        self._inst = NULL
        self._free_inst = True

        # cached property defaults


    def FluencePoint(self):
        """"""
        self._inst = new cpp_fluence_point.FluencePoint()
    
    

    # attributes
    property F:
        """Fluence value itself (float).  In units of [neutrons/kilobarn], abbr [n/kb]."""
        def __get__(self):
            return float((<cpp_fluence_point.FluencePoint *> self._inst).F)
    
        def __set__(self, value):
            (<cpp_fluence_point.FluencePoint *> self._inst).F = <double> value
    
    
    property f:
        """Index (int) of fluence immediately lower than the value of F."""
        def __get__(self):
            return int((<cpp_fluence_point.FluencePoint *> self._inst).f)
    
        def __set__(self, value):
            (<cpp_fluence_point.FluencePoint *> self._inst).f = value
    
    
    property m:
        """The slope (float) dBU/dF between points f and f+1.  Has the odd units of [MWd kb / kgIHM n]"""
        def __get__(self):
            return float((<cpp_fluence_point.FluencePoint *> self._inst).m)
    
        def __set__(self, value):
            (<cpp_fluence_point.FluencePoint *> self._inst).m = <double> value
    
    
    # methods

