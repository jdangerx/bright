################################################
#                 WARNING!                     #
# This file has been auto-generated by xdress. #
# Do not modify!!!                             #
#                                              #
#                                              #
#                    Come on, guys. I mean it! #
################################################
"""Python wrapper for Storage.
"""
cimport fccomp
cimport pyne.stlcontainers
from bright cimport cpp_fccomp
from libc.stdlib cimport free
from libcpp.map cimport map as cpp_map
from libcpp.string cimport string as std_string
from pyne cimport cpp_material
from pyne cimport material

from pyne import material
import fccomp
import pyne.stlcontainers



cdef class decay_nuc:
    """no docstring for decay_nuc, please file a bug report!"""



    # constuctors
    def __cinit__(self, *args, **kwargs):
        self._inst = NULL
        self._free_inst = True

        # cached property defaults


    def __init__(self):
        """decay_nuc(self)
        """
        self._inst = new cpp_storage.decay_nuc()
    
    
    def __dealloc__(self):
        if self._free_inst:
            free(self._inst)

    # attributes
    property branchratio:
        """no docstring for branchratio, please file a bug report!"""
        def __get__(self):
            return float((<cpp_storage.decay_nuc *> self._inst).branchratio)
    
        def __set__(self, value):
            (<cpp_storage.decay_nuc *> self._inst).branchratio = <double> value
    
    
    property decayconst:
        """no docstring for decayconst, please file a bug report!"""
        def __get__(self):
            return float((<cpp_storage.decay_nuc *> self._inst).decayconst)
    
        def __set__(self, value):
            (<cpp_storage.decay_nuc *> self._inst).decayconst = <double> value
    
    
    property fromiso:
        """no docstring for fromiso, please file a bug report!"""
        def __get__(self):
            return int((<cpp_storage.decay_nuc *> self._inst).fromiso)
    
        def __set__(self, value):
            (<cpp_storage.decay_nuc *> self._inst).fromiso = <int> value
    
    
    property halflife:
        """no docstring for halflife, please file a bug report!"""
        def __get__(self):
            return float((<cpp_storage.decay_nuc *> self._inst).halflife)
    
        def __set__(self, value):
            (<cpp_storage.decay_nuc *> self._inst).halflife = <double> value
    
    
    property toiso:
        """no docstring for toiso, please file a bug report!"""
        def __get__(self):
            return int((<cpp_storage.decay_nuc *> self._inst).toiso)
    
        def __set__(self, value):
            (<cpp_storage.decay_nuc *> self._inst).toiso = <int> value
    
    
    # methods


    pass





cdef class Storage(fccomp.FCComp):
    """Storage Fuel Cycle Component Class.  Daughter of FCComp class.
    
    Parameters
    ----------
    n : str, optional
        The name of the storage fuel cycle component instance.
    
    """



    # constuctors
    def __cinit__(self, *args, **kwargs):
        self._inst = NULL
        self._free_inst = True

        # cached property defaults


    def __init__(self, n=""):
        """Storage(self, n="")
        """
        cdef char * n_proxy
        n_bytes = n.encode()
        self._inst = new cpp_storage.Storage(std_string(<char *> n_bytes))
    
    

    # attributes
    property decay_time:
        """This the float (double) attribute that represents how long an input fuel 
        mass should be stored for.  This time is represented in seconds, so be sure 
        to convert to the proper units before using.  Consider using pyne.utils.to_sec() 
        for second conversions.
        """
        def __get__(self):
            return float((<cpp_storage.Storage *> self._inst).decay_time)
    
        def __set__(self, value):
            (<cpp_storage.Storage *> self._inst).decay_time = <double> value
    
    
    # methods
    def _storage_calc_0(self):
        """calc(self)
        As usual, calc sets up the Storage component's input stream and calculates 
        the corresponding output Material.  Here, this amounts to calling bateman() 
        for every nuclide in mat_feed, for each chain that ends with a nuclide in track_nucs.
        
        Parameters
        ----------
        input : dict or Material or None, optional 
            If input is present, it set as the component's mat_feed.  If input is a 
            isotopic dictionary (zzaaam keys, float values), this dictionary is first 
            converted into a Material before being set as mat_feed.
        decay_time : float or None, optional 
            decay_time is set to the time value here prior to any other calculations.  
            This time has units of seconds.
        
        Returns
        -------
        output : Material
            mat_prod
        
        """
        cdef cpp_material.Material rtnval
        cdef material._Material rtnval_proxy
        rtnval = (<cpp_fccomp.FCComp *> self._inst).calc()
        rtnval_proxy = material.Material()
        rtnval_proxy.mat_pointer[0] = rtnval
        return rtnval_proxy
    
    
    def _storage_calc_1(self, incomp):
        """calc(self, incomp)
        As usual, calc sets up the Storage component's input stream and calculates 
        the corresponding output Material.  Here, this amounts to calling bateman() 
        for every nuclide in mat_feed, for each chain that ends with a nuclide in track_nucs.
        
        Parameters
        ----------
        input : dict or Material or None, optional 
            If input is present, it set as the component's mat_feed.  If input is a 
            isotopic dictionary (zzaaam keys, float values), this dictionary is first 
            converted into a Material before being set as mat_feed.
        decay_time : float or None, optional 
            decay_time is set to the time value here prior to any other calculations.  
            This time has units of seconds.
        
        Returns
        -------
        output : Material
            mat_prod
        
        """
        cdef pyne.stlcontainers._MapIntDouble incomp_proxy
        cdef cpp_material.Material rtnval
        cdef material._Material rtnval_proxy
        incomp_proxy = pyne.stlcontainers.MapIntDouble(incomp, not isinstance(incomp, pyne.stlcontainers._MapIntDouble))
        rtnval = (<cpp_fccomp.FCComp *> self._inst).calc(incomp_proxy.map_ptr[0])
        rtnval_proxy = material.Material()
        rtnval_proxy.mat_pointer[0] = rtnval
        return rtnval_proxy
    
    
    def _storage_calc_2(self, incomp, t):
        """calc(self, incomp, t)
        As usual, calc sets up the Storage component's input stream and calculates 
        the corresponding output Material.  Here, this amounts to calling bateman() 
        for every nuclide in mat_feed, for each chain that ends with a nuclide in track_nucs.
        
        Parameters
        ----------
        input : dict or Material or None, optional 
            If input is present, it set as the component's mat_feed.  If input is a 
            isotopic dictionary (zzaaam keys, float values), this dictionary is first 
            converted into a Material before being set as mat_feed.
        decay_time : float or None, optional 
            decay_time is set to the time value here prior to any other calculations.  
            This time has units of seconds.
        
        Returns
        -------
        output : Material
            mat_prod
        
        """
        cdef pyne.stlcontainers._MapIntDouble incomp_proxy
        cdef cpp_material.Material rtnval
        cdef material._Material rtnval_proxy
        incomp_proxy = pyne.stlcontainers.MapIntDouble(incomp, not isinstance(incomp, pyne.stlcontainers._MapIntDouble))
        rtnval = (<cpp_storage.Storage *> self._inst).calc(incomp_proxy.map_ptr[0], <double> t)
        rtnval_proxy = material.Material()
        rtnval_proxy.mat_pointer[0] = rtnval
        return rtnval_proxy
    
    
    def _storage_calc_3(self, mat):
        """calc(self, mat)
        As usual, calc sets up the Storage component's input stream and calculates 
        the corresponding output Material.  Here, this amounts to calling bateman() 
        for every nuclide in mat_feed, for each chain that ends with a nuclide in track_nucs.
        
        Parameters
        ----------
        input : dict or Material or None, optional 
            If input is present, it set as the component's mat_feed.  If input is a 
            isotopic dictionary (zzaaam keys, float values), this dictionary is first 
            converted into a Material before being set as mat_feed.
        decay_time : float or None, optional 
            decay_time is set to the time value here prior to any other calculations.  
            This time has units of seconds.
        
        Returns
        -------
        output : Material
            mat_prod
        
        """
        cdef material._Material mat_proxy
        cdef cpp_material.Material rtnval
        cdef material._Material rtnval_proxy
        mat_proxy = material.Material(mat, free_mat=not isinstance(mat, material._Material))
        rtnval = (<cpp_fccomp.FCComp *> self._inst).calc(mat_proxy.mat_pointer[0])
        rtnval_proxy = material.Material()
        rtnval_proxy.mat_pointer[0] = rtnval
        return rtnval_proxy
    
    
    def _storage_calc_4(self, mat, t):
        """calc(self, mat, t)
        As usual, calc sets up the Storage component's input stream and calculates 
        the corresponding output Material.  Here, this amounts to calling bateman() 
        for every nuclide in mat_feed, for each chain that ends with a nuclide in track_nucs.
        
        Parameters
        ----------
        input : dict or Material or None, optional 
            If input is present, it set as the component's mat_feed.  If input is a 
            isotopic dictionary (zzaaam keys, float values), this dictionary is first 
            converted into a Material before being set as mat_feed.
        decay_time : float or None, optional 
            decay_time is set to the time value here prior to any other calculations.  
            This time has units of seconds.
        
        Returns
        -------
        output : Material
            mat_prod
        
        """
        cdef material._Material mat_proxy
        cdef cpp_material.Material rtnval
        cdef material._Material rtnval_proxy
        mat_proxy = material.Material(mat, free_mat=not isinstance(mat, material._Material))
        rtnval = (<cpp_storage.Storage *> self._inst).calc(mat_proxy.mat_pointer[0], <double> t)
        rtnval_proxy = material.Material()
        rtnval_proxy.mat_pointer[0] = rtnval
        return rtnval_proxy
    
    
    def _storage_calc_5(self, t=0.0):
        """calc(self, t=0.0)
        As usual, calc sets up the Storage component's input stream and calculates 
        the corresponding output Material.  Here, this amounts to calling bateman() 
        for every nuclide in mat_feed, for each chain that ends with a nuclide in track_nucs.
        
        Parameters
        ----------
        input : dict or Material or None, optional 
            If input is present, it set as the component's mat_feed.  If input is a 
            isotopic dictionary (zzaaam keys, float values), this dictionary is first 
            converted into a Material before being set as mat_feed.
        decay_time : float or None, optional 
            decay_time is set to the time value here prior to any other calculations.  
            This time has units of seconds.
        
        Returns
        -------
        output : Material
            mat_prod
        
        """
        cdef cpp_material.Material rtnval
        cdef material._Material rtnval_proxy
        rtnval = (<cpp_storage.Storage *> self._inst).calc(<double> t)
        rtnval_proxy = material.Material()
        rtnval_proxy.mat_pointer[0] = rtnval
        return rtnval_proxy
    
    
    def _storage_calc_6(self, t):
        """calc(self, t)
        As usual, calc sets up the Storage component's input stream and calculates 
        the corresponding output Material.  Here, this amounts to calling bateman() 
        for every nuclide in mat_feed, for each chain that ends with a nuclide in track_nucs.
        
        Parameters
        ----------
        input : dict or Material or None, optional 
            If input is present, it set as the component's mat_feed.  If input is a 
            isotopic dictionary (zzaaam keys, float values), this dictionary is first 
            converted into a Material before being set as mat_feed.
        decay_time : float or None, optional 
            decay_time is set to the time value here prior to any other calculations.  
            This time has units of seconds.
        
        Returns
        -------
        output : Material
            mat_prod
        
        """
        cdef cpp_material.Material rtnval
        cdef material._Material rtnval_proxy
        rtnval = (<cpp_storage.Storage *> self._inst).calc(<double> t)
        rtnval_proxy = material.Material()
        rtnval_proxy.mat_pointer[0] = rtnval
        return rtnval_proxy
    
    
    _storage_calc_0_argtypes = frozenset()
    _storage_calc_1_argtypes = frozenset(((0, pyne.stlcontainers.MapIntDouble), ("incomp", pyne.stlcontainers.MapIntDouble)))
    _storage_calc_2_argtypes = frozenset(((0, pyne.stlcontainers.MapIntDouble), (1, float), ("incomp", pyne.stlcontainers.MapIntDouble), ("t", float)))
    _storage_calc_3_argtypes = frozenset(((0, material.Material), ("mat", material.Material)))
    _storage_calc_4_argtypes = frozenset(((0, material.Material), (1, float), ("mat", material.Material), ("t", float)))
    _storage_calc_5_argtypes = frozenset(((0, float), ("t", float)))
    _storage_calc_6_argtypes = frozenset(((0, float), ("t", float)))
    
    def calc(self, *args, **kwargs):
        """calc(self, t)
        As usual, calc sets up the Storage component's input stream and calculates 
        the corresponding output Material.  Here, this amounts to calling bateman() 
        for every nuclide in mat_feed, for each chain that ends with a nuclide in track_nucs.
        
        Parameters
        ----------
        input : dict or Material or None, optional 
            If input is present, it set as the component's mat_feed.  If input is a 
            isotopic dictionary (zzaaam keys, float values), this dictionary is first 
            converted into a Material before being set as mat_feed.
        decay_time : float or None, optional 
            decay_time is set to the time value here prior to any other calculations.  
            This time has units of seconds.
        
        Returns
        -------
        output : Material
            mat_prod
        
        """
        types = set([(i, type(a)) for i, a in enumerate(args)])
        types.update([(k, type(v)) for k, v in kwargs.items()])
        # vtable-like dispatch for exactly matching types
        if types <= self._storage_calc_0_argtypes:
            return self._storage_calc_0(*args, **kwargs)
        if types <= self._storage_calc_1_argtypes:
            return self._storage_calc_1(*args, **kwargs)
        if types <= self._storage_calc_3_argtypes:
            return self._storage_calc_3(*args, **kwargs)
        if types <= self._storage_calc_5_argtypes:
            return self._storage_calc_5(*args, **kwargs)
        if types <= self._storage_calc_6_argtypes:
            return self._storage_calc_6(*args, **kwargs)
        if types <= self._storage_calc_2_argtypes:
            return self._storage_calc_2(*args, **kwargs)
        if types <= self._storage_calc_4_argtypes:
            return self._storage_calc_4(*args, **kwargs)
        # duck-typed dispatch based on whatever works!
        try:
            return self._storage_calc_0(*args, **kwargs)
        except (RuntimeError, TypeError, NameError):
            pass
        try:
            return self._storage_calc_1(*args, **kwargs)
        except (RuntimeError, TypeError, NameError):
            pass
        try:
            return self._storage_calc_3(*args, **kwargs)
        except (RuntimeError, TypeError, NameError):
            pass
        try:
            return self._storage_calc_5(*args, **kwargs)
        except (RuntimeError, TypeError, NameError):
            pass
        try:
            return self._storage_calc_6(*args, **kwargs)
        except (RuntimeError, TypeError, NameError):
            pass
        try:
            return self._storage_calc_2(*args, **kwargs)
        except (RuntimeError, TypeError, NameError):
            pass
        try:
            return self._storage_calc_4(*args, **kwargs)
        except (RuntimeError, TypeError, NameError):
            pass
        raise RuntimeError('method calc() could not be dispatched')
    
    def calc_params(self):
        """calc_params(self)
        Here the parameters for Storage are set.  For storage, this amounts to just
        a "Mass" parameter::
        
            self.params_prior_calc["Mass"]  = self.mat_feed.mass
            self.params_after_calc["Mass"] = self.mat_prod.mass
        
        """
        (<cpp_fccomp.FCComp *> self._inst).calc_params()
    
    

    pass





