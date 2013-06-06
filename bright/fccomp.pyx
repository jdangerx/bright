################################################
#                 WARNING!                     #
# This file has been auto-generated by xdress. #
# Do not modify!!!                             #
#                                              #
#                                              #
#                    Come on, guys. I mean it! #
################################################
"""Python wrapper for fccomp.
"""
cimport pyne.stlcontainers
from libc.stdlib cimport free
from libcpp.map cimport map as cpp_map
from libcpp.set cimport set as cpp_set
from libcpp.string cimport string as std_string
from pyne cimport cpp_material
from pyne cimport material

from pyne import material
import pyne.stlcontainers

cdef class FCComp:
    """Abstract Base Fuel Cycle Component Class.
    
    Parameters
    ----------
    paramtrack : sequence of str, optional
        A set of parameter names (str) that the component will track.
    n : str, optional
        The name of the fuel cycle component instance.
        
    """

    # constuctors
    def __cinit__(self, *args, **kwargs):
        self._inst = NULL
        self._free_inst = True

        # cached property defaults
        self._mat_feed = None
        self._mat_prod = None
        self._params_after_calc = None
        self._params_prior_calc = None
        self._track_params = None

    def _fccomp_fccomp_0(self, n=""):
        """FCComp(self, n="")
        """
        cdef char * n_proxy
        n_bytes = n.encode()
        self._inst = new cpp_fccomp.FCComp(std_string(<char *> n_bytes))
    
    
    def _fccomp_fccomp_1(self, ptrack, n=""):
        """FCComp(self, ptrack, n="")
        """
        cdef pyne.stlcontainers._SetStr ptrack_proxy
        cdef char * n_proxy
        ptrack_proxy = pyne.stlcontainers.SetStr(ptrack, not isinstance(ptrack, pyne.stlcontainers._SetStr))
        n_bytes = n.encode()
        self._inst = new cpp_fccomp.FCComp(ptrack_proxy.set_ptr[0], std_string(<char *> n_bytes))
    
    
    _fccomp_fccomp_0_argtypes = frozenset(((0, str), ("n", str)))
    _fccomp_fccomp_1_argtypes = frozenset(((0, pyne.stlcontainers.SetStr), (1, str), ("ptrack", pyne.stlcontainers.SetStr), ("n", str)))
    
    def __init__(self, *args, **kwargs):
        """FCComp(self, ptrack, n="")
        """
        types = set([(i, type(a)) for i, a in enumerate(args)])
        types.update([(k, type(v)) for k, v in kwargs.items()])
        # vtable-like dispatch for exactly matching types
        if types <= self._fccomp_fccomp_0_argtypes:
            self._fccomp_fccomp_0(*args, **kwargs)
            return
        if types <= self._fccomp_fccomp_1_argtypes:
            self._fccomp_fccomp_1(*args, **kwargs)
            return
        # duck-typed dispatch based on whatever works!
        try:
            self._fccomp_fccomp_0(*args, **kwargs)
            return
        except (RuntimeError, TypeError, NameError):
            pass
        try:
            self._fccomp_fccomp_1(*args, **kwargs)
            return
        except (RuntimeError, TypeError, NameError):
            pass
        raise RuntimeError('method __init__() could not be dispatched')
    
    def __dealloc__(self):
        if self._free_inst:
            free(self._inst)

    # attributes
    property mat_feed:
        """A pyne.material.Material object that represents the flow of material into 
        this component for this pass."""
        def __get__(self):
            cdef material._Material mat_feed_proxy
            if self._mat_feed is None:
                mat_feed_proxy = material.Material(free_mat=False)
                mat_feed_proxy.mat_pointer = &(<cpp_fccomp.FCComp *> self._inst).mat_feed
                self._mat_feed = mat_feed_proxy
            return self._mat_feed
    
        def __set__(self, value):
            cdef material._Material value_proxy
            value_proxy = material.Material(value, free_mat=not isinstance(value, material._Material))
            (<cpp_fccomp.FCComp *> self._inst).mat_feed = value_proxy.mat_pointer[0]
            self._mat_feed = None
    
    
    property mat_prod:
        """A pyne.material.Material object that represents the flow of material out of
        this component for this pass.  Calling the FCComp.calc() method should calculate 
        mat_prod from the mat_feed value.
        """
        def __get__(self):
            cdef material._Material mat_prod_proxy
            if self._mat_prod is None:
                mat_prod_proxy = material.Material(free_mat=False)
                mat_prod_proxy.mat_pointer = &(<cpp_fccomp.FCComp *> self._inst).mat_prod
                self._mat_prod = mat_prod_proxy
            return self._mat_prod
    
        def __set__(self, value):
            cdef material._Material value_proxy
            value_proxy = material.Material(value, free_mat=not isinstance(value, material._Material))
            (<cpp_fccomp.FCComp *> self._inst).mat_prod = value_proxy.mat_pointer[0]
            self._mat_prod = None
    
    
    property name:
        """The string identifier for the component.  Defaults to an empty string."""
        def __get__(self):
            return bytes(<char *> (<cpp_fccomp.FCComp *> self._inst).name.c_str()).decode()
    
        def __set__(self, value):
            cdef char * value_proxy
            value_bytes = value.encode()
            (<cpp_fccomp.FCComp *> self._inst).name = std_string(<char *> value_bytes)
    
    
    property natural_name:
        """The natural name string identifier for the component. Computed from name value."""
        def __get__(self):
            return bytes(<char *> (<cpp_fccomp.FCComp *> self._inst).natural_name.c_str()).decode()
    
        def __set__(self, value):
            cdef char * value_proxy
            value_bytes = value.encode()
            (<cpp_fccomp.FCComp *> self._inst).natural_name = std_string(<char *> value_bytes)
    
    
    property params_after_calc:
        """A dictionary (or C++ map) that represents component-specific parameters at output 
        for this pass. The keys are restricted to strings while their associated values are 
        floats (doubles).  For example, a reactor component may have a "BUd" key that 
        represents the discharge burnup that the input fuel achieved.  This attribute does 
        not have a meaningful value until the calc_params() method is called."""
        def __get__(self):
            cdef pyne.stlcontainers._MapStrDouble params_after_calc_proxy
            if self._params_after_calc is None:
                params_after_calc_proxy = pyne.stlcontainers.MapStrDouble(False, False)
                params_after_calc_proxy.map_ptr = &(<cpp_fccomp.FCComp *> self._inst).params_after_calc
                self._params_after_calc = params_after_calc_proxy
            return self._params_after_calc
    
        def __set__(self, value):
            cdef pyne.stlcontainers._MapStrDouble value_proxy
            value_proxy = pyne.stlcontainers.MapStrDouble(value, not isinstance(value, pyne.stlcontainers._MapStrDouble))
            (<cpp_fccomp.FCComp *> self._inst).params_after_calc = value_proxy.map_ptr[0]
            self._params_after_calc = None
    
    
    property params_prior_calc:
        """A dictionary (or C++ map) that represents component-specific parameters at input 
        for this pass. The keys are restricted to strings while their associated values are 
        floats (doubles).  For example, a reprocessing component may record the total mass 
        input using a "Mass" key.  This attribute does not have a meaningful value until the 
        calc_params() method is called."""
        def __get__(self):
            cdef pyne.stlcontainers._MapStrDouble params_prior_calc_proxy
            if self._params_prior_calc is None:
                params_prior_calc_proxy = pyne.stlcontainers.MapStrDouble(False, False)
                params_prior_calc_proxy.map_ptr = &(<cpp_fccomp.FCComp *> self._inst).params_prior_calc
                self._params_prior_calc = params_prior_calc_proxy
            return self._params_prior_calc
    
        def __set__(self, value):
            cdef pyne.stlcontainers._MapStrDouble value_proxy
            value_proxy = pyne.stlcontainers.MapStrDouble(value, not isinstance(value, pyne.stlcontainers._MapStrDouble))
            (<cpp_fccomp.FCComp *> self._inst).params_prior_calc = value_proxy.map_ptr[0]
            self._params_prior_calc = None
    
    
    property pass_num:
        """An integer representing the number of passes this component has been cycled 
        through. It starts at zero and is incremented by one each time the write_mat_pass() 
        method is called."""
        def __get__(self):
            return int((<cpp_fccomp.FCComp *> self._inst).pass_num)
    
        def __set__(self, value):
            (<cpp_fccomp.FCComp *> self._inst).pass_num = <int> value
    
    
    property track_params:
        """A set of strings that holds the keys of params_prior_calc and params_after_calc.
        Every component type has its own set of parameters it is able to track.  This is why
        track_params is a component-specific attribute, while track_nucs is a bright-level
        variable."""
        def __get__(self):
            cdef pyne.stlcontainers._SetStr track_params_proxy
            if self._track_params is None:
                track_params_proxy = pyne.stlcontainers.SetStr(False, False)
                track_params_proxy.set_ptr = &(<cpp_fccomp.FCComp *> self._inst).track_params
                self._track_params = track_params_proxy
            return self._track_params
    
        def __set__(self, value):
            cdef pyne.stlcontainers._SetStr value_proxy
            value_proxy = pyne.stlcontainers.SetStr(value, not isinstance(value, pyne.stlcontainers._SetStr))
            (<cpp_fccomp.FCComp *> self._inst).track_params = value_proxy.set_ptr[0]
            self._track_params = None
    
    
    # methods
    def _fccomp_calc_0(self):
        """calc(self)
        This method is used to determine a component's output material from its input 
        material.  Therefore, this is typically where the bulk of a fuel cycle component's 
        algorithm lies.  As each component type has a distinct methodology, the calc() 
        method needs to be overridden child classes.
        
        This method should return mat_prod so that component calculations may be easily 
        daisy-chained together.
        
        Parameters
        ----------
        input : dict or Material, optional
            If input is present, it set as mat_feed.  If input is a nuclide dictionary 
            (zzaaam keys, float values), this dictionary is first converted to a Material 
            before being set as mat_feed.
        
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
    
    
    def _fccomp_calc_1(self, incomp):
        """calc(self, incomp)
        This method is used to determine a component's output material from its input 
        material.  Therefore, this is typically where the bulk of a fuel cycle component's 
        algorithm lies.  As each component type has a distinct methodology, the calc() 
        method needs to be overridden child classes.
        
        This method should return mat_prod so that component calculations may be easily 
        daisy-chained together.
        
        Parameters
        ----------
        input : dict or Material, optional
            If input is present, it set as mat_feed.  If input is a nuclide dictionary 
            (zzaaam keys, float values), this dictionary is first converted to a Material 
            before being set as mat_feed.
        
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
    
    
    def _fccomp_calc_2(self, mat):
        """calc(self, mat)
        This method is used to determine a component's output material from its input 
        material.  Therefore, this is typically where the bulk of a fuel cycle component's 
        algorithm lies.  As each component type has a distinct methodology, the calc() 
        method needs to be overridden child classes.
        
        This method should return mat_prod so that component calculations may be easily 
        daisy-chained together.
        
        Parameters
        ----------
        input : dict or Material, optional
            If input is present, it set as mat_feed.  If input is a nuclide dictionary 
            (zzaaam keys, float values), this dictionary is first converted to a Material 
            before being set as mat_feed.
        
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
    
    
    _fccomp_calc_0_argtypes = frozenset()
    _fccomp_calc_1_argtypes = frozenset(((0, pyne.stlcontainers.MapIntDouble), ("incomp", pyne.stlcontainers.MapIntDouble)))
    _fccomp_calc_2_argtypes = frozenset(((0, material.Material), ("mat", material.Material)))
    
    def calc(self, *args, **kwargs):
        """calc(self, mat)
        This method is used to determine a component's output material from its input 
        material.  Therefore, this is typically where the bulk of a fuel cycle component's 
        algorithm lies.  As each component type has a distinct methodology, the calc() 
        method needs to be overridden child classes.
        
        This method should return mat_prod so that component calculations may be easily 
        daisy-chained together.
        
        Parameters
        ----------
        input : dict or Material, optional
            If input is present, it set as mat_feed.  If input is a nuclide dictionary 
            (zzaaam keys, float values), this dictionary is first converted to a Material 
            before being set as mat_feed.
        
        Returns
        -------
        output : Material
            mat_prod
        
        """
        types = set([(i, type(a)) for i, a in enumerate(args)])
        types.update([(k, type(v)) for k, v in kwargs.items()])
        # vtable-like dispatch for exactly matching types
        if types <= self._fccomp_calc_0_argtypes:
            return self._fccomp_calc_0(*args, **kwargs)
        if types <= self._fccomp_calc_1_argtypes:
            return self._fccomp_calc_1(*args, **kwargs)
        if types <= self._fccomp_calc_2_argtypes:
            return self._fccomp_calc_2(*args, **kwargs)
        # duck-typed dispatch based on whatever works!
        try:
            return self._fccomp_calc_0(*args, **kwargs)
        except (RuntimeError, TypeError, NameError):
            pass
        try:
            return self._fccomp_calc_1(*args, **kwargs)
        except (RuntimeError, TypeError, NameError):
            pass
        try:
            return self._fccomp_calc_2(*args, **kwargs)
        except (RuntimeError, TypeError, NameError):
            pass
        raise RuntimeError('method calc() could not be dispatched')
    
    def calc_params(self):
        """calc_params(self)
        By calling this method, all parameter values are calculated and set for the fuel 
        cycle component.  This should be done following a calc() calculation but before data 
        is written out.  If a component has important parameters associated with it, this 
        function must be overridden and called.
        
        Note that this is called first thing when write_params_pass() is called.  For 
        example, reprocessing only has a "Mass" parameter.  Translated into Python, 
        calc_params() here looks like the following::
        
            def calc_params(self):
                self.params_prior_calc["Mass"] = self.mat_feed.mass
                self.params_after_calc["Mass"] = self.mat_prod.mass
                return
        
        """
        (<cpp_fccomp.FCComp *> self._inst).calc_params()
    
    
    def write(self):
        """write(self)
        This is a convenience function that first increments up pass_num.
        Then, it checks to see if there are any parameters for this component.
        If there are, it sets the current values using calc_params().
        
        If bright.bright_conf.write_hdf5 is True, then write_hdf5() is called.
        
        If bright.bright_conf.write_text is True, then write_text() is called.
        
        This is what is most often used to write Bright output.  Therefore it is
        typically the last step for every component in each pass.
        """
        (<cpp_fccomp.FCComp *> self._inst).write()
    
    
    def write_hdf5(self):
        """write_hdf5(self)
        This method writes out the isotopic pass data to an HDF5 file. 
        Then, if available, it also writes parameter data as well.  
        Using write() instead is recommended.
        """
        (<cpp_fccomp.FCComp *> self._inst).write_hdf5()
    
    
    def write_mat_pass(self):
        """write_mat_pass(self)
        This method is responsible for adding the current material data for this pass to 
        the output text and hdf5 files for this component.  Further calculations for this 
        pass should not be performed after write_mat_pass() has been called.
        
        This function has one very important subtletywhen writing *text* output. It does 
        not write out material data directly. Rather, input columns are given as normalized 
        isotopic vectors.  As weight fractions, input columns are in units of 
        [kg_mat_feed[nuc]/kg_mat_feed.mass].  Moreover, the output columns are given in 
        terms relative to the mass of the input mass, [kg_mat_prod[nuc]/kg_mat_feed.mass].  
        These are calculated via the following expressions.
        
        .. math::
        
            \\mbox{inpcol[nuc]} = \\mbox{mat\_feed.comp[nuc]}
        
            \\mbox{outcol[nuc]} = \\mbox{mat\_prod.comp[nuc]} \\times \\frac{\\mbox{mat\_prod.mass}}{\\mbox{mat\_feed.mass}}
        
        Because of the units of these two columns, total mass flow data may often only be 
        recovered via the "Mass" parameter in the parameter file.  Here is a sample 
        LWRIsos.txt file for a light water reactor for the first pass::
        
            Isotope 1in             1out    
            H1      0.000000E+00    0.000000E+00
            H3      0.000000E+00    8.568522E-08
            HE4     0.000000E+00    4.421615E-07
            B10     0.000000E+00    0.000000E+00
            B11     0.000000E+00    0.000000E+00
            C14     0.000000E+00    4.015091E-11
            O16     0.000000E+00    0.000000E+00
            SR90    0.000000E+00    8.221283E-04
            TC99    0.000000E+00    1.112580E-03
            CS137   0.000000E+00    1.821226E-03
            U234    0.000000E+00    2.807466E-06
            U235    4.773292E-02    8.951725E-03
            U236    0.000000E+00    6.155297E-03
            U237    0.000000E+00    1.719458E-05
            U238    9.522671E-01    9.211956E-01
            U239    0.000000E+00    6.953862E-07
            NP237   0.000000E+00    8.057270E-04
            PU238   0.000000E+00    2.842232E-04
            PU239   0.000000E+00    5.353362E-03
            PU240   0.000000E+00    2.114728E-03
        
        Note that HDF5 ouput perisists materials themseleves and do not present this 
        normalization problem.
        """
        (<cpp_fccomp.FCComp *> self._inst).write_mat_pass()
    
    
    def write_params_pass(self):
        """write_params_pass(self)
        What write_mat_pass() does for a component's input and output isotopics, 
        this method does for the components parameters.  To ensure that meaningful 
        data is available, calc_params() should be called elsewhere in the program  
        prior to this method.  Note that to get the pass numbering correct, 
        pass_num should always be incremented prior to this method.  The 
        following is an example of a text file output for a light water 
        reactor spent fuel reprocessing facility::
        
            Param   1in             1out    
            Mass    9.985828E-01    9.975915E-01
        
        """
        (<cpp_fccomp.FCComp *> self._inst).write_params_pass()
    
    
    def write_text(self):
        """write_text(self)
        This method calls write_mat_pass() and then, if available, calls 
        write_params_pass().  This is convience function for producing 
        text-based output.  However, using write() is recommended.
        """
        (<cpp_fccomp.FCComp *> self._inst).write_text()
    
    





