from bright.apigen import cygen as cg

from nose.tools import assert_equal

toaster_desc = {
    'name': 'Toaster',
    'header_filename': 'toaster.h',
    'namespace': 'bright',
    'library_docstring': "I am the Toaster lib! Hear me sizzle!", 
    'parents': ['FCComp'],
    'attrs': {
        'nslices': 'uint',
        'toastiness': 'str',
        'rate': 'float',
        },
    'methods': {
        ('Toaster',): None,
        ('~Toaster',): None, 
        ('make_toast', ('when', 'str'), ('nslices', 'unit', 1)): 'int',
        },
    }


exp_cpppxd = \
"""################################################
#                 Warning!                     #
# This file has been auto-generated by Bright. #
# Do not modify!!!                             #
#                                              #
#                                              #
#                    Come on, guys. I mean it! #
################################################
from libcpp.string cimport string as std_string
from pyne cimport extra_types

cdef extern from "toaster.h" namespace "bright":

    cdef cppclass Toaster(FCComp):
        # constructors
        Toaster() except +
        ~Toaster()

        # attributes
        extra_types.uint nslices
        double rate
        std_string toastiness

        # methods
        int make_toast(std_string) except +
        int make_toast(std_string, extra_types.uint) except +


"""

def test_gencpppxd():
    obs = cg.gencpppxd(toaster_desc)
    assert_equal(obs, exp_cpppxd)
