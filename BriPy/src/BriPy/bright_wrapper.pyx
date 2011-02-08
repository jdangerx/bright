"""Python wrapper for isoname library."""
# Cython imports
from libcpp.map cimport map as cpp_map
from libcpp.set cimport set as cpp_set
from cython cimport pointer
from cython.operator cimport dereference as deref
from cython.operator cimport preincrement as inc

# local imports 
cimport std
cimport cpp_bright

import isoname


#
# FCComps Configuration namespace
#

class isos2track(object):
    def __get_value__(self):
        value = set()
        cdef cpp_set[int].iterator isos_iter = cpp_bright.isos2track.begin()
        while isos_iter != cpp_bright.isos2track.end():
            value.add(deref(isos_iter))
            inc(isos_iter)
        return value

    def __set_value__(self, value):
        cdef int iso_zz
        cpp_bright.isos2track.clear()
        for iso in value:
            iso_zz = isoname.mixed_2_zzaaam(iso)
            cpp_bright.isos2track.insert(iso_zz)

    value = property(__get_value__, __set_value__)

# Make isos2track a singleton
isos2track = isos2track().value

# Load isos2track from file functions
def load_isos2track_hdf5(char * filename, char * datasetname="", bint clear=False):
    """This convience function tries to load the isos2track set from a dataset 
    in an HDF5 file.  The dataset *must* be of integer type.  String-based
    nuclide names are currently not supported. 

    Args:
        * filename (str): Path to the data library.
        * dataset (str):  Dataset name to grab nuclides from.
        * clear (bool):   Flag that if set removes the currrent entries
          from isos2track prior to loading in new values.

    If the dataset argument is not provided or empty, the function tries to 
    load from various default datasets in the following order::

        "/isos2track"  
        "/Isos2Track"
        "/isostrack"   
        "/IsosTrack"
        "/isotrack"   
        "/IsoTrack"    
        "/ToIso"
        "/ToIsos"
        "/ToIso_zz"
        "/ToIso_MCNP"
        "/FromIso"  
        "/FromIsos"  
        "/FromIso_zz" 
        "/FromIso_MCNP"
    """
    cpp_bright.load_isos2track_hdf5(std.string(filename), std.string(datasetname), clear)


def load_isos2track_text(char * filename, bint clear=False):
    """This convience function tries to load the isos2track set from a text
    file.  The nuclide names may use any naming convention.  Mixing different
    conventions in the same file is allowed.  Whitespace is required between
    isotopic names.

    Args:
        * filename (str): Path to the data library.
        * clear (bool):   Flag that if set removes the currrent entries
          from isos2track prior to loading in new values.
    """
    cpp_bright.load_isos2track_text(std.string(filename), clear)


# Simple settings
verbosity = cpp_bright.verbosity
write_hdf5 = cpp_bright.write_hdf5
write_text = cpp_bright.write_text

class output_filename(object):
    def __get_value__(self):
        cdef std.string value = cpp_bright.output_filename
        return value.c_str()

    def __set_value__(self, char * value):
        cpp_bright.output_filename = std.string(value)

    value = property(__get_value__, __set_value__)

# Make isos2track a singleton
output_filename = output_filename().value

