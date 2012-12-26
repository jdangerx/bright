from bright.apigen import typesystem as ts

from nose.tools import assert_equal, with_setup

# setup and teardown new refinement cases
new_refined = {
    'comp_map': ['map', 'nucid', 'float64'],
    ('intrange', ('low', 'int32'), ('high', 'int32')): 'int32',
    ('nucrange', ('low', 'nucid'), ('high', 'nucid')): 'nucid',
    ('range', 'vtype', ('low', 'vtype'), ('high', 'vtype')): 'vtype',
    }
add_new_refined = lambda: ts.refined_types.update(new_refined)
del_new_refined = lambda: [ts.refined_types.pop(key) for key in new_refined]



def check_canon(t, exp):
    obs = ts.canon(t)
    assert_equal(obs, exp)

@with_setup(add_new_refined, del_new_refined)
def test_canon():
    cases = (
        ('str', 'str'),
        (['str'], ('str', 0)),
        ('f4', 'float32'),
        ('nucid', ('int32', 'nucid')),
        (['nucid'], (('int32', 'nucid'), 0)),
        (['set', 'complex'], ('set', 'complex128', 0)),
        (['map', 'nucid', 'float'], ('map', ('int32', 'nucid'), 'float64', 0)),
        ('comp_map', (('map', ('int32', 'nucid'), 'float64', 0), 'comp_map')),
        (['char', '*'], ('char', '*')),
        (['char', 42], ('char', 42)),
        (['map', 'nucid', ['set', 'nucname']], 
            ('map', ('int32', 'nucid'), ('set', ('str', 'nucname'), 0), 0)),
        (['intrange', 1, 2], ('int32', ('intrange', 
                                            ('low', 'int32', 1), 
                                            ('high', 'int32', 2)))), 
        (['nucrange', 92000, 93000], (('int32', 'nucid'), 
                                        ('nucrange', 
                                            ('low', ('int32', 'nucid'), 92000), 
                                            ('high', ('int32', 'nucid'), 93000)))), 
        (['range', 'int32', 1, 2], ('int32', 
                                        ('range', 'int32',
                                            ('low', 'int32', 1), 
                                            ('high', 'int32', 2)))), 
        (['range', 'nucid', 92000, 93000], (('int32', 'nucid'), 
                                        ('range', ('int32', 'nucid'),
                                            ('low', ('int32', 'nucid'), 92000), 
                                            ('high', ('int32', 'nucid'), 93000)))), 
    )
    for t, exp in cases:
        yield check_canon, t, exp            # Check that the case works,
        yield check_canon, ts.canon(t), exp  # And that it is actually canonical.


def check_cython_ctype(t, exp):
    obs = ts.cython_ctype(t)
    assert_equal(obs, exp)

@with_setup(add_new_refined, del_new_refined)
def test_cython_ctype():
    cases = (
        ('str', 'std_string'),
        (['str'], 'std_string'),
        ('f4', 'float'),
        ('nucid', 'int'),
        (['nucid'], 'int'), 
        (['set', 'complex'], 'cpp_set[extra_types.complex_t]'),
        (['map', 'nucid', 'float'], 'cpp_map[int, double]'),
        ('comp_map', 'cpp_map[int, double]'),
        (['char', '*'], 'char *'),
        (['char', 42], 'char [42]'),
        (['map', 'nucid', ['set', 'nucname']], 'cpp_map[int, cpp_set[std_string]]'),
        (['intrange', 1, 2], 'int'), 
        (['nucrange', 92000, 93000], 'int'),
        (['range', 'int32', 1, 2], 'int'), 
        (['range', 'nucid', 92000, 93000], 'int'), 
    )
    for t, exp in cases:
        yield check_cython_ctype, t, exp  # Check that the case works,


def check_cython_cimport_tuples_no_cy(t, exp):
    obs = ts.cython_cimport_tuples(t, inc=set(['c']))
    assert_equal(obs, exp)

@with_setup(add_new_refined, del_new_refined)
def test_cython_cimport_tuples_no_cy():
    cases = (
        ('str', set([('libcpp.string', 'string', 'std_string')])),
        (['str'], set([('libcpp.string', 'string', 'std_string')])),
        ('f4', set()),
        ('nucid', set()),
        (['nucid'], set()), 
        (['set', 'complex'], set([('libcpp.set', 'set', 'cpp_set'), 
                                  ('pyne', 'extra_types')])),
        (['map', 'nucid', 'float'], set([('libcpp.map', 'map', 'cpp_map')])),
        ('comp_map', set([('libcpp.map', 'map', 'cpp_map')])),
        (['char', '*'], set()),
        (['char', 42], set()),
        (['map', 'nucid', ['set', 'nucname']], set([('libcpp.set', 'set', 'cpp_set'),
                                                    ('libcpp.map', 'map', 'cpp_map'),
                                        ('libcpp.string', 'string', 'std_string')])),
        (['intrange', 1, 2], set()), 
        (['nucrange', 92000, 93000], set()),
        (['range', 'int32', 1, 2], set()), 
        (['range', 'nucid', 92000, 93000], set()), 
    )
    for t, exp in cases:
        yield check_cython_cimport_tuples_no_cy, t, exp  # Check that the case works


def check_cython_cimport_tuples_with_cy(t, exp):
    obs = ts.cython_cimport_tuples(t)
    assert_equal(obs, exp)

@with_setup(add_new_refined, del_new_refined)
def test_cython_cimport_tuples_with_cy():
    cases = (
        ('str', set([('libcpp.string', 'string', 'std_string')])),
        (['str'], set([('libcpp.string', 'string', 'std_string')])),
        ('f4', set()),
        ('nucid', set()),
        (['nucid'], set()), 
        (['set', 'complex'], set([('libcpp.set', 'set', 'cpp_set'), 
                                  ('pyne', 'stlconverters', 'conv'),
                                  ('pyne', 'extra_types')])),
        (['map', 'nucid', 'float'], set([('libcpp.map', 'map', 'cpp_map'), 
                                         ('pyne', 'stlconverters', 'conv')])),
        ('comp_map', set([('libcpp.map', 'map', 'cpp_map'), 
                          ('pyne', 'stlconverters', 'conv')])),
        (['char', '*'], set()),
        (['char', 42], set()),
        (['map', 'nucid', ['set', 'nucname']], set([('libcpp.set', 'set', 'cpp_set'),
                                                    ('libcpp.map', 'map', 'cpp_map'),
                                                    ('pyne', 'stlconverters', 'conv'),
                                        ('libcpp.string', 'string', 'std_string')])),
        (['intrange', 1, 2], set()), 
        (['nucrange', 92000, 93000], set()),
        (['range', 'int32', 1, 2], set()), 
        (['range', 'nucid', 92000, 93000], set()), 
    )
    for t, exp in cases:
        yield check_cython_cimport_tuples_with_cy, t, exp  # Check that the case works


def check_cython_cimports(t, exp):
    obs = ts.cython_cimports(t)
    assert_equal(obs, exp)

def test_cython_cimports():
    cases = (
        # type checks
        ('str', set(['from libcpp.string cimport string as std_string',])),
        ('f4', set()),
        # seen checks
        (set([('orly',)]), set(['cimport orly',])),
        (set([('orly','yarly')]), set(['from orly cimport yarly',])),
        (set([('orly','yarly','nowai')]), set(['from orly cimport yarly as nowai',])),
        (set([('orly',), ('orly','yarly')]), 
            set(['cimport orly', 'from orly cimport yarly'])),
    )
    for t, exp in cases:
        yield check_cython_cimports, t, exp  # Check that the case works,


def check_cython_cytype(t, exp):
    obs = ts.cython_cytype(t)
    assert_equal(obs, exp)

@with_setup(add_new_refined, del_new_refined)
def test_cython_cytype():
    cases = (
        ('str', 'char *'),
        (['str'], 'char *'),
        ('f4', 'float'),
        ('nucid', 'int'),
        (['nucid'], 'int'), 
        (['set', 'complex'], 'conv._SetComplex'),
        (['map', 'nucid', 'float'], 'conv._MapIntDouble'),
        ('comp_map', 'conv._MapIntDouble'),
        (['char', '*'], 'char *'),
        (['char', 42], 'char [42]'),
        (['map', 'nucid', ['set', 'nucname']], 'conv._MapIntSetStr'),
        (['intrange', 1, 2], 'int'), 
        (['nucrange', 92000, 93000], 'int'),
        (['range', 'int32', 1, 2], 'int'), 
        (['range', 'nucid', 92000, 93000], 'int'), 
    )
    for t, exp in cases:
        yield check_cython_cytype, t, exp  # Check that the case works,


def check_cython_pytype(t, exp):
    obs = ts.cython_pytype(t)
    assert_equal(obs, exp)

@with_setup(add_new_refined, del_new_refined)
def test_cython_pytype():
    cases = (
        ('str', 'str'),
        (['str'], 'str'),
        ('f4', 'float'),
        ('nucid', 'int'),
        (['nucid'], 'int'), 
        (['set', 'complex'], 'conv.SetComplex'),
        (['map', 'nucid', 'float'], 'conv.MapIntDouble'),
        ('comp_map', 'conv.MapIntDouble'),
        (['char', '*'], 'str'),
        (['char', 42], 'str'),
        (['map', 'nucid', ['set', 'nucname']], 'conv.MapIntSetStr'),
        (['intrange', 1, 2], 'int'), 
        (['nucrange', 92000, 93000], 'int'),
        (['range', 'int32', 1, 2], 'int'), 
        (['range', 'nucid', 92000, 93000], 'int'), 
    )
    for t, exp in cases:
        yield check_cython_pytype, t, exp  # Check that the case works,


def check_cython_c2py(name, t, exp):
    obs = ts.cython_c2py(name, t)
    assert_equal(obs, exp)

@with_setup(add_new_refined, del_new_refined)
def test_cython_c2py():
    cases = (
        (('llama', 'str'), (None, None, 'str(<char *> llama.c_str())')),
        (('llama', ['str']), (None, None, 'str(<char *> llama.c_str())')),
        (('llama', 'f4'), (None, None, 'float(llama)')),
        (('llama', 'nucid'), (None, None, 'int(llama)')),
        (('llama', ['nucid']), (None, None, 'int(llama)')), 
        (('llama', ['set', 'complex']), ('cdef conv._SetComplex llama_proxy', 
            ('if self._llama is None:\n'
             '    llama_proxy = conv.SetComplex(False, False)\n'
             '    llama_proxy.set_ptr = &self._inst.llama\n'
             '    self._llama = llama_proxy\n'), 'self._llama')),
        (('llama', ['map', 'nucid', 'float']), ('cdef conv._MapIntDouble llama_proxy', 
            ('if self._llama is None:\n'
             '    llama_proxy = conv.MapIntDouble(False, False)\n'
             '    llama_proxy.map_ptr = &self._inst.llama\n'
             '    self._llama = llama_proxy\n'), 'self._llama')),
        (('llama', 'comp_map'), ('cdef conv._MapIntDouble llama_proxy', 
            ('if self._llama is None:\n'
             '    llama_proxy = conv.MapIntDouble(False, False)\n'
             '    llama_proxy.map_ptr = &self._inst.llama\n'
             '    self._llama = llama_proxy\n'), 'self._llama')),
        (('llama', ['char', '*']), (None, None, 'str(<char *> llama)')),
        (('llama', ['char', 42]), (None, None, 'str(<char *> llama)')),
        (('llama', ['map', 'nucid', ['set', 'nucname']]), 
            ('cdef conv._MapIntSetStr llama_proxy', 
            ('if self._llama is None:\n'
             '    llama_proxy = conv.MapIntSetStr(False, False)\n'
             '    llama_proxy.map_ptr = &self._inst.llama\n'
             '    self._llama = llama_proxy\n'), 'self._llama')),
        (('llama', ['intrange', 1, 2]), (None, None, 'int(llama)')), 
        (('llama', ['nucrange', 92000, 93000]), (None, None, 'int(llama)')),
        (('llama', ['range', 'int32', 1, 2]), (None, None, 'int(llama)')), 
        (('llama', ['range', 'nucid', 92000, 93000]), (None, None, 'int(llama)')), 
    )
    for (name, t), exp in cases:
        yield check_cython_c2py, name, t, exp  # Check that the case works,

