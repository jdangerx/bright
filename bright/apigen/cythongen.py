"""Generates a Cython wrapper for various Bright classes from either
description dictionaries or from header files.
"""
import math
from copy import deepcopy

from bright.apigen.utils import indent, expand_default_args
from bright.apigen import typesystem as ts
from bright.apigen.typesystem import cython_ctype, cython_cimport_tuples, \
    cython_cimports, register_class, cython_cytype, cython_pytype, cython_c2py, \
    cython_py2c, cython_import_tuples, cython_imports, isrefinement

AUTOGEN_WARNING = \
"""################################################
#                 WARNING!                     #
# This file has been auto-generated by Bright. #
# Do not modify!!!                             #
#                                              #
#                                              #
#                    Come on, guys. I mean it! #
################################################
"""

_cpppxd_template = AUTOGEN_WARNING + \
"""{cimports}

cdef extern from "{header_filename}" namespace "{namespace}":

    cdef cppclass {name}{parents}:
        # constructors
{constructors_block}

        # attributes
{attrs_block}

        # methods
{methods_block}

{extra}
"""


def gencpppxd(desc, exception_type='+'):
    """Generates a cpp_*.pxd Cython header file for exposing C/C++ data from to 
    other Cython wrappers based off of a dictionary (desc)ription.
    """
    pars = ', '.join([cython_ctype(p) for p in desc['parents'] or ()])
    d = {'parents': pars if 0 == len(pars) else '('+pars+')'}
    copy_from_desc = ['name', 'namespace', 'header_filename']
    for key in copy_from_desc:
        d[key] = desc[key]
    inc = set(['c'])

    cimport_tups = set()
    for parent in desc['parents'] or ():
        cython_cimport_tuples(parent, cimport_tups, inc)

    alines = []
    attritems = sorted(desc['attrs'].items())
    for aname, atype in attritems:
        if aname.startswith('_'):
            continue
        alines.append("{0} {1}".format(cython_ctype(atype), aname))
        cython_cimport_tuples(atype, cimport_tups, inc)
    d['attrs_block'] = indent(alines, 8)

    mlines = []
    clines = []
    estr = str() if exception_type is None else  ' except {0}'.format(exception_type)
    methitems = sorted(expand_default_args(desc['methods'].items()))
    for mkey, mrtn in methitems:
        mname, margs = mkey[0], mkey[1:]
        if mname.startswith('_') or mname.startswith('~'):
            continue  # private or destructor
        argfill = ", ".join([cython_ctype(a[1]) for a in margs])
        for a in margs:
            cython_cimport_tuples(a[1], cimport_tups, inc)
        line = "{0}({1}){2}".format(mname, argfill, estr)
        if mrtn is None:
            # this must be a constructor
            if line not in clines:
                clines.append(line)
        else:
            # this is a normal method
            rtype = cython_ctype(mrtn)
            cython_cimport_tuples(mrtn, cimport_tups, inc)
            line = rtype + " " + line
            if line not in mlines:
                mlines.append(line)
    d['methods_block'] = indent(mlines, 8)
    d['constructors_block'] = indent(clines, 8)

    d['cimports'] = "\n".join(sorted(cython_cimports(cimport_tups)))
    d['extra'] = desc.get('extra', {}).get('cpppxd', '')
    cpppxd = _cpppxd_template.format(**d)
    if 'cpppxd_filename' not in desc:
        desc['cpppxd_filename'] = 'cpp_{0}.pxd'.format(d['name'].lower())
    return cpppxd
    


_pxd_template = AUTOGEN_WARNING + \
"""{cimports}

cdef class {name}{parents}:
{body}

{extra}
"""


def genpxd(desc):
    """Generates a *.pxd Cython header file for exposing C/C++ data from to 
    other Cython wrappers based off of a dictionary (desc)ription.
    """
    if 'pxd_filename' not in desc:
        desc['pxd_filename'] = '{0}.pxd'.format(desc['name'].lower())
    pars = ', '.join([cython_cytype(p) for p in desc['parents'] or ()])
    d = {'parents': pars if 0 == len(pars) else '('+pars+')'}
    copy_from_desc = ['name',]
    for key in copy_from_desc:
        d[key] = desc[key]

    cimport_tups = set()
    for parent in desc['parents'] or ():
        cython_cimport_tuples(parent, cimport_tups, set(['cy']))

    from_cpppxd = desc['cpppxd_filename'].rsplit('.', 1)[0]
    # This is taken care of in main!
    #register_class(desc['name'], cython_cimport=from_cpppxd,
    #               cython_c_type="{0}.{1}".format(from_cpppxd, desc['name']),)
    d['name_type'] = cython_ctype(desc['name'])
    cython_cimport_tuples(desc['name'], cimport_tups, set(['c']))

    parentless_body = ['cdef void * _inst', 'cdef public bint _free_inst'] 
    body = parentless_body if desc['parents'] is None else []
    attritems = sorted(desc['attrs'].items())
    for aname, atype in attritems:
        if aname.startswith('_'):
            continue  # skip private
        _, _, cachename, iscached = cython_c2py(aname, atype, cache_prefix=None)
        if iscached:
            cython_cimport_tuples(atype, cimport_tups)
            cyt = cython_cytype(atype)
            decl = "cdef public {0} {1}".format(cyt, cachename)
            body.append(decl)

    d['cimports'] = "\n".join(sorted(cython_cimports(cimport_tups)))
    d['body'] = indent(body or ['pass'])
    d['extra'] = desc.get('extra', {}).get('pxd', '')
    pxd = _pxd_template.format(**d)
    return pxd
    

_pyx_template = AUTOGEN_WARNING + \
'''"""{module_docstring}
"""
{cimports}

{imports}

cdef class {name}{parents}:
{class_docstring}

    # constuctors
    def __cinit__(self, *args, **kwargs):
        self._inst = NULL
        self._free_inst = True

        # cached property defaults
{property_defaults}

{constructor_block}

    # attributes
{attrs_block}
    # methods
{methods_block}

{extra}
'''

def _gen_property_get(name, t, cached_names=None, inst_name="self._inst"):
    """This generates a Cython property getter for a variable of a given 
    name and type."""
    lines = ['def __get__(self):']
    decl, body, rtn, iscached = cython_c2py(name, t, inst_name=inst_name)
    if decl is not None: 
        lines += indent(decl, join=False)
    if body is not None:
        lines += indent(body, join=False)
    if iscached and cached_names is not None:
        cached_names.append(rtn)
    lines += indent("return {0}".format(rtn), join=False)
    return lines


def _gen_property_set(name, t, inst_name="self._inst", cached_name=None):
    """This generates a Cython property setter for a variable of a given 
    name and type."""
    lines = ['def __set__(self, value):']
    decl, body, rtn = cython_py2c('value', t)
    if decl is not None: 
        lines += indent(decl, join=False)
    if body is not None:
        lines += indent(body, join=False)
    lines += indent("{0}.{1} = {2}".format(inst_name, name, rtn), join=False)
    if cached_name is not None:
        lines += indent("{0} = None".format(cached_name), join=False)
    return lines


def _gen_property(name, t, doc=None, cached_names=None, inst_name="self._inst"):
    """This generates a Cython property for a variable of a given name and type."""
    lines  = ['property {0}:'.format(name)] 
    lines += [] if doc is None else indent('\"\"\"{0}\"\"\"'.format(doc), join=False)
    oldcnlen = 0 if cached_names is None else len(cached_names)
    lines += indent(_gen_property_get(name, t, cached_names=cached_names, 
                                      inst_name=inst_name), join=False)
    lines += ['']
    newcnlen = 0 if cached_names is None else len(cached_names)
    cached_name = cached_names[-1] if newcnlen == 1 + oldcnlen else None
    lines += indent(_gen_property_set(name, t, inst_name=inst_name, 
                                      cached_name=cached_name), join=False)
    lines += ['', ""]
    return lines


def _gen_method(name, name_mangled, args, rtn, doc=None, inst_name="self._inst"):
    argfill = ", ".join(['self'] + [a[0] for a in args if 2 == len(a)] + \
                        ["{0}={1}".format(a[0], a[2]) for a in args if 3 == len(a)])
    lines  = ['def {0}({1}):'.format(name_mangled, argfill)]
    lines += [] if doc is None else indent('\"\"\"{0}\"\"\"'.format(doc), join=False)
    decls = []
    argbodies = []
    argrtns = {}
    for a in args:
        adecl, abody, artn = cython_py2c(a[0], a[1])
        if adecl is not None: 
            decls += indent(adecl, join=False)
        if abody is not None:
            argbodies += indent(abody, join=False)
        argrtns[a[0]] = artn
    rtype = cython_ctype(rtn)
    hasrtn = rtype not in set(['None', None, 'NULL', 'void'])
    argvals = ', '.join([argrtns[a[0]] for a in args])
    fcall = '{0}.{1}({2})'.format(inst_name, name, argvals)
    if hasrtn:
        fcdecl, fcbody, fcrtn, fccached = cython_c2py('rtnval', rtn, cached=False)
        decls += indent("cdef {0} {1}".format(rtype, 'rtnval'), join=False)
        func_call = indent('rtnval = {0}'.format(fcall), join=False)
        if fcdecl is not None: 
            decls += indent(fcdecl, join=False)
        if fcbody is not None:
            func_call += indent(fcbody, join=False)
        func_rtn = indent("return {0}".format(fcrtn), join=False)
    else:
        func_call = indent(fcall, join=False)
        func_rtn = []
    lines += decls
    lines += argbodies
    lines += func_call
    lines += func_rtn
    lines += ['', ""]
    return lines


def _gen_constructor(name, name_mangled, classname, args, doc=None, 
                     cpppxd_filename=None, inst_name="self._inst"):
    argfill = ", ".join(['self'] + [a[0] for a in args if 2 == len(a)] + \
                        ["{0}={1}".format(a[0], a[2]) for a in args if 3 == len(a)])
    lines  = ['def {0}({1}):'.format(name_mangled, argfill)]
    lines += [] if doc is None else indent('\"\"\"{0}\"\"\"'.format(doc), join=False)
    decls = []
    argbodies = []
    argrtns = {}
    for a in args:
        adecl, abody, artn = cython_py2c(a[0], a[1])
        if adecl is not None: 
            decls += indent(adecl, join=False)
        if abody is not None:
            argbodies += indent(abody, join=False)
        argrtns[a[0]] = artn
    argvals = ', '.join([argrtns[a[0]] for a in args])
    classname = classname if cpppxd_filename is None else \
                    "{0}.{1}".format(cpppxd_filename.rsplit('.', 1)[0], classname)
    fcall = 'self._inst = new {0}({1})'.format(classname, argvals)
    func_call = indent(fcall, join=False)
    lines += decls
    lines += argbodies
    lines += func_call
    lines += ['', ""]
    return lines

def _gen_dispatcher(name, name_mangled, doc=None, hasrtn=True):
    argfill = ", ".join(['self', '*args', '**kwargs'])
    lines  = ['def {0}({1}):'.format(name, argfill)]
    lines += [] if doc is None else indent('\"\"\"{0}\"\"\"'.format(doc), join=False)
    types = ["types = set([(i, type(a)) for i, a in enumerate(args)])",
             "types.update([(k, type(v)) for k, v in kwargs.iteritems()])",]
    lines += indent(types, join=False)
    refinenum = lambda x: (sum([int(isrefinement(a[1])) for a in x[0][1:]]), len(x[0]), x[1])
    mangitems = sorted(name_mangled.items(), key=refinenum)
    mtypeslines = []
    lines += indent("# vtable-like dispatch for exactly matching types", join=False)
    for key, mangled_name in mangitems:
        cargs = key[1:]
        arang = range(len(cargs))
        anames = [ca[0] for ca in cargs]
        pytypes = [cython_pytype(ca[1]) for ca in cargs]
        mtypes = ", ".join(
            ["({0}, {1})".format(i, pyt) for i, pyt in zip(arang, pytypes)] + \
            ['("{0}", {1})'.format(n, pyt) for n, pyt in zip(anames, pytypes)])
        mtups = '(' + mtypes + ')' if 0 < len(mtypes) else mtypes
        mtypeslines.append(mangled_name + "_argtypes = frozenset(" + mtups + ")")
        cond = ["if types <= self.{0}_argtypes:".format(mangled_name),]
        if hasrtn:
            rline = "return self.{0}(*args, **kwargs)".format(mangled_name)
        else:
            rline = ["self.{0}(*args, **kwargs)".format(mangled_name), "return"]
        cond += indent(rline, join=False)
        lines += indent(cond, join=False)
    lines = sorted(mtypeslines) + [''] +  lines
    lines += indent("# duck-typed dispatch based on whatever works!", join=False)
    refineopp = lambda x: (-1*sum([int(isrefinement(a[1])) for a in x[0][1:]]), len(x[0]), x[1])
    mangitems = sorted(name_mangled.items(), key=refineopp)
    for key, mangled_name in mangitems:
        lines += indent('try:', join=False)
        if hasrtn:
            rline = "return self.{0}(*args, **kwargs)".format(mangled_name)
        else:
            rline = ["self.{0}(*args, **kwargs)".format(mangled_name), "return"]
        lines += indent(indent(rline, join=False), join=False)
        lines += indent(["except (RuntimeError, TypeError, NameError):",
                         indent("pass", join=False)[0],], join=False)
    errmsg = "raise RuntimeError('method {0}() could not be dispatched')".format(name)
    lines += indent(errmsg, join=False)
    lines += ['']
    return lines


def _method_instance_names(desc, env, key, rtn):
    classnames = (desc['parents'] or []) + [desc['name']]
    for classname in classnames:
        classrtn = env.get(classname, {}).get('methods', {}).get(key, NotImplemented)
        if rtn != classrtn:
            continue
        #class_ctype = cython_ctype(desc['name'])
        class_ctype = cython_ctype(classname)
        inst_name = "(<{0} *> self._inst)".format(class_ctype)
        return inst_name, classname
    return "(<{0} *> self._inst)".format(desc['name']), desc['name']


def _count0(x):
    c = {}
    for v in x:
        v0 = v[0]
        c[v0] = c.get(v0, 0) + 1
    return c

def genpyx(desc, env=None):
    """Generates a *.pyx Cython wrapper implementation for exposing a C/C++ 
    class based off of a dictionary (desc)ription.  The (env)ironment is a 
    dictionary of all class names known to their descriptions.
    """
    if env is None:
        env = {desc['name']: desc}
    nodocmsg = "no docstring for {0}, please file a bug report!"
    pars = ', '.join([cython_cytype(p) for p in desc['parents'] or ()])
    d = {'parents': pars if 0 == len(pars) else '('+pars+')'}
    copy_from_desc = ['name', 'namespace', 'header_filename']
    for key in copy_from_desc:
        d[key] = desc[key]
    d['module_docstring'] = desc.get('docstrings', {})\
                                .get('module', nodocmsg.format(desc['name'].lower()))
    class_doc = desc.get('docstrings', {}).get('class', nodocmsg.format(desc['name']))
    d['class_docstring'] = indent('\"\"\"{0}\"\"\"'.format(class_doc))

    class_ctype = cython_ctype(desc['name'])
    inst_name = "(<{0} *> self._inst)".format(class_ctype)

    import_tups = set()
    cimport_tups = set()
    for parent in desc['parents'] or ():
        cython_import_tuples(parent, import_tups)
        cython_cimport_tuples(parent, cimport_tups)

    alines = []
    cached_names = []
    attritems = sorted(desc['attrs'].items())
    for aname, atype in attritems:
        if aname.startswith('_'):
            continue  # skip private
        adoc = desc.get('docstrings', {}).get('attrs', {})\
                                         .get(aname, nodocmsg.format(aname))
        alines += _gen_property(aname, atype, adoc, cached_names=cached_names, 
                                inst_name=inst_name)
        cython_import_tuples(atype, import_tups)
        cython_cimport_tuples(atype, cimport_tups)
    d['attrs_block'] = indent(alines)
    pd = ["{0} = None".format(n) for n in cached_names]
    d['property_defaults'] = indent(indent(pd, join=False))

    mlines = []
    clines = []
    methcounts = _count0(desc['methods'])
    currcounts = {k: 0 for k in methcounts}
    mangled_mnames = {}
    methitems = sorted(desc['methods'].items())
    for mkey, mrtn in methitems:
        mname, margs = mkey[0], mkey[1:]
        if mname.startswith('_'):
            continue  # skip private
        if 1 < methcounts[mname]:
            mname_mangled = "_{0}_{1}_{2:0{3}}".format(desc['name'], mname, 
                    currcounts[mname], int(math.log(methcounts[mname], 10)+1)).lower()
        else:
            mname_mangled = mname
        currcounts[mname] += 1
        mangled_mnames[mkey] = mname_mangled
        for a in margs:
            cython_import_tuples(a[1], import_tups)
            cython_cimport_tuples(a[1], cimport_tups)
        minst_name, mcname = _method_instance_names(desc, env, mkey, mrtn)
        if mcname != desc['name']:
            cython_import_tuples(mcname, import_tups)
            cython_cimport_tuples(mcname, cimport_tups)
        if mrtn is None:
            # this must be a constructor
            if mname not in (desc['name'], '__init__'):
                continue  # skip destuctors
            if 1 == methcounts[mname]:
                mname_mangled = '__init__'
                mangled_mnames[mkey] = mname_mangled
            mdoc = desc.get('docstrings', {}).get('methods', {}).get(mname, '')
            clines += _gen_constructor(mname, mname_mangled, 
                                       desc['name'], margs, doc=mdoc, 
                                       cpppxd_filename=desc['cpppxd_filename'],
                                       inst_name=minst_name)
            if 1 < methcounts[mname] and currcounts[mname] == methcounts[mname]:
                # write dispatcher
                nm = {k: v for k, v in mangled_mnames.iteritems() if k[0] == mname}
                clines += _gen_dispatcher('__init__', nm, doc=mdoc, hasrtn=False)
        else:
            # this is a normal method
            cython_import_tuples(mrtn, import_tups)
            cython_cimport_tuples(mrtn, cimport_tups)
            mdoc = desc.get('docstrings', {}).get('methods', {})\
                                             .get(mname, nodocmsg.format(mname))
            mlines += _gen_method(mname, mname_mangled, margs, mrtn, mdoc, 
                                  inst_name=minst_name)
            if 1 < methcounts[mname] and currcounts[mname] == methcounts[mname]:
                # write dispatcher
                nm = {k: v for k, v in mangled_mnames.iteritems() if k[0] == mname}
                mlines += _gen_dispatcher(mname, nm, doc=mdoc)
    if desc['parents'] is None:
        clines += ["def __dealloc__(self):"]
        clines += indent("if self._free_inst:", join=False)
        clines += indent(indent("free(self._inst)", join=False), join=False)
        cimport_tups.add(('libc.stdlib', 'free'))

    d['methods_block'] = indent(mlines)
    d['constructor_block'] = indent(clines)

    d['imports'] = "\n".join(sorted(cython_imports(import_tups)))
    d['cimports'] = "\n".join(sorted(cython_cimports(cimport_tups)))
    if 'numpy' in d['cimports']:
        d['imports'] += "\n\nnp.import_array()"
    d['extra'] = desc.get('extra', {}).get('pyx', '')
    pyx = _pyx_template.format(**d)
    if 'pyx_filename' not in desc:
        desc['pyx_filename'] = '{0}.pyx'.format(d['name'].lower())
    return pyx
