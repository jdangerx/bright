from traits.api import HasTraits, Str, Dict


class ClassModel(HasTraits):

    name = Str
    var = Str
    class_name = Str
    extra_data_parameter = Dict
    
    
    import_line = Str("from bright import {name}")
    instance_line = Str("{var} = {classname}()")
    calc_line = Str("{var}.calc({ms})")
    
    def add_import(self):
        return self.import_line.format(name=self.name)

    def add_instance(self):
        return self.instance_line.format(var = self.var, classname = self.class_name, mass_dict = self.extra_data_parameter)

    def add_calc(self, msname):
        return self.calc_line.format(var = self.var, ms = msname)
