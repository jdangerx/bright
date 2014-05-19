from __future__ import print_function

import os
import re
from lxml import etree
import pyne.nucname as nucname

class NCodeOpenMC(object):
    def __init__(self, rc):
        self.rc = rc
        self.write_functions = {"lwr": self.write_lwr}
        self.get_possible_xs()

    def write(self):
        self.write_functions[self.rc.reactor]()

    def dict_to_xml(self, d, filter_ls = None, elt = None):
        xml = ""
        for key in d:
            if filter_ls and not (key in filter_ls):
                continue
            if elt == None:
                key_elt = etree.Element(key)
            else:
                key_elt = etree.SubElement(elt, key)
            if isinstance(d[key], dict):
                self.dict_to_xml(d[key], elt=key_elt)
            else:
                key_elt.text=str(d[key])
            key_elt_str = etree.tostring(key_elt, pretty_print=True)
            xml = "\n".join((xml, key_elt_str))
        return xml.strip()

    def write_lwr(self):
        self.write_geometry_lwr()

    def write_geometry_lwr(self):
        geom = etree.Element("geometry")
        cell_1 = etree.SubElement(geom, "cell", id="1", fill="5", surfaces="1 -2 3 -4")
        cell_101 = etree.SubElement(geom, "cell", id="101", universe="1", material="1", surfaces="-5")
        cell_102 = etree.SubElement(geom, "cell", id="102", universe="1", material="2", surfaces="5")
        cell_201 = etree.SubElement(geom, "cell", id="201", universe="2", material="1", surfaces="-6")
        cell_202 = etree.SubElement(geom, "cell", id="202", universe="2", material="2", surfaces="6")
        cell_301 = etree.SubElement(geom, "cell", id="301", universe="3", material="1", surfaces="-7")
        cell_302 = etree.SubElement(geom, "cell", id="302", universe="3", material="2", surfaces="7")
        lattice = etree.SubElement(geom, "lattice", id="5")
        lattice_type = etree.SubElement(lattice, "type")
        lattice_type.text = "rectangular"
        dimension = etree.SubElement(lattice, "dimension")
        dimension.text = "4 4"
        lower_left = etree.SubElement(lattice, "lower_left")
        lower_left.text = "-2.0 -2.0"
        width = etree.SubElement(lattice, "width")
        if self.rc.get("unit_cell_pitch") != None:
            width.text = "{} {}".format(self.rc.get("unit_cell_pitch"))
        else:
            width.text = "1.0 1.0"
        universes = etree.SubElement(lattice, "universes")
        universes.text="""
        1 2 1 2
        2 3 2 3
        1 2 1 2
        2 3 2 3
        """
        surface_1 = etree.SubElement(geom, "surface", id="1", type="x-plane", coeffs="-2.0", boundary="vacuum")
        surface_2 = etree.SubElement(geom, "surface", id="2", type="x-plane", coeffs="2.0", boundary="vacuum")
        surface_3 = etree.SubElement(geom, "surface", id="3", type="y-plane", coeffs="-2.0", boundary="vacuum")
        surface_4 = etree.SubElement(geom, "surface", id="4", type="y-plane", coeffs="2.0", boundary="vacuum")
        surface_5 = etree.SubElement(geom, "surface", id="5", type="z-cylinder", coeffs="0.0 0.0 0.4")
        surface_6 = etree.SubElement(geom, "surface", id="6", type="z-cylinder", coeffs="0.0 0.0 0.3")
        surface_7 = etree.SubElement(geom, "surface", id="7", type="z-cylinder", coeffs="0.0 0.0 0.2")
        # are these the right places to be putting the cell radii?
        if self.rc.get("fuel_cell_radius"):
            surface_7.set("coeffs", "0.0 0.0 0.0 {}".format(self.rc.get("fuel_cell_radius")))
        if self.rc.get("void_cell_radius"):
            surface_6.set("coeffs", "0.0 0.0 0.0 {}".format(self.rc.get("void_cell_radius")))
        if self.rc.get("clad_cell_radius"):
            surface_5.set("coeffs", "0.0 0.0 0.0 {}".format(self.rc.get("clad_cell_radius")))
        self.geometry_xml = etree.tostring(geom, pretty_print=True)
        print(self.geometry_xml)

    def get_possible_xs(self):
        if os.path.isfile(os.environ['CROSS_SECTIONS']):
            xs_etree = etree.parse(os.environ['CROSS_SECTIONS'])
            self.xs_data = xs_etree.findall(".//*[@temperature]")
        # TODO: alternate xs locations

    def diff_temps(self, xs, exp_temp):
        mev_to_kelvin = 293.6/2.53e-8
        return abs(float(xs.get("temperature"))*mev_to_kelvin - exp_temp)

    def get_closest_xs_to_temp(self, name, temp):
        only_correct_name = [xs for xs in self.xs_data if xs.get("name").split(".")[0] == name]
        print(name)
        tol = min([self.diff_temps(xs, temp) for xs in only_correct_name])
        return [xs.get("name") for xs in only_correct_name if self.diff_temps(xs, temp) <= tol][0]

    def write_materials_lwr(self):
        materials = etree.Element("materials")
        def_temp = self.rc.get("temperature")

        fuel = etree.SubElement(materials, "material", id="1")
        fuel_density = etree.SubElement(fuel, "density", value="4.5", units="g/cc")
        ihm = self.rc.get("initial_heavy_metal")
        for nuc in ihm:
            hyphenated_name = re.sub(r"(\d+)", r"-\1", nucname.name(nuc))
            nuc_zaid = str(nucname.zzaaam(nuc)/10)
            xs_id = self.get_closest_xs_to_temp(nuc_zaid, def_temp).split(".")[-1]
            nuc_tag = etree.SubElement(fuel, "nuclide", name=hyphenated_name, mo=str(ihm[nuc]), xs=xs_id)

        coolant = etree.SubElement(materials, "material", id="2")
        coolant_comp = self.rc.get("coolant")
        coolant_density = etree.SubElement(coolant, "density", value="1.0", units="g/cc")
        for nuc in coolant_comp:
            hyphenated_name = re.sub(r"(\d+)", r"-\1", nucname.name(nuc))
            nuc_zaid = str(nucname.zzaaam(nuc)/10)
            xs_id = self.get_closest_xs_to_temp(nuc_zaid, def_temp).split(".")[-1]
            nuc_tag = etree.SubElement(fuel, "nuclide", name=hyphenated_name, ao=str(coolant_comp[nuc]), xs=xs_id)

        coolant_sab_name = self.rc.get("coolant_sab")
        coolant_sab_xs_id = self.get_closest_xs_to_temp(nuc_zaid, def_temp).split(".")[-1]
        coolant_sab = etree.SubElement(coolant, "sab", name=coolant_sab_name, xs=coolant_sab_xs_id)
        self.materials_xml = etree.tostring(materials, pretty_print=True)
        print(self.materials_xml)

    def make_settings(self):
        opts = ["confidence_intervals", "cross_sections", "cutoff", "eigenvalue",
                "energy_grid", "entropy", "fixed_source", "no_reduce", "output",
                "output_path", "ptables", "run_cmfd", "seed", "source",
                "state_point", "source_point", "survival_biasing", "threads",
                "trace", "track", "uniform_fs", "verbosity"]
        self.settings_xml = self.dict_to_xml(self.rc, filter_ls = opts)

    def write_tallies(self):
        pass

    def write_plots(self):
        pass

    def write_cmfd(self):
        pass
