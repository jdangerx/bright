from __future__ import print_function

import os
import re
import subprocess as sp
from lxml import etree
import pyne.nucname as nucname

class NCodeOpenMC(object):
    def __init__(self, rc = {"reactor":"lwr"}):
        self.rc = rc
        default_rcs = {"lwr": "defxsgen_lwr.py"}
        reactor_rc = {}
        execfile(default_rcs[self.rc["reactor"]], reactor_rc)
        self.rc.update(reactor_rc)
        self.write_functions = {"lwr": self.write_lwr}
        self.get_possible_xs()

    def write(self):
        self.write_functions[self.rc["reactor"]]()
        self.write_tallies()
        self.write_settings()

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
        self.write_materials_lwr()

    def write_geometry_lwr(self):
        geom = etree.Element("geometry")
        cell_1 = etree.SubElement(geom, "cell", id="1", fill="5", surfaces="1 -2 3 -4")

        # universe 1 is a universe with a fuel pin
        fuel_pin = etree.SubElement(geom, "cell", id="101", universe="1", material="1", surfaces="-5")
        void = etree.SubElement(geom, "cell", id="201", universe="1", material="void", surfaces="5 -6")
        cladding = etree.SubElement(geom, "cell", id="301", universe="1", material="2", surfaces="6 -7")
        coolant = etree.SubElement(geom, "cell", id="401", universe="1", material="3", surfaces="7")

        #universe 2 is a universe with just coolant everywhere!
        coolant_cell = etree.SubElement(geom, "cell", id="501", universe="2", material="3", surfaces="1 -1")

        # lattice description
        lattice = etree.SubElement(geom, "lattice", id="5")
        lattice_type = etree.SubElement(lattice, "type")
        lattice_type.text = "rectangular"
        dimension = etree.SubElement(lattice, "dimension")
        dimension.text = self.rc["dimension"]
        lower_left = etree.SubElement(lattice, "lower_left")
        lower_left.text = self.rc["lower_left"]
        width = etree.SubElement(lattice, "width")
        width.text = "{0} {0}".format(self.rc["unit_cell_pitch"])
        universes = etree.SubElement(lattice, "universes")
        universes.text= self.rc["lattice_str"]


        xmin, ymin = self.rc["lower_left"].split()
        xwidth, ywidth = self.rc["dimension"].split()
        xmax = str(float(xmin)+float(xwidth))
        ymax = str(float(ymin)+float(ywidth))

        x_low_surface = etree.SubElement(geom, "surface", id="1", type="x-plane", coeffs=xmin, boundary="vacuum")
        x_high_surface = etree.SubElement(geom, "surface", id="2", type="x-plane", coeffs=xmax, boundary="vacuum")
        y_low_surface = etree.SubElement(geom, "surface", id="3", type="y-plane", coeffs=ymin, boundary="vacuum")
        y_high_surface = etree.SubElement(geom, "surface", id="4", type="y-plane", coeffs=ymax, boundary="vacuum")

        fuel_void = etree.SubElement(geom, "surface", id="5", type="z-cylinder")
        fuel_void.set("coeffs", "0.0 0.0 0.0 {}".format(self.rc["fuel_cell_radius"]))
        void_clad = etree.SubElement(geom, "surface", id="6", type="z-cylinder")
        void_clad.set("coeffs", "0.0 0.0 0.0 {}".format(self.rc["void_cell_radius"]))
        clad_cool = etree.SubElement(geom, "surface", id="7", type="z-cylinder")
        clad_cool.set("coeffs", "0.0 0.0 0.0 {}".format(self.rc["clad_cell_radius"]))
        self.geometry_xml = etree.tostring(geom, pretty_print=True)
        # with open("geometry.xml", "w") as f:
            # f.write(self.geometry_xml)
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
        tol = min([self.diff_temps(xs, temp) for xs in only_correct_name])
        return [xs.get("name") for xs in only_correct_name if self.diff_temps(xs, temp) <= tol][0]

    def write_materials_lwr(self):
        materials = etree.Element("materials")
        temp = self.rc["temperature"]

        fuel = etree.SubElement(materials, "material", id="1")
        fuel_density = etree.SubElement(fuel, "density", value="4.5", units="g/cm3")
        ihm = self.rc["initial_heavy_metal"]
        for nuc in ihm:
            hyphenated_name = re.sub(r"(\d+)", r"-\1", nucname.name(nuc))
            nuc_zaid = str(nucname.zzaaam(nuc)/10)
            xs_id = self.get_closest_xs_to_temp(nuc_zaid, temp).split(".")[-1]
            nuc_tag = etree.SubElement(fuel, "nuclide", name=hyphenated_name, mo=str(ihm[nuc]), xs=xs_id)

        cladding = etree.SubElement(materials, "material", id="2")
        cladding_comp = self.rc["cladding"]
        cladding_density = etree.SubElement(cladding, "density", value=str(self.rc["clad_density"]), units="g/cm3")
        for nuc in cladding_comp:
            hyphenated_name = re.sub(r"(\d+)", r"-\1", nucname.name(nuc))
            nuc_zaid = str(nucname.zzaaam(nuc)/10)
            xs_id = self.get_closest_xs_to_temp(nuc_zaid, temp).split(".")[-1]
            nuc_tag = etree.SubElement(cladding, "nuclide", name=hyphenated_name, ao=str(cladding_comp[nuc]), xs=xs_id)

        coolant = etree.SubElement(materials, "material", id="3")
        coolant_comp = self.rc["coolant"]
        coolant_density = etree.SubElement(coolant, "density", value=str(self.rc["cool_density"]), units="g/cm3")
        for nuc in coolant_comp:
            hyphenated_name = re.sub(r"(\d+)", r"-\1", nucname.name(nuc))
            nuc_zaid = str(nucname.zzaaam(nuc)/10)
            xs_id = self.get_closest_xs_to_temp(nuc_zaid, temp).split(".")[-1]
            nuc_tag = etree.SubElement(coolant, "nuclide", name=hyphenated_name, ao=str(coolant_comp[nuc]), xs=xs_id)

        coolant_sab_name = self.rc["coolant_sab"]
        coolant_sab_xs_id = self.get_closest_xs_to_temp(nuc_zaid, temp).split(".")[-1]
        coolant_sab = etree.SubElement(coolant, "sab", name=coolant_sab_name, xs=coolant_sab_xs_id)

        self.materials_xml = etree.tostring(materials, pretty_print=True)
        print(self.materials_xml)

    def write_settings(self):
        opts = ["confidence_intervals", "cross_sections", "cutoff", "eigenvalue",
                "energy_grid", "entropy", "fixed_source", "no_reduce", "output",
                "output_path", "ptables", "run_cmfd", "seed", "source",
                "state_point", "source_point", "survival_biasing", "threads",
                "trace", "track", "uniform_fs", "verbosity"]
        self.settings_xml = self.dict_to_xml(self.rc, filter_ls = opts)
        print(self.settings_xml)

    def write_tallies(self):
        tally = etree.Element("tally")
        nuclides = etree.SubElement(tally, "nuclides")
        nuclides.text = " ".join(self.rc["core_load_isos"])
        self.tallies_xml = etree.tostring(tally, pretty_print=True)
        print(self.tallies_xml)

    def write_plots(self):
        pass

    def write_cmfd(self):
        pass
