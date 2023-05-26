import argparse
import os
import sys
import time
from time import sleep
from typing import Union

import numpy as np
from lxml.etree import tostring
from matplotlib import pyplot as plt
from pyPanair.preprocess import wgs_creator
from scipy.optimize import root_scalar

from utils.ISACalculator import ISACalculator
from utils.XMLutils import DataHandler
from utils.read_output import read_output
from utils.transform_airfoils import transform_airfoils


class WingGeometryConditions:
    """
    Class to generate and store wing geometry parameters and flight conditions
    """
    def __init__(self, var_Sref, var_AR, var_lambda, var_Lambda, var_tc, var_h, var_M, var_WT, var_Theta, var_chord_kink=0, var_sweep_outboard=0, var_b_kink=0):
        feettom = 0.3048

        #Set initial parameters, metric units
        self.Sref = var_Sref*(feettom**2)
        self.AR = var_AR
        self.taper = var_lambda
        self.sweep = var_Lambda
        self.tc_ratio = var_tc
        self.twist = var_Theta

        self.h = var_h*feettom
        self.M = var_M
        self.WT = var_WT

        #Determine wing planform
        self.Swing = self.Sref/2
        self.b = np.sqrt(self.AR * self.Sref)
        self.half_span = self.b/2
        self.root_chord = 2 * self.Sref / (self.b * (1 + self.taper))
        self.tip_chord = self.root_chord * self.taper
        self.MAC = 2 / 3 * self.root_chord * ((1 + self.taper + self.taper ** 2) / (1 + self.taper))

        #Determine atmospheric conditions and cruise lift coefficient
        T, p, rho, a = ISACalculator(self.h)
        self.T = T
        self.p = p
        self.rho = rho
        self.a = a

        self.V = self.a * self.M
        self.cl = self.WT/(0.5*self.rho*self.V**2*self.Sref)

        #kink parameters
        if not (var_chord_kink == 0 and var_b_kink == 0 and var_sweep_outboard == 0):
            self.b_kink = var_b_kink*self.half_span
            self.sweep_outboard = var_sweep_outboard
            self.chord_kink = var_chord_kink*self.root_chord
            self.HasKink = True
        else:
            self.HasKink = False

    def prepare_input(self, airfoils, path_title_dict):
        self.airfoils = airfoils

        for k, v in path_title_dict.items():
            setattr(self, k, v)

        ### Paths and names ###
        airfoil_folder_path = os.path.join(self.panair_folder_path, 'Airfoils')

        ### Create Unix paths for panair ###
        path_elems = self.panair_folder_path.split(os.sep)
        self.panair_folder_path_unix = '/'.join(path_elems)
        self.panair_exec_path_unix = '/'.join([self.panair_folder_path_unix, 'panair_exec.sh'])

        ### Title ###
        wgs = wgs_creator.LaWGS(self.title)

        ### Airfoils ###

        if all(airfoil[:4] =='naca' and len(airfoil) == 8 for airfoil in airfoils.values()):
            '''
            Currently only used for test case: simple wing with only root and tip airfoil.
            '''
            for airfoil in airfoils.values():
                root_airfoil_wgs = wgs_creator.naca4digit(airfoil[4:], num = 25, chord=self.root_chord, y_coordinate=0)
                tip_airfoil_wgs  = wgs_creator.naca4digit(airfoil[4:], num = 25, chord=self.tip_chord, y_coordinate=self.half_span)
        else:
            '''
            Airfoils are read from csv files with their xz coordinates, i.e. their size and rotation.
            Scaling, translation, and rotation must thus be done beforehand and written to the csv file to be read. 
            '''
            transform_airfoils(airfoils=airfoils,
                               airfoil_folder_path=airfoil_folder_path,
                               root_chord=self.root_chord,
                               tip_chord=self.tip_chord,
                               sweep=self.sweep,
                               half_span=self.half_span,
                               tc_ratio=self.tc_ratio,
                               twist=self.twist)

            root_airfoil_csv = os.path.join(airfoil_folder_path, 'root_airfoil.csv')
            tip_airfoil_csv = os.path.join(airfoil_folder_path, 'tip_airfoil.csv')

            root_airfoil_wgs = wgs_creator.read_airfoil(root_airfoil_csv, y_coordinate=0.)
            tip_airfoil_wgs = wgs_creator.read_airfoil(tip_airfoil_csv, y_coordinate=self.half_span)

        ### Wing ###
        wing = root_airfoil_wgs.linspace(tip_airfoil_wgs, num=20)
        wgs.append_network("wing", wing, boun_type=1)

        ### Create wingtip surface ###
        wingtip_upper, wingtip_lower = tip_airfoil_wgs.split_half()
        wingtip_lower = wingtip_lower.flip()
        wingtip = wingtip_upper.linspace(wingtip_lower, num=5)
        wgs.append_network("wingtip", wingtip, 1)

        # plot_wireframe(wing)
        # plot_wireframe(wingtip)

        ### Make wing wake ###
        wingwake = wing.make_wake(edge_number=3, wake_length=50 * self.root_chord)
        wgs.append_network("wingwake", wingwake, 18)

        wgs.create_wgs(filename=os.path.join(self.panair_folder_path, f'{self.title}.wgs'))

        self.wgs = wgs
        self.wing = wing
    def evaluate_AoA(self, AoA: Union[float,int]):
        """
        Method to evaluate aerodynamic coefficients at a given angle of attack using PanAir.
        :param AoA: Angle of attack to be evaluated
        :param output:
        :return:
        """
        ### Clean up temporary files, must be done for new run ###
        cleanup_command = os.path.join(self.panair_folder_path, 'clean502.bat')
        os.system(cleanup_command)

        self.wgs.create_aux(alpha=AoA, mach=self.M, cbar=self.MAC, span=self.b,
                            sref=self.Sref, xref=self.root_chord / 4, zref=0.,
                            filename=os.path.join(self.panair_folder_path,
                                                  f'{self.title}.aux'))

        run_panair_command = f"{self.cygwin_path} --login -c '{self.panair_exec_path_unix}" \
                             f" {self.panair_folder_path_unix} {self.title}.aux'"
        os.system(run_panair_command)

        #time to write files before reading
        sleep(0.1)

        ### Read results based on output settings###
        results_file = os.path.join(self.panair_folder_path, 'ffmf')
        self.results_df = read_output(results_file)
        return float(self.results_df['cl'][0])
    def find_cruise_AoA(self):
        def opt_func(AoA):
            cl = self.evaluate_AoA(AoA)
            diff = cl - self.cl
            #print(f'AoA of {AoA} gives cl of {cl}. Difference is {diff}.')
            return diff

        #AoA_cruise = root(opt_func, 0, tol=1e-4)
        AoA_cruise = root_scalar(opt_func, x0=-1, x1=1, method='secant', rtol=1e-5)
        return AoA_cruise.root

    def convert_coefficients_to_forces(self):
        cl = float(self.results_df['cl'][0])
        cd = float(self.results_df['cdi'][0])

        L = cl*0.5*self.rho*self.V**2*self.Sref
        D = cd*0.5*self.rho*self.V**2*self.Sref
        return L,D

    def plot_wing(self):
        plot_wireframe(self.wing)
def plot_wireframe(wing, show_corners=True, show_edges=True, show_normvec=True):
    """ plot the Network as a wireframe. Adapted from wgs_creator to resolve versioning issues.
    :param show_corners: display the corner numbers of the Network
    :param show_edges: display the edge numbers of the Network
    :param show_normvec: show a vector pointing out from the front side of the Network
                        (will not work for skewed Networks and Networks consisting only a pair of Lines)"""
    fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
    ax.set_aspect("equal")
    x, y, z = (wing[:, :, i] for i in range(3))
    ax.plot_wireframe(x, y, z)
    # place invisible markers to set the aspect ratio of each axis to be equal
    max_range = np.array([x.max() - x.min(), y.max() - y.min(), z.max() - z.min()]).max()
    xb = 0.5 * max_range * np.mgrid[-1:2:2, -1:2:2, -1:2:2][0].ravel() + 0.5 * (x.max() + x.min())
    yb = 0.5 * max_range * np.mgrid[-1:2:2, -1:2:2, -1:2:2][1].ravel() + 0.5 * (y.max() + y.min())
    zb = 0.5 * max_range * np.mgrid[-1:2:2, -1:2:2, -1:2:2][2].ravel() + 0.5 * (z.max() + z.min())
    for xxb, yyb, zzb in zip(xb, yb, zb):
        ax.plot([xxb], [yyb], [zzb], 'w')
    ax.set_xlabel("$x$", fontsize=16)
    ax.set_ylabel("$y$", fontsize=16)
    ax.set_zlabel("$z$", fontsize=16)
    rowmid, colmid = (i // 2 for i in wing.shape[:2])
    if show_corners:
        ax.text(*wing[0, 0, :], s="1")
        ax.text(*wing[-1, 0, :], s="2")
        ax.text(*wing[-1, -1, :], s="3")
        ax.text(*wing[0, -1, :], s="4")
    if show_edges:
        ax.text(*wing[rowmid, 0, :], s="edge1")
        ax.text(*wing[-1, colmid, :], s="edge2")
        ax.text(*wing[rowmid, -1, :], s="edge3")
        ax.text(*wing[0, colmid, :], s="edge4")
    if show_normvec:
        v_edge4 = wing[rowmid, colmid, :] - wing[rowmid, colmid - 1, :]
        v_edge1 = wing[rowmid, colmid, :] - wing[rowmid - 1, colmid, :]
        normvec = np.cross(v_edge1, v_edge4)
        normvec = normvec / np.linalg.norm(normvec) * max_range / 10.
        a = ([pos, vec + pos] for (vec, pos) in zip(normvec, wing[rowmid, colmid, :]))
        # a = Arrow3D(*a, color="k", mutation_scale=20, arrowstyle="-|>")
        # ax.add_artist(a)
    plt.show()
    return

def read_input_from_XML(in_file, HF_tool_info):
    input_obj = DataHandler.create_objectify_from_XML_file(in_file)

    vardict = dict()
    for var_path in HF_tool_info['inputs']:
        var = input_obj.xpath(var_path)[0]
        varkey = f'var_{var.tag}'
        varval = float(var.text)
        vardict.update({varkey:varval})
    return vardict

def write_output_to_XML(outputdict: dict, out_file, HF_tool_info):
    output_obj = DataHandler.create_objectify_from_XML_file(out_file)

    for var_path in HF_tool_info['outputs']:
        var = output_obj.xpath(var_path)[0]
        val = outputdict[var.tag]
        var._setText(str(val))

    with open(out_file, 'wb') as doc:
        doc.write(tostring(output_obj, xml_declaration=True, encoding='utf-8', pretty_print=True))

def run_tool(in_file, out_file):
    HF_tool_info = dict(
        name="Panair",
        description="Panel code solver for sub- and supersonic flow",
        inputs=[
            '/dataSchema/aircraft/geometry/lambda',
            '/dataSchema/aircraft/geometry/Lambda',
            '/dataSchema/aircraft/geometry/Sref',
            '/dataSchema/aircraft/geometry/tc',
            '/dataSchema/aircraft/geometry/AR',
            '/dataSchema/aircraft/geometry/Theta',
            '/dataSchema/aircraft/weight/WT',
            '/dataSchema/reference/h',
            '/dataSchema/reference/M',
        ],
        outputs=[
            '/dataSchema/aircraft/other/L',
            '/dataSchema/aircraft/other/D',
            '/dataSchema/aircraft/other/fin'
        ],
        version=0.1,
        tool_to_replace='Aerodynamics'
    )

    input_dict = read_input_from_XML(in_file, HF_tool_info)

    start = time.time()

    ### Airfoils ###
    root_airfoil = 'naca64206-il'
    kink_airfoil = 'naca64206-il'
    tip_airfoil = 'naca64206-il'


    airfoils = dict(
        root_airfoil=root_airfoil,
        kink_airfoil=kink_airfoil,
        tip_airfoil=tip_airfoil
    )

    ### Title and paths ###
    title = 'SSBJ_wing'
    panair_folder_path = os.getcwd()
    cygwin_path = r"C:\cygwin64\bin\bash"

    path_title_dict = dict(
        title=title,
        panair_folder_path=panair_folder_path,
        cygwin_path=cygwin_path
    )

    ### Generate wing parameters ###
    SSBJ_wing = WingGeometryConditions(**input_dict)

    SSBJ_wing.prepare_input(airfoils=airfoils,
                            path_title_dict=path_title_dict)

    SSBJ_wing.plot_wing()

    AoA_cruise = SSBJ_wing.find_cruise_AoA()

    ### Obtain results from stored dataframe (last run done in root finding) ###
    L, D = SSBJ_wing.convert_coefficients_to_forces()

    print(f'Area of single wing: {(SSBJ_wing.root_chord + SSBJ_wing.tip_chord) * SSBJ_wing.half_span / 2} \n'
          f'Angle of attack at cruise: {AoA_cruise} \n'
          f'Lift: {L} \n'
          f'Drag: {D} \n'
          f'L/D: {L / D}')

    end = time.time()
    print(f'Elapsed time: {end - start}')

    outputdict = dict(
        L=L,
        D=D,
        fin=L / D
    )
    write_output_to_XML(outputdict, out_file, HF_tool_info)
    return



if __name__ == '__main__':

    if len(sys.argv) > 1:
        parser = argparse.ArgumentParser(description="next design step tool")
        parser.add_argument('-i', '--infile',
                            help="first input XML file")
        parser.add_argument('-o', '--outfile',
                            help="output XML file")
        args = parser.parse_args()
        in_file = args.infile
        out_file = args.outfile
    else:
        in_file = os.path.join(os.getcwd(), 'ToolInput', 'ToolInput.xml')
        out_file = os.path.join(os.getcwd(), 'ToolOutput', 'ToolOutput.xml')

    run_tool(in_file, out_file)

