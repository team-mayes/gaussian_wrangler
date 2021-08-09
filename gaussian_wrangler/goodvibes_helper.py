# !/usr/bin/env python


"""
Uses GoodVibes to check input and calculate G at a range of temperatures
"""

import argparse
import os
import re
import sys
import numpy as np
import gaussian_wrangler.goodvibes_hm
import jpype
import jpype.imports
from pathlib import Path
from collections import defaultdict
from common_wrangler.common import (InvalidDataError, warning, RG, KB, PLANCK_CONST_JS, EHPART_TO_KCAL_MOL,
                                    GOOD_RET, INPUT_ERROR, IO_ERROR, INVALID_DATA,
                                    write_csv, create_out_fname, make_fig, parse_stoich, capture_stdout, list_to_file,
                                    round_sig_figs)
from gaussian_wrangler import __version__

__author__ = 'hmayes'


# Constants #

# Config keys
DEF_OUT_FILE_NAME = 'aea_out.csv'

# For data processing; standard hartree fieldnames below
# '"File Name","Solvent type","Stoichiometry","Charge","Mult","Functional","Basis Set","Energy (A.U.)","dipole",
# "ZPE (Hartrees)","H298 (Hartrees)","G298 (Hartrees)","Freq 1","Freq 2","BSSE (Hartrees)"'
FUNCTIONAL = 'Functional'
BASIS_SET = 'Basis Set'
SOLV = 'Solvent type'
STOICH = 'Stoichiometry'
FREQS = 'Frequencies 1 and 2'
GAUSS_VER_PAT = re.compile(r"Gaussian.*:.*Rev.*")
GOODVIBES_ERROR_PAT = re.compile(r"x .*")
GOODVIBES_DATA_PAT = re.compile(r"Structure .*")
REACT_PROD_SEP = 'TS'

HARTREE_OUT = 'hartree_output'
GOODVIBES_OUT = 'goodvibes_output'

# for printing
FILE1 = 'file1'
FILE2 = 'file2'
FILE3 = 'file3'
FILE4 = 'file4'
FILE5 = 'file5'
A = 'A (1/s if uni)'
EA = 'Ea (kcal/mol)'
DELTA_G_TEMP = '\u0394G temp (K)'
RATE_COEFF_AT_G_TEMP = 'Rate coefficient (k) at \u0394G temp (1/s if unimolecular)'
DELTA_G_TS = '\u0394G\u2021 (kcal/mol)'
DELTA_G_RXN = '\u0394G_rxn (kcal/mol)'
QH_A = 'qh_A (1/s if uni)'
QH_EA = 'qh_Ea (kcal/mol)'
QH_RATE_COEFF_AT_G_TEMP = 'qh_Rate coefficient (k) at \u0394G temp (1/s if unimolecular)'
QH_DELTA_G_TS = 'qh_\u0394G\u2021 (kcal/mol)'
QH_DELTA_G_RXN = 'qh_\u0394G_rxn (kcal/mol)'

OUTPUT_HEADERS = [FILE1, FILE2, FILE3, FILE4, FILE5, A, EA, DELTA_G_TEMP, RATE_COEFF_AT_G_TEMP, DELTA_G_TS, DELTA_G_RXN,
                  QH_A, QH_EA, QH_RATE_COEFF_AT_G_TEMP, QH_DELTA_G_TS, QH_DELTA_G_RXN]


class HartreeWrapper:
    def __init__(self):
        jar_path = Path(__file__).parent / "hartree"
        found_jars = sorted(jar_path.glob('*.jar'))
        if len(found_jars) == 0:
            raise Exception("Could not find any JARs in dir {}".format(jar_path))

        for found_jar in found_jars:
            jpype.addClassPath(found_jar)

        print("Classpath: ", jpype.getClassPath())

        jpype.startJVM(convertStrings=False)

        # TODO: Figure out how to scope these imports for the class
        # noinspection PyUnresolvedReferences
        from org.cmayes.hartree.loader.gaussian import SnapshotLoader

        self.loader = SnapshotLoader()

    def read_all_gaussian(self, files):
        # noinspection PyUnresolvedReferences
        from java.io import FileReader
        mapped_results = {}
        for cur_file in files:
            mapped_results[cur_file] = self.loader.load(cur_file, FileReader(cur_file))

        return mapped_results

    def read_gaussian(self, tgt_file):
        # noinspection PyUnresolvedReferences
        from java.io import FileReader
        return self.loader.load(tgt_file, FileReader(tgt_file))

    def __del__(self):
        if jpype.isJVMStarted():
            jpype.shutdownJVM()


# Jpype can't restart the JVM, so we make this global.
hartree = HartreeWrapper()


def parse_cmdline(argv):
    """
    Returns the parsed argument list and return code.
    `argv` is a list of arguments, or `None` for ``sys.argv[1:]``.
    """
    if argv is None:
        argv = sys.argv[1:]

    # initialize the parser object:
    parser = argparse.ArgumentParser(description='Calculates A and Ea from Gaussian output files using GoodVibes. '
                                                 'List files to be analyzed, reactant(s) first and ending with the '
                                                 'transition structure. These can be listed on the command line or in '
                                                 'a file (each line listing a set of reactant(s) and transition '
                                                 'structure).')
    parser.add_argument("-d", "--out_dir", help="A directory where output files should be saved. The default location "
                                                "is the current working directory.", default=None)
    parser.add_argument("-f", dest="freq_cutoff", help="Cut-off frequency for both entropy and enthalpy (wavenumbers) "
                                                       "(default = 0)", default="0")
    parser.add_argument("-l", "--list", help="The location of the list of Gaussian output files. "
                                             "The default file name.", default=None)
    parser.add_argument("-q", "--quasiharmonic", help="Use the '-q' option in GoodVibes, which turns on turns on "
                                                      "quasi-harmonic corrections to both entropy and enthalpy in the "
                                                      "Gibbs free energy (qh-G(T)) output from GoodVibes. ",
                        action='store_true')
    parser.add_argument("--temp", help="Temperature in K for calculating \u0394G. The default is the first "
                                       "temperature in 'temp_range' (if specified). If a value is given, the program "
                                       "will use the temperature closest to it in the temp_range.", default=None)
    parser.add_argument("-ti", "--temp_range", help="Initial temp, final temp, (and optionally) step size (K) for "
                                                    "thermochemistry calculations. The default range is 300,600,30",
                        default="300,600,30")
    parser.add_argument("-v", "--vib_scale", help="Scaling factor to be used for vibrational frequencies. If not "
                                                  "provided, the GoodVibes default value will be used.",
                        default=None)
    parser.add_argument("-p", "--plot", help="Make a \u0394G plot at the specified temp. The default is False.",
                        action='store_true')
    parser.add_argument("-pl", "--plot_labels", help="Optional labels for \u0394G plot. Enter as a list.",
                        default=None)
    parser.add_argument("-c", "--vibes_check", help="In addition to standard checks always run (matching solvent, "
                                                    "level of theory, stoichiometry, charge, multiplicity, and "
                                                    "Gaussian versions), run files through GoodVibes '--check' before "
                                                    "performing calculations. The default is False.",
                        action='store_true')
    parser.add_argument("-o", "--output_fname", help="The name of the output file to be created. The default is the "
                                                     "list name with the extension '.csv', or '{}' if no list name "
                                                     "provided.".format(DEF_OUT_FILE_NAME), default=None)

    parser.add_argument("-s", "--save_vibes", help="Save the output from running GoodVibes in separate files, "
                                                   "named with the Gaussian log file prefix and '.dat'. "
                                                   "The default is False.",
                        action='store_true')
    parser.add_argument("-t", "--tog_vibes", help="Save the output from running GoodVibes in one file, "
                                                  "renamed with the output file prefix and '.dat'. "
                                                  "The default is False.",
                        action='store_true')

    args = None
    try:
        args = parser.parse_known_args(argv)
        options = args[0]
        if not options.out_dir:
            options.out_dir = os.getcwd()
        # user can define a new directory as the output directory
        if not os.path.exists(options.out_dir):
            os.makedirs(options.out_dir)

        if options.output_fname:
            options.output_fname = os.path.abspath(os.path.join(options.out_dir, options.output_fname))
        elif options.list:
            options.output_fname = create_out_fname(options.list, ext='.csv', base_dir=options.out_dir)
        else:
            options.output_fname = create_out_fname(DEF_OUT_FILE_NAME, ext='.csv', base_dir=options.out_dir)

        if options.plot_labels:
            options.plot_labels = options.plot_labels.split(',')
        else:
            options.plot_labels = ['']

        if options.vib_scale:
            options.vib_scale = float(options.vib_scale)

    except (SystemExit, ValueError) as e:
        if hasattr(e, 'code') and e.code == 0:
            return args, GOOD_RET
        warning(e)
        parser.print_help()
        return args, INPUT_ERROR

    return args, GOOD_RET


def get_gauss_results(options, unique_fnames):
    """
    Run hartree and goodvibes only once per file name
    :param options: user-specified options, used here to determine goodvibes input
    :param unique_fnames: a set of unique file names (really, file locations)
    :return: results_dict: dictionary of results from running hartree and goodvibes
    """
    results_dict = defaultdict(dict)
    for fname in unique_fnames:
        if fname != REACT_PROD_SEP:
            base_name = os.path.basename(fname)
            gauss_results = hartree.read_gaussian(fname)
            solvent = gauss_results.getSolvent()

            results_dict[base_name][HARTREE_OUT] = gauss_results
            results_dict[base_name][SOLV] = solvent
            results_dict[base_name][FREQS] = gauss_results.getFrequencyValues()
            # later, a regex will be performed on STOICH, and it will expect a standard string, not a java.lang.String
            results_dict[base_name][STOICH] = str(gauss_results.getStoichiometry())
            vibes_input = [fname, "--ti", options.temp_range, "-f", options.freq_cutoff]
            if solvent:
                vibes_input += ["-c", "1"]
            if options.quasiharmonic:
                vibes_input += ["-q"]
            if options.vib_scale:
                vibes_input += ["-v", str(options.vib_scale)]
            with capture_stdout(gaussian_wrangler.goodvibes_hm.main, vibes_input) as output:
                results_dict[base_name][GOODVIBES_OUT] = output.split('\n')
    return results_dict


def check_gausslog_fileset(file_set, good_vibes_check, results_dict):
    """
    checks include:
       using hartree to get info to check for:
           the correct number of imaginary freq
           the stoichiometry adds up
           same implicit solvent (or lack thereof)
           same functional and basis set
        awk for same versions of Gaussian
        made GoodVibes checks optional to save run time
    :param file_set: list of reactant file(s) and TS file
    :param good_vibes_check: boolean to run goodvibes checking; will slow down calculations
    :param results_dict: dictionary of results from running hartree and goodvibes
    :return: reaction_type: integer for molecularity of reaction
    """
    total_react_charge = 0
    ts_charge = np.nan
    ts_index = None
    total_product_charge = 0
    multiplicities = np.full([len(file_set)], np.nan)

    # initialize the variables below to make my IDE happy
    solvent = None
    func = None
    basis = None
    gauss_ver = None
    react_stoich_dict = {}  # empty dict to make IDE happy
    ts_stoich_dict = {}  # empty dict to make IDE happy
    prod_stoich_dict = {}  # empty dict to make IDE happy
    reading_reactants = True

    for index, fname in enumerate(file_set):
        base_name = os.path.basename(fname)
        if fname == REACT_PROD_SEP:
            ts_index = index
            reading_reactants = False
            # make it so that later multiplicity check works
            if index == 0:
                raise InvalidDataError("Did not expect to read find reactant/product separator ('{}') as the first "
                                       "entry".format(REACT_PROD_SEP))
            else:
                multiplicities[index] = multiplicities[index - 1]
            continue

        # now start checks by getting info from hartree
        gauss_result = results_dict[base_name][HARTREE_OUT]
        freq_vals = results_dict[base_name][FREQS]
        stoich = results_dict[base_name][STOICH]

        # exit effort if there files with more than one imaginary frequency
        if freq_vals[0] < 0 and freq_vals[1] < 0:
            raise InvalidDataError("The first two frequencies are both imaginary in file: {}".format(fname))
        # First see if it is the TS
        if freq_vals[0] < 0:
            if not reading_reactants:
                raise InvalidDataError("In each set of files, only one file with an imaginary frequency is expected\n"
                                       "    (or none if reactant/product separator is used). Unexpectedly found an "
                                       "imaginary frequency in\n    file: {}\n    after already finding the TS to be: "
                                       "{}".format(fname, file_set[ts_index]))
            reading_reactants = False
            ts_index = index
            ts_charge = gauss_result.getCharge()
            ts_stoich_dict = parse_stoich(stoich)
        elif reading_reactants:
            total_react_charge += gauss_result.getCharge()
            if len(react_stoich_dict) == 0:
                react_stoich_dict = parse_stoich(stoich)
            else:
                react_stoich_dict = parse_stoich(stoich, add_to_dict=react_stoich_dict)
        else:
            total_product_charge += gauss_result.getCharge()
            if len(prod_stoich_dict) == 0:
                prod_stoich_dict = parse_stoich(stoich)
            else:
                prod_stoich_dict = parse_stoich(stoich, add_to_dict=prod_stoich_dict)

        # additional checks on all files as we go...
        multiplicities[index] = int(gauss_result.getMult())
        file_gauss_ver = None
        i = 0
        with open(fname) as f:
            for line in f:
                s_line = line.strip()
                if GAUSS_VER_PAT.match(s_line):
                    file_gauss_ver = s_line.split()[:3]
                    break
                i += 1
                # just in case not caught... don't read the whole file
                if i > 160:
                    break
        if index == 0:
            # make all lower case to remove chance of flagging this insignificant difference
            solvent = str(results_dict[base_name][SOLV]).lower()
            func = str(gauss_result.getFunctional()).lower()
            # ignore differences between restricted and unrestricted versions of the functional
            if func.startswith("r") or func.startswith("u"):
                func = func[1:]
            basis = str(gauss_result.getBasisSet()).lower()
            gauss_ver = file_gauss_ver
        else:
            if str(gauss_result.getSolvent()).lower() != solvent:
                raise InvalidDataError("Different solvents ({}, {}) found for file set: "
                                       "{}".format(solvent, gauss_result.getSolvent(), file_set))
            # ignore differences between restricted and unrestricted versions of the functional
            current_func = str(gauss_result.getFunctional()).lower()
            if current_func.startswith("u") or current_func.startswith("r"):
                current_func = current_func[1:]
            if current_func != func:
                raise InvalidDataError("Different functionals ({}, {}) found for file set: "
                                       "{}".format(func, gauss_result.getFunctional(), file_set))
            if str(gauss_result.getBasisSet()).lower() != basis:
                raise InvalidDataError("Different basis sets ({}, {}) found for file set: "
                                       "{}".format(basis, gauss_result.getBasisSet(), file_set))
            if gauss_ver != file_gauss_ver:
                warning("Different Gaussian versions ({}, {}) found for file set: {}".
                        format(gauss_ver, file_gauss_ver, file_set))

    # Now overall checks
    file_set_str = "\n              ".join([""] + [os.path.relpath(f) for f in file_set])
    if len(ts_stoich_dict) > 0:
        if react_stoich_dict != ts_stoich_dict:
            raise InvalidDataError("Check stoichiometries of reactant(s) and transition state for set: {}\n"
                                   "reactants: {}, products: {}".format(file_set_str, react_stoich_dict,
                                                                        ts_stoich_dict))
        if total_react_charge != ts_charge:
            raise InvalidDataError("Check charge of reactant(s) and transition state for set: {}\n"
                                   "Found {} and {}, respectively".format(file_set_str, total_react_charge, ts_charge))
    if len(prod_stoich_dict) > 0:
        if react_stoich_dict != prod_stoich_dict:
            raise InvalidDataError("Check stoichiometries of reactant(s) and product(s) for set: {}\n"
                                   "reactants: {}, products: {}".format(file_set_str, react_stoich_dict,
                                                                        prod_stoich_dict))
        if total_react_charge != total_product_charge:
            raise InvalidDataError("Check charge of reactant(s) and product(s) for set: {}\nFound {} and {}, "
                                   "respectively".format(file_set_str, total_react_charge, total_product_charge))

    mult_check = np.sum(multiplicities - multiplicities[0])
    if mult_check != 0:
        raise InvalidDataError("Check multiplicities in set: {}\nFound: {}".format(file_set, multiplicities))

    if good_vibes_check:
        if solvent:
            file_set = file_set + ["-c", "1"]
        if REACT_PROD_SEP in file_set:
            file_set.remove(REACT_PROD_SEP)
        with capture_stdout(gaussian_wrangler.goodvibes_hm.main, file_set + ["--check"]) as output:
            vibes_out = output.split("\n")
        for line in vibes_out:
            if GOODVIBES_ERROR_PAT.match(line):
                if 'Different charge and multiplicity' in line:
                    # already checked this, so we can skip this error
                    continue
                raise InvalidDataError("See GoodVibes error checking report: 'Goodvibes_output.dat'")

    return solvent, ts_index


def get_thermochem(file_set, results_dict, save_vibes, out_dir, tog_output_fname, qh_h_opt, write_mode):
    """
    Calls GoodVibes to get thermochem at a range of temps
    :param file_set: list of reactant file(s), TS file (or separator), and optionally products
    :param results_dict: dictionary of results from running hartree and goodvibes
    :param save_vibes: boolean to determine whether to save each GoodVibes output separately
    :param out_dir: directory to save GoodVibes output files (if requested)
    :param tog_output_fname: None or string (file name) if saving each GoodVibes output together
    :param qh_h_opt: boolean to use the '-q' option in GoodVibes (corrections to both entropy and enthalpy)
    :param write_mode: boolean to start a new to add to an all-together goodvibes output file
    :return: nothing
    """
    h = []
    qh_h = []
    gt = []
    qh_gt = []
    temps = []
    for index, file in enumerate(file_set):
        base_name = os.path.basename(file)
        if file == REACT_PROD_SEP:
            h.append(np.full([len(temps)], np.nan))
            qh_h.append(np.full([len(temps)], np.nan))
            gt.append(np.full([len(temps)], np.nan))
            qh_gt.append(np.full([len(temps)], np.nan))
            continue
        vibes_out = results_dict[base_name][GOODVIBES_OUT]
        found_structure = False
        skip_line = True
        h.append([])
        qh_h.append([])
        gt.append([])
        qh_gt.append([])
        # we know the last line should be dropped, and at least the first 10
        for line in vibes_out[10:-2]:
            if GOODVIBES_ERROR_PAT.match(line):
                raise InvalidDataError("See GoodVibes output: {}".format(vibes_out))
            if not found_structure:
                if GOODVIBES_DATA_PAT.match(line):
                    found_structure = True
                    continue
            elif skip_line:
                skip_line = False
                continue
            else:
                vals = line.split()
                if index == 0:
                    temps.append(float(vals[1]))
                h[index].append(float(vals[2]))
                if qh_h_opt:
                    qh_h[index].append(float(vals[3]))
                gt[index].append(float(vals[-2]))
                qh_gt[index].append(float(vals[-1]))
        if save_vibes:
            vibes_out_fname = os.path.relpath(create_out_fname(file, suffix='_vibes', base_dir=out_dir, ext='.dat'))
            list_to_file(vibes_out, vibes_out_fname, print_message=False)
            print('Saved GoodVibes output as: {}'.format(vibes_out_fname))
        if tog_output_fname:
            list_to_file(vibes_out, tog_output_fname, mode=write_mode, print_message=False)
            if write_mode == 'w':
                print("Adding all GoodVibes output to: {}".format(tog_output_fname))
                write_mode = "a"

    temps = np.asarray(temps)
    # for each molecule, multiply the array to convert to kcal/mol
    for index in range(len(gt)):
        h[index] = np.asarray(h[index]) * EHPART_TO_KCAL_MOL
        if qh_h_opt:
            qh_h[index] = np.asarray(qh_h[index]) * EHPART_TO_KCAL_MOL
        gt[index] = np.asarray(gt[index]) * EHPART_TO_KCAL_MOL
        qh_gt[index] = np.asarray(qh_gt[index]) * EHPART_TO_KCAL_MOL

    return temps, h, qh_h, gt, qh_gt


def get_deltas(temps, vals, ts_index):
    """
    Calculate the difference in values (e.g. Gibbs free energy) at each temperature in temps
    :param temps: np array of temperatures (in K) to be evaluated
    :param vals: list of np arrays with quasi-harmonic treatment of G(T) for each output file at each temp
    :param ts_index: int, index of array with values for ts
    :return: delta_ts, delta_rxn: np arrays of delta_gibbs
    """
    nan_array = np.full([len(temps)], np.nan)
    vals_react = nan_array
    vals_ts = nan_array
    vals_prod = nan_array
    # any files before the ts will be a reactant file
    for index in range(len(vals)):
        if index == 0:
            vals_react = vals[index]
        elif index < ts_index:
            vals_react += vals[index]
        elif index == ts_index:
            vals_ts = vals[index]
        elif index == ts_index + 1:
            vals_prod = vals[index]
        else:
            vals_prod += vals[index]
    delta_ts = vals_ts - vals_react
    delta_rxn = vals_prod - vals_react
    return delta_ts, delta_rxn


def get_kt(temps, delta_gibbs_ts):
    """
    Calculate the rate coefficient at each temperature in temps
    :param temps: np array of temperatures (in K) to be evaluated
    :param delta_gibbs_ts: np array of temperatures (in K) to be evaluated
    :return: kt: np array of rate coefficients at each temp
    """
    # rate coefficient from Eyring equation
    return KB / PLANCK_CONST_JS * temps * np.exp(-delta_gibbs_ts / RG / temps)  # [1/s] if unimolecular


def fit_arrhenius(temps, kt):
    """
    Fits a straight line of 1/temps vs. ln(kt) to get A (1/s), Ea (kcal/mol)
    :param temps: numpy array of temps in K
    :param kt: numpy array of rate coefficients in 1/s (or appropriate units)
    :return: calcs A (1/s), Ea (kcal/mol)
    """
    inv_temp = 1 / temps
    ln_kt = np.log(kt)
    fit = np.polyfit(inv_temp, ln_kt, 1)
    slope = fit[0]
    intercept = fit[1]
    a = np.exp(intercept)  # [1/s]
    ea = -slope * RG  # [kcal/mol]
    return round_sig_figs(a), round_sig_figs(ea)


def print_results(a, ea, qh_a, qh_ea, g_temp, k_temp, g_ts, g_rxn, qh_k_temp, qh_g_ts, qh_g_rxn, file_set, out_fname,
                  print_mode, print_message=False):
    file_names = []
    for index in range(-5, 0):
        try:
            file_names.append(os.path.basename(file_set[index]))
        except IndexError:
            file_names.append('')
    output_dict = {FILE1: file_names[0], FILE2: file_names[1], FILE3: file_names[2], FILE4: file_names[3],
                   FILE5: file_names[4], A: a, EA: ea, QH_A: qh_a, QH_EA: qh_ea, DELTA_G_TEMP: g_temp,
                   RATE_COEFF_AT_G_TEMP: k_temp, DELTA_G_TS: g_ts, DELTA_G_RXN: g_rxn,
                   QH_RATE_COEFF_AT_G_TEMP: qh_k_temp, QH_DELTA_G_TS: qh_g_ts, QH_DELTA_G_RXN: qh_g_rxn}
    write_csv([output_dict], out_fname, OUTPUT_HEADERS, extrasaction="ignore", mode=print_mode,
              print_message=print_message)


def get_temp_index(temp, temps):
    """
    Find temperature index to use for returning single values of quantities from arrays
    :param temp: float, temp in K to calculate delta values
    :param temps: list of floats, temps for which we have calculated reaction rates, delta G's, etc.
    :return: temp_index: integer, array location of temp to use to return specific values from arrays
    """
    if temp:
        temp_index = (np.abs(temps - float(temp))).argmin()
    else:
        temp_index = 0

    return temp_index


def plot_delta(fname, temp, delta_ts_list, delta_rxn_list, labels, var='G'):
    """
    Makes a plot of delta G at the specified temp
    :param fname: string, to save plot
    :param temp: float, temp at which delta Gs were calculated
    :param delta_ts_list: list of floats
    :param delta_rxn_list: list of floats
    :param labels: list of strings
    :param var: string, which Delta variable is being plotted
    :return: nothing, just save
    """
    max_y_lines = 5
    x_axis = np.array([0, 1, 3, 4, 6, 7])
    y_axis = []
    y_labels = []
    # y_min = np.floor(np.array(g_rxn_list).min())
    # y_max = np.ceil(np.array(g_ts_list).max())
    for index in range(max_y_lines):
        try:
            y_axis.append(np.array([0.0, 0.0, delta_ts_list[index], delta_ts_list[index],
                                    delta_rxn_list[index], delta_rxn_list[index]]))
        except IndexError:
            y_axis.append(None)
        try:
            y_labels.append(labels[index])
        except IndexError:
            y_labels.append(None)

    make_fig(fname, x_axis, y_axis[0],
             x_label='reaction coordinate', y_label='\u0394' + var + ' at {} K (kcal/mol)'.format(temp),
             y1_label=y_labels[0], y2_label=y_labels[1], y3_label=y_labels[2], y4_label=y_labels[3],
             y5_label=y_labels[4], y2_array=y_axis[1], y3_array=y_axis[2], y4_array=y_axis[3], y5_array=y_axis[4],
             ls2='-', ls3='-', ls4='-', ls5='-',
             # y_lima=y_min, y_limb=y_max,
             hide_x=True,
             )


def process_file_set(file_set, options, print_mode, results_dict, tog_fname):
    solvent, ts_index = check_gausslog_fileset(file_set, options.vibes_check, results_dict)
    temps, h, qh_h, gt, qh_gt = get_thermochem(file_set, results_dict, options.save_vibes,
                                               options.out_dir, tog_fname, options.quasiharmonic, print_mode)
    delta_h_ts, delta_h_rxn = get_deltas(temps, h, ts_index)
    if options.quasiharmonic:
        qh_delta_h_ts, qh_delta_h_rxn = get_deltas(temps, qh_h, ts_index)
    else:
        qh_delta_h_ts, qh_delta_h_rxn = 0, 0  # Just to make IDE happy...
    delta_gibbs_ts, delta_gibbs_rxn = get_deltas(temps, gt, ts_index)
    qh_delta_gibbs_ts, qh_delta_gibbs_rxn = get_deltas(temps, qh_gt, ts_index)
    if REACT_PROD_SEP in file_set:
        kt, qh_kt, a, ea, qh_a, qh_ea = '', '', '', '', '', ''
    else:
        kt = get_kt(temps, delta_gibbs_ts)
        qh_kt = get_kt(temps, qh_delta_gibbs_ts)
        a, ea = fit_arrhenius(temps, kt)
        qh_a, qh_ea = fit_arrhenius(temps, qh_kt)
    return (temps, a, ea, kt, delta_h_ts, delta_h_rxn, delta_gibbs_ts, delta_gibbs_rxn,
            qh_a, qh_ea, qh_kt, qh_delta_h_ts, qh_delta_h_rxn, qh_delta_gibbs_ts, qh_delta_gibbs_rxn)


def main(argv=None):
    print(f"Running GaussianWrangler script goodvibes_helper version {__version__}")
    # Read input
    args, ret = parse_cmdline(argv)
    if ret != GOOD_RET or args is None:
        return ret

    try:
        # Make a list of lists; each inner list a set of reactant file(s) with TS
        # Include anything in the "list" file as well as entered on the command line
        options = args[0]
        if options.list:
            with open(options.list) as f:
                row_list = [row.strip().split() for row in f.readlines()]
                row_list = list(filter(None, row_list))
        else:
            row_list = []
        if len(args[1]) > 0:
            row_list.append(args[1])
        if len(row_list) == 0:
            raise InvalidDataError("No files or list of files found")

        # now a quick first check that all files exist, and get unique names
        missing_files = set()
        unique_fnames = set()
        for file_set in row_list:
            for file in file_set:
                if file != REACT_PROD_SEP:
                    if os.path.isfile(file):
                        unique_fnames.add(file)
                    else:
                        missing_files.add(file)
        if len(missing_files) > 0:
            raise IOError(missing_files)

        # Initialization to make IDE happy; used for plotting
        g_ts_list, g_rxn_list, qh_g_ts_list, qh_g_rxn_list = [], [], [], []
        g_temp = None
        h_ts_list, h_rxn_list, qh_h_ts_list, qh_h_rxn_list = [], [], [], []
        # now the calculations and printing
        print_mode = 'w'  # for the AEa output, so only prints header once, and then appends to file
        print_message = True
        if options.tog_vibes:
            tog_fname = os.path.relpath(create_out_fname(options.output_fname, suffix='_vibes', ext='.dat'))
        else:
            tog_fname = None
        results_dict = get_gauss_results(options, unique_fnames)
        for file_set in row_list:
            # the called method returns values needed for printing and plotting
            temps, a, ea, kt, delta_h_ts, delta_h_rxn, delta_gibbs_ts, delta_gibbs_rxn, qh_a, qh_ea, qh_kt, \
                qh_delta_h_ts, qh_delta_h_rxn, qh_delta_gibbs_ts, qh_delta_gibbs_rxn = \
                process_file_set(file_set, options, print_mode, results_dict, tog_fname)

            temp_index = get_temp_index(options.temp, temps)
            if REACT_PROD_SEP in file_set:
                k_temp = ""
                qh_k_temp = ""
            else:
                k_temp = round_sig_figs(kt[temp_index])
                qh_k_temp = round_sig_figs(qh_kt[temp_index])
            g_temp = temps[temp_index]
            g_ts = round_sig_figs(delta_gibbs_ts[temp_index])
            g_rxn = round_sig_figs(delta_gibbs_rxn[temp_index])
            qh_g_ts = round_sig_figs(qh_delta_gibbs_ts[temp_index])
            qh_g_rxn = round_sig_figs(qh_delta_gibbs_rxn[temp_index])
            h_ts = round_sig_figs(delta_h_ts[temp_index])
            h_rxn = round_sig_figs(delta_h_rxn[temp_index])
            if options.quasiharmonic:
                qh_h_ts = round_sig_figs(qh_delta_h_ts[temp_index])
                qh_h_rxn = round_sig_figs(qh_delta_h_rxn[temp_index])
            else:
                qh_h_ts, qh_h_rxn = 0, 0  # So don't use an undefined variable below

            print_results(a, ea, qh_a, qh_ea, g_temp, k_temp, g_ts, g_rxn, qh_k_temp, qh_g_ts, qh_g_rxn,
                          file_set, options.output_fname, print_mode, print_message=print_message)
            if options.plot:
                g_ts_list.append(g_ts)
                g_rxn_list.append(g_rxn)
                qh_g_ts_list.append(qh_g_ts)
                qh_g_rxn_list.append(qh_g_rxn)
                h_ts_list.append(h_ts)
                h_rxn_list.append(h_rxn)
                if options.quasiharmonic:
                    qh_h_ts_list.append(qh_h_ts)
                    qh_h_rxn_list.append(qh_h_rxn)

            print_mode = 'a'
            print_message = False

        if options.plot:
            g_fname = create_out_fname(options.output_fname, suffix='_g', ext='.png')
            plot_delta(g_fname, g_temp, g_ts_list, g_rxn_list, options.plot_labels)
            qh_g_fname = create_out_fname(options.output_fname, suffix='_g_qh', ext='.png')
            plot_delta(qh_g_fname, g_temp, qh_g_ts_list, qh_g_rxn_list, options.plot_labels)
            h_fname = create_out_fname(options.output_fname, suffix='_h', ext='.png')
            plot_delta(h_fname, g_temp, h_ts_list, h_rxn_list, options.plot_labels, var='H')
            if options.quasiharmonic:
                qh_h_fname = create_out_fname(options.output_fname, suffix='_h_qh', ext='.png')
                plot_delta(qh_h_fname, g_temp, qh_h_ts_list, qh_h_rxn_list, options.plot_labels, var='H')

    except IOError as e:
        warning("Problems reading file:", e)
        return IO_ERROR
    except InvalidDataError as e:
        warning("Problems reading data:", e)
        return INVALID_DATA

    return GOOD_RET  # success


if __name__ == '__main__':
    status = main()
    sys.exit(status)
