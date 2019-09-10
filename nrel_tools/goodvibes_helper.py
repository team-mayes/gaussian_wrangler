#!/usr/bin/env python
"""
Uses GoodVibes to check input and calculate G at a range of temperatures
"""

from __future__ import print_function

import csv
import os
import subprocess
import sys
import argparse
import re
import numpy as np
# from goodvibes import GoodVibes
from nrel_tools.common import (InvalidDataError, warning, RG, KB, H, EHPART_TO_KCALMOL,
                               GOOD_RET, INPUT_ERROR, IO_ERROR, INVALID_DATA,
                               write_csv, silent_remove, create_out_fname)

try:
    # noinspection PyCompatibility
    from ConfigParser import ConfigParser, MissingSectionHeaderError
except ImportError:
    # noinspection PyCompatibility
    from configparser import ConfigParser, MissingSectionHeaderError

__author__ = 'hmayes'

# Constants #


# Config keys
DEF_HARTREE_LOC = '/Users/hmayes/.local/bin/hartree-cli-1.2.4.jar'
DEF_OUT_FILE_NAME = 'aea_out.csv'

# For data processing; standard hartree fieldnames below
# '"File Name","Solvent type","Stoichiometry","Charge","Mult","Functional","Basis Set","Energy (A.U.)","dipole",
# "ZPE (Hartrees)","H298 (Hartrees)","G298 (Hartrees)","Freq 1","Freq 2","BSSE (Hartrees)"'
FUNCTIONAL = 'Functional'
BASIS_SET = 'Basis Set'
STOICH = 'Stoichiometry'
SOLV = 'Solvent type'
CHARGE = 'Charge'
MULT = 'Mult'
FREQ1 = 'Freq 1'
FREQ2 = 'Freq 2'
AWK_GRAB_GAUSS_VER = ['awk', '/\*\*\*/{getline; print; exit}']
GAUSSIAN_SEPARATOR = '******************************************'
GOODVIBES_OUT_FNAME = "Goodvibes_output.dat"
GOODVIBES_ERROR_PAT = re.compile(r"x .*")
GOODVIBES_DATA_PAT = re.compile(r" {3}Structure .*")
REACT_PROD_SEP = 'TS'

# for printing
FILE1 = 'file1'
FILE2 = 'file2'
FILE3 = 'file3'
FILE4 = 'file4'
FILE5 = 'file5'
A = 'A (1/s if uni)'
EA = 'Ea (kcal/mol)'
DELTA_G_TEMP = '\u0394G temp (K)'
DELTA_G_TS = '\u0394G\u2021 (kcal/mol)'
DELTA_G_RXN = '\u0394G_rxn (kcal/mol)'
QH_A = 'qh_A (1/s if uni)'
QH_EA = 'qh_Ea (kcal/mol)'
QH_DELTA_G_TS = '\u0394G\u2021 (kcal/mol)'
QH_DELTA_G_RXN = '\u0394G_rxn (kcal/mol)'

OUTPUT_HEADERS = [FILE1, FILE2, FILE3, FILE4, FILE5, A, EA, DELTA_G_TEMP, DELTA_G_TS, DELTA_G_RXN,
                  QH_A, QH_EA, QH_DELTA_G_TS, QH_DELTA_G_RXN]


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
    parser.add_argument("-hl", "--hartree_location", help="Optional: the call to invoke hartree. The default is: "
                                                          "{}".format(DEF_HARTREE_LOC), default=DEF_HARTREE_LOC)
    parser.add_argument("-l", "--list", help="The location of the list of Gaussian output files. "
                                             "The default file name.", default=None)
    parser.add_argument("--temp", help="Temperature in K for calculating \u0394G. The default is the first "
                                       "temperature in 'temp_range' (if specified).", default=np.nan)
    parser.add_argument("-ti", "--temp_range", help="Initial temp, final temp, (and optionally) step size (K) for "
                                                    "thermochemistry calculations. The default range is 300,600,30",
                        default="300,600,30")
    parser.add_argument("-o", "--output_fname", help="The name of the output file to be created. The default is {}"
                                                     "".format(DEF_OUT_FILE_NAME), default=DEF_OUT_FILE_NAME)
    parser.add_argument("-q", "--quasiharmonic", help="Use the '-q' option in GoodVibes, which turns on turns on "
                                                      "quasi-harmonic corrections to both entropy and enthalpy in the "
                                                      "Gibbs free energy (qh-G(T)) output from GoodVibes. ",
                        action='store_true')
    parser.add_argument("-s", "--save_vibes", help="Save the output from running GoodVibes in separate files, "
                                                   "renamed with the Gaussian log file prefix and '.dat'. "
                                                   "The default is False.",
                        action='store_true')
    parser.add_argument("-t", "--tog_vibes", help="Save the output from running GoodVibes in one file, "
                                                  "renamed with the output file prefix and '.dat'. "
                                                  "The default is False.",
                        action='store_true')
    parser.add_argument("-v", "--vibes_check", help="In addition to standard checks always run (matching solvent, "
                                                    "level of theory, stoichiometry, charge, multiplicity, and "
                                                    "Gaussian versions), run files through GoodVibes '--check' before "
                                                    "performing calculations. The default is False.",
                        action='store_true')

    args = None
    try:
        args = parser.parse_known_args(argv)
        if not args[0].out_dir:
            args[0].out_dir = os.getcwd()
        # user can define a new directory as the output directory
        if not os.path.exists(args[0].out_dir):
            os.makedirs(args[0].out_dir)
        args[0].output_fname = os.path.abspath(os.path.join(args[0].out_dir, args[0].output_fname))

    except SystemExit as e:
        if hasattr(e, 'code') and e.code == 0:
            return args, GOOD_RET
        warning(e)
        parser.print_help()
        return args, INPUT_ERROR

    return args, GOOD_RET


def check_gausslog_fileset(file_set, hartree_loc, good_vibes_check):
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
    :param hartree_loc: location where hartree can be found to run
    :param good_vibes_check: boolean to run goodvibes checking; will slow down calculations
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
        if fname == REACT_PROD_SEP:
            ts_index = index
            reading_reactants = False
            # make it so that later multiplicity check works
            if index == 0:
                raise InvalidDataError("Did not expect to read find reactant/product separator ('{}') as the first "
                                       "entry".format(REACT_PROD_SEP))
            else:
                multiplicities[index] = multiplicities[index-1]
            continue

        # now start checks by getting info from hartree
        hartree_output = subprocess.check_output(["java", "-jar", hartree_loc,
                                                  "snap", "-f", fname]).decode("utf-8").strip().split("\n")
        hartree_list = list(csv.DictReader(hartree_output, quoting=csv.QUOTE_NONNUMERIC))
        for row in hartree_list:
            # exclude any crazy two imaginary frequency files
            if float(row[FREQ1]) < 0 and float(row[FREQ2]) < 0:
                raise InvalidDataError("The first two frequencies are both imaginary in file: {}".format(fname))
            # First see if it is the TS
            if float(row[FREQ1]) < 0:
                if not reading_reactants:
                    raise InvalidDataError("Expected only one file with an imaginary frequency (or none if "
                                           "reactant/product separator is used). Unexpectedly found an imaginary "
                                           "frequency in file: {}".format(fname))
                reading_reactants = False
                ts_index = index
                ts_charge = int(row[CHARGE])
                ts_stoich_dict = parse_stoich(row[STOICH])
            elif reading_reactants:
                total_react_charge += int(row[CHARGE])
                if len(react_stoich_dict) == 0:
                    react_stoich_dict = parse_stoich(row[STOICH])
                else:
                    react_stoich_dict = parse_stoich(row[STOICH], add_to_dict=react_stoich_dict)
            else:
                total_product_charge += int(row[CHARGE])
                if len(prod_stoich_dict) == 0:
                    prod_stoich_dict = parse_stoich(row[STOICH])
                else:
                    prod_stoich_dict = parse_stoich(row[STOICH], add_to_dict=react_stoich_dict)

            # additional checks on all files as we go...
            multiplicities[index] = int(row[MULT])
            file_gauss_ver = subprocess.check_output(AWK_GRAB_GAUSS_VER + [fname]).strip().split()[:3]
            if index == 0:
                solvent = row[SOLV]
                func = row[FUNCTIONAL]
                basis = row[BASIS_SET]
                gauss_ver = file_gauss_ver
            else:
                if row[SOLV] != solvent:
                    raise InvalidDataError("Different solvents ({}, {}) found for file set: "
                                           "{}".format(solvent, row[SOLV], file_set))
                if row[FUNCTIONAL] != func:
                    raise InvalidDataError("Different functionals ({}, {}) found for file set: "
                                           "{}".format(func, row[FUNCTIONAL], file_set))
                if row[BASIS_SET] != basis:
                    raise InvalidDataError("Different basis sets ({}, {}) found for file set: "
                                           "{}".format(basis, row[BASIS_SET], file_set))
                if gauss_ver != file_gauss_ver:
                    raise InvalidDataError("Different Gaussian versions ({}, {}) found for file set: "
                                           "{}".format(gauss_ver, file_gauss_ver, file_set))

    if len(ts_stoich_dict) > 0:
        if react_stoich_dict != ts_stoich_dict:
            raise InvalidDataError("Check stoichiometries of reactant(s) and transition state for "
                                   "set: {}".format(file_set))
        if total_react_charge != ts_charge:
            raise InvalidDataError("Check charge of reactant(s) and transition state for "
                                   "set: {}".format(file_set))
    if len(prod_stoich_dict) > 0:
        if react_stoich_dict != prod_stoich_dict:
            raise InvalidDataError("Check stoichiometries of reactant(s) and product(s) for set: {}".format(file_set))
        if total_react_charge != total_product_charge:
            raise InvalidDataError("Check charge of reactant(s) and product(s) for set: {}".format(file_set))

    mult_check = np.sum(multiplicities - multiplicities[0])
    if mult_check != 0:
        raise InvalidDataError("Check multiplicities in set: {}".format(file_set))

    if good_vibes_check:
        if solvent:
            file_set = file_set + ["-c", "1"]
        if REACT_PROD_SEP in file_set:
            file_set.remove(REACT_PROD_SEP)
        vibes_out = subprocess.check_output(["python", "-m", "goodvibes"] + file_set +
                                            ["--check"]).decode("utf-8").strip().split("\n")
        for line in vibes_out:
            if GOODVIBES_ERROR_PAT.match(line):
                if 'Different charge and multiplicity' in line:
                    # already checked this, so we can skip this error
                    continue
                raise InvalidDataError("See GoodVibes error checking report: 'Goodvibes_output.dat'")

    return solvent, ts_index


def parse_stoich(stoich_string, add_to_dict=None):
    raw_list = re.findall(r'([A-Z][a-z]*)(\d*)', stoich_string)
    stoich_dict = {}
    for atom_tuple in raw_list:
        if atom_tuple[1] == '':
            stoich_dict[atom_tuple[0]] = 1
        else:
            stoich_dict[atom_tuple[0]] = int(atom_tuple[1])
    if add_to_dict:
        for key, val in add_to_dict.items():
            if key in stoich_dict:
                stoich_dict[key] += add_to_dict[key]
            else:
                stoich_dict[key] = add_to_dict[key]
    return stoich_dict


def get_thermochem(file_set, temp_range, solvent, save_vibes, out_dir, tog_output_fname, qs_h_s):
    """
    Calls GoodVibes to get thermochem at a range of temps
    :param file_set: list of reactant file(s), TS file (or separator), and optionally products
    :param temp_range: string with range of temperatures at which to calculate thermochem
    :param solvent: boolean to decide whether to include
    :param save_vibes: boolean to determine whether to save each GoodVibes output separately
    :param out_dir: directory to save GoodVibes output files (if requested)
    :param tog_output_fname: None or string (file name) if saving each GoodVibes output together
    :param qs_h_s: boolean to use the '-q' option in GoodVibes (corrections to both entropy and enthalpy)
    :return: nothing
    """
    gt = []
    qh_gt = []
    temps = []
    for index, file in enumerate(file_set):
        if file == REACT_PROD_SEP:
            gt[index].append(np.full([len(temps)], np.nan))
            qh_gt[index].append(np.full([len(temps)], np.nan))
        vibes_input = ["python", "-m", "goodvibes", file, "--ti", temp_range]
        if solvent:
            vibes_input += ["-c", "1"]
        if qs_h_s:
            vibes_input += ["-q"]
        vibes_out = subprocess.check_output(vibes_input).decode("utf-8").strip().split("\n")
        found_structure = False
        skip_line = True
        gt.append([])
        qh_gt.append([])
        # we know the last line should be dropped, and at least the first 10
        for line in vibes_out[10:-1]:
            if GOODVIBES_ERROR_PAT.match(line):
                raise InvalidDataError("See GoodVibes error checking report: {}".format(GOODVIBES_OUT_FNAME))
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
                    temps.append(float(vals[2]))
                gt[index].append(float(vals[-2]))
                qh_gt[index].append(float(vals[-1]))
        if save_vibes:
            vibes_out_fname = create_out_fname(file, suffix='_vibes', base_dir=out_dir, ext='.dat')
            os.rename(GOODVIBES_OUT_FNAME, vibes_out_fname)
        if tog_output_fname:
            with open(tog_output_fname, 'a') as f:
                with open(GOODVIBES_OUT_FNAME) as infile:
                    f.write(infile.read())

    temps = np.asarray(temps)
    # for each molecule, multiply the array to convert to kcal/mol
    for index in range(len(gt)):
        gt[index] = np.asarray(gt[index]) * EHPART_TO_KCALMOL
        qh_gt[index] = np.asarray(qh_gt[index]) * EHPART_TO_KCALMOL

    return temps, gt, qh_gt


def get_delta_gibbs(temps, gt, ts_index):
    """
    Calculate the difference in Gibbs free energy at each temperature in temps
    :param temps: np array of temperatures (in K) to be evaluated
    :param gt: list of np arrays with quasi-harmonic treatment of G(T) for each output file at each temp
    :param ts_index: int, index of array with values for ts
    :return: delta_gibbs_ts, delta_gibbs_rxn: np arrays of delta_gibbs
    """
    nan_array = np.full([len(temps)], np.nan)
    gibbs_react = nan_array
    gibbs_ts = nan_array
    gibbs_prod = nan_array
    # any files before the ts will be a reactant file
    for index in range(len(gt)):
        if index == 0:
            gibbs_react = gt[index]
        elif index < ts_index:
            gibbs_react += gt[index]
        elif index == ts_index:
            gibbs_ts = gt[index]
        elif index == ts_index + 1:
            gibbs_prod = gt[index]
        else:
            gibbs_prod += gt[index]
    delta_gibbs_ts = gibbs_ts - gibbs_react
    delta_gibbs_rxn = gibbs_prod - gibbs_react
    return delta_gibbs_ts, delta_gibbs_rxn


def get_kt(temps, delta_gibbs_ts):
    """
    Calculate the rate coefficient at each temperature in temps
    :param temps: np array of temperatures (in K) to be evaluated
    :param delta_gibbs_ts: np array of temperatures (in K) to be evaluated
    :return: kt: np array of rate coefficients at each temp
    """
    # rate coefficient from Eyring equation
    return KB / H * temps * np.exp(-delta_gibbs_ts / RG / temps)  # [1/s]


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
    return a, ea


def print_results(a, ea, qh_a, qh_ea, g_temp, g_ts, g_rxn, qh_g_ts, qh_g_rxn, file_set, out_fname, print_mode):
    file_names = []
    for index in range(-5, 0):
        try:
            file_names.append(os.path.basename(file_set[index]))
        except IndexError:
            file_names.append('')
    output_dict = {FILE1: file_names[0], FILE2: file_names[1], FILE3: file_names[2], FILE4: file_names[3],
                   FILE5: file_names[4],
                   A: a, EA: ea, QH_A: qh_a, QH_EA: qh_ea, DELTA_G_TEMP: g_temp, DELTA_G_TS: g_ts,
                   DELTA_G_RXN: g_rxn, QH_DELTA_G_TS: qh_g_ts, QH_DELTA_G_RXN: qh_g_rxn}
    write_csv([output_dict], out_fname, OUTPUT_HEADERS, extrasaction="ignore", mode=print_mode,
              print_message=False)


def get_delta_g(temp, temps, delta_g_ts, delta_g_rxn):
    """
    Grab delta G's at given temp
    :param temp: float, temp in K to calculate delta G
    :param temps: list of floats, temps for which we have G values
    :param  delta_g_ts: np_array, delta Gibbs free energies between react(s) and TS (kcal/mol) at temps
    :param delta_g_rxn: np_array, delta Gibbs free energies between react(s) and prod(s) (kcal/mol) at temps
    :return: delta_g_ts, delta_g_rxn at requested temp (default to first temp)
    """
    if temp:
        temp_index = (np.abs(temps - temp)).argmin()
    else:
        temp_index = 0

    return temps[temp_index], delta_g_ts[temp_index], delta_g_rxn[temp_index]


def main(argv=None):
    # Read input
    args, ret = parse_cmdline(argv)
    if ret != GOOD_RET or args is None:
        return ret

    try:
        # Make a list of lists; each inner list a set of reactant file(s) with TS
        # Include anything in the "list" file as well as entered on the command line
        if args[0].list:
            with open(args[0].list) as f:
                row_list = [row.strip().split() for row in f.readlines()]
                row_list = list(filter(None, row_list))
        else:
            row_list = []
        if len(args[1]) > 0:
            row_list.append(args[1])
        if len(row_list) == 0:
            raise InvalidDataError("No files or list of files found")

        # now a quick first check that all files exist
        missing_files = []
        for file_set in row_list:
            for file in file_set:
                if not os.path.isfile(file) and file != REACT_PROD_SEP:
                    missing_files.append(file)
        if len(missing_files) > 0:
            raise IOError(missing_files)

        # now the calculations and printing
        print_mode = 'w'  # for the AEa output, so only prints header once, and then appends to file
        if args[0].tog_vibes:
            tog_fname = create_out_fname(args[0].output_fname, suffix='_vibes', ext='.dat')
            # delete if exits because program always appends to it
            silent_remove(tog_fname)
        else:
            tog_fname = None
        for file_set in row_list:
            solvent, ts_index = check_gausslog_fileset(file_set, args[0].hartree_location, args[0].vibes_check)
            temps, gt, qh_gt = get_thermochem(file_set, args[0].temp_range, solvent, args[0].save_vibes,
                                              args[0].out_dir, tog_fname, args[0].quasiharmonic)
            delta_gibbs_ts, delta_gibbs_rxn = get_delta_gibbs(temps, gt, ts_index)
            qh_delta_gibbs_ts, qh_delta_gibbs_rxn = get_delta_gibbs(temps, qh_gt, ts_index)
            if REACT_PROD_SEP in file_set:
                a, ea, qh_a, qh_ea = '', '', '', ''
            else:
                kt = get_kt(temps, delta_gibbs_ts)
                qh_kt = get_kt(temps, qh_delta_gibbs_ts)
                a, ea = fit_arrhenius(temps, kt)
                qh_a, qh_ea = fit_arrhenius(temps, qh_kt)
            g_temp, g_ts, g_rxn = get_delta_g(args[0].temp, temps, delta_gibbs_ts, delta_gibbs_rxn)
            g_temp, qh_g_ts, qh_g_rxn = get_delta_g(args[0].temp, temps, qh_delta_gibbs_ts, qh_delta_gibbs_rxn)
            print_results(a, ea, qh_a, qh_ea, g_temp, g_ts, g_rxn, qh_g_ts, qh_g_rxn,
                          file_set, args[0].output_fname, print_mode)
            print_mode = 'a'

        # clean up GoodVibes detritus
        silent_remove(GOODVIBES_OUT_FNAME)

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
