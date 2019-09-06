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
GAUSSIAN_SEPARATOR = '******************************************'
GOODVIBES_OUT_FNAME = "Goodvibes_output.dat"
GOODVIBES_ERROR_PAT = re.compile(r"x .*")
GOODVIBES_DATA_PAT = re.compile(r" {3}Structure .*")

# for printing
FILE1 = 'file1'
FILE2 = 'file2'
FILE3 = 'file3'
FILE4 = 'file4'
A = 'A (1/s if uni)'
EA = 'Ea (kcal/mol)'
OUTPUT_HEADERS = [FILE1, FILE2, FILE3, FILE4, A, EA, ]


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
                                                "is the current working directory.", default=False)
    parser.add_argument("-hl", "--hartree_location", help="Optional: the call to invoke hartree. The default is: "
                                                          "{}".format(DEF_HARTREE_LOC), default=DEF_HARTREE_LOC)
    parser.add_argument("-l", "--list", help="The location of the list of Gaussian output files. "
                                             "The default file name.", default=None)
    parser.add_argument("-ti", "--temp_range", help="Initial temp, final temp, (and optionally) step size (K) for "
                                                    "thermochemistry calculations. The default range is 300,600,30",
                        default="300,600,30")
    parser.add_argument("-o", "--output_fname", help="The name of the output file to be created. The default is {}"
                                                     "".format(DEF_OUT_FILE_NAME), default=DEF_OUT_FILE_NAME)
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
    # Runs check that there is one terminal TS file, then calls GoodVibes for its checks
    if len(file_set) < 2 or len(file_set) > 4:
        raise InvalidDataError("Expected 2,3, or 4 files in a set, but found {}: {}".format(len(file_set), file_set))

    total_react_charge = 0
    ts_charge = np.nan
    multiplicities = np.full([len(file_set)], np.nan)

    # initialize the variables below to make my IDE happy
    solvent = None
    func = None
    basis = None
    gauss_ver = None
    react_stoich_dict = {}  # empty dict to make IDE happy
    ts_stoich_dict = {}  # empty dict to make IDE happy
    for index, fname in enumerate(file_set):
        hartree_output = subprocess.check_output(["java", "-jar", hartree_loc,
                                                  "snap", "-f", fname]).decode("utf-8").strip().split("\n")

        hartree_list = list(csv.DictReader(hartree_output, quoting=csv.QUOTE_NONNUMERIC))
        for row in hartree_list:
            multiplicities[index] = int(row[MULT])
            if index < len(file_set) - 1:
                if float(row[FREQ1]) < 0 or float(row[FREQ2]) < 0:
                    raise InvalidDataError("Expected no imaginary frequencies in file: {}".format(fname))
                total_react_charge += int(row[CHARGE])
                if index == 0:
                    react_stoich_dict = parse_stoich(row[STOICH])
                else:
                    react_stoich_dict = parse_stoich(row[STOICH], add_to_dict=react_stoich_dict)
            else:
                if float(row[FREQ1]) > 0 or float(row[FREQ2]) < 0:
                    raise InvalidDataError("Expected one imaginary frequency in file: {}".format(fname))
                ts_charge = int(row[CHARGE])
                ts_stoich_dict = parse_stoich(row[STOICH])

            # additional checks as we go...
            awk_command = ['awk', '/\*\*\*/{getline; print; exit}', fname]
            if index == 0:
                solvent = row[SOLV]
                func = row[FUNCTIONAL]
                basis = row[BASIS_SET]
                gauss_ver = subprocess.check_output(awk_command).strip().split()[:3]
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
                this_gauss_ver = subprocess.check_output(awk_command).strip().split()[:3]
                if gauss_ver != this_gauss_ver:
                    raise InvalidDataError("Different Gaussian versions ({}, {}) found for file set: "
                                           "{}".format(gauss_ver, this_gauss_ver, file_set))

    if react_stoich_dict != ts_stoich_dict:
        raise InvalidDataError("Check stoichiometries of reactant(s) and transition state for set: {}".format(file_set))

    if good_vibes_check:
        if solvent:
            file_set = file_set + ["-c", "1"]
        vibes_out = subprocess.check_output(["python", "-m", "goodvibes"] + file_set +
                                            ["--check"]).decode("utf-8").strip().split("\n")
        for line in vibes_out:
            if GOODVIBES_ERROR_PAT.match(line):
                if 'Different charge and multiplicity' in line:
                    # check if multiplicities all the same
                    mult_check = np.sum(multiplicities - multiplicities[0])
                    if mult_check == 0 and total_react_charge == ts_charge:
                        continue
                raise InvalidDataError("See GoodVibes error checking report: 'Goodvibes_output.dat'")

    return solvent


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


def get_thermochem(file_set, temp_range, solvent, save_vibes, out_dir, tog_output_fname):
    """
    Calls GoodVibes to get thermochem at a range of temps
    :param file_set: list of reactant file(s) and TS file
    :param temp_range: string with range of temperatures at which to calculate thermochem
    :param solvent: boolean to decide whether to include
    :param save_vibes: boolean to determine whether to save each GoodVibes output separately
    :param out_dir: directory to save GoodVibes output files (if requested)
    :param tog_output_fname: None or string (file name) if saving each GoodVibes output together
    :return: nothing
    """
    qh_gt = []
    temps = []
    for index, file in enumerate(file_set):
        vibes_input = ["python", "-m", "goodvibes", file, "--ti", temp_range]
        if solvent:
            # Todo: see why this doesn't change the answers
            #  Check whether I need to add the factor Rkcal*temp*LN(Ratm*temp)
            vibes_input += ["-c", "1"]
        vibes_out = subprocess.check_output(vibes_input).decode("utf-8").strip().split("\n")
        found_structure = False
        skip_line = True
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
                qh_gt[index].append(float(vals[-1]))
        if save_vibes:
            vibes_out_fname = create_out_fname(file, suffix='_vibes', base_dir=out_dir, ext='.dat')
            os.rename(GOODVIBES_OUT_FNAME, vibes_out_fname)
        if tog_output_fname:
            with open(tog_output_fname, 'a') as f:
                with open(GOODVIBES_OUT_FNAME) as infile:
                    f.write(infile.read())

    temps = np.asarray(temps)
    for index in range(len(qh_gt)):
        qh_gt[index] = np.asarray(qh_gt[index]) * EHPART_TO_KCALMOL

    return temps, qh_gt


def get_kt(temps, qh_gt):
    """
    Calculate the rate coefficient at each temperature in temps
    :param temps: np array of temperatures (in K) to be evaluated
    :param qh_gt: list of np arrays with quasi-harmonic treatment of G(T) for each output file at each temp
    :return: kt: np array of rate coefficients at each temp
    """
    gibbs_react = qh_gt[0]
    # any files before the last will be a reactant file
    for index in range(len(qh_gt) - 1):
        if index == 0:
            gibbs_react = qh_gt[0]
        else:
            gibbs_react += qh_gt[index]
    gibbs_ts = qh_gt[-1]
    delta_gibbs = gibbs_ts - gibbs_react

    # rate coefficient from Eyring equation
    return KB / H * temps * np.exp(-delta_gibbs / RG / temps)  # [1/s]


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


def print_results(a, ea, file_set, out_fname, print_mode):
    file_names = []
    for index in range(-4, 0):
        try:
            file_names.append(os.path.basename(file_set[index]))
        except IndexError:
            file_names.append('')
    output_dict = {FILE1: file_names[0], FILE2: file_names[1], FILE3: file_names[2], FILE4: file_names[3],
                   A: a, EA: ea, }
    write_csv([output_dict], out_fname, OUTPUT_HEADERS, extrasaction="ignore", mode=print_mode,
              print_message=False)


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
                if not os.path.isfile(file):
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
            solvent = check_gausslog_fileset(file_set, args[0].hartree_location, args[0].vibes_check)
            temps, qh_gt = get_thermochem(file_set, args[0].temp_range, solvent, args[0].save_vibes, args[0].out_dir,
                                          tog_fname)
            kt = get_kt(temps, qh_gt)
            a, ea = fit_arrhenius(temps, kt)
            print_results(a, ea, file_set, args[0].output_fname, print_mode)
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
