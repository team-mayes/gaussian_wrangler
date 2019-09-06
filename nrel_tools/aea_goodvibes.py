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
                               write_csv, silent_remove)

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
STOICH = 'Stoichiometry'
SOLV = 'Solvent type'
CHARGE = 'Charge'
MULT = 'Mult'
FREQ1 = 'Freq 1'
FREQ2 = 'Freq 2'
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
                                             "The default file name.", default=False)
    parser.add_argument("-ti", "--temp_range", help="Initial temp, final temp, (and optionally) step size (K) for "
                                                    "thermochemistry calculations. The default range is 300,600,30",
                        default="300,600,30")
    parser.add_argument("-o", "--output_fname", help="The name of the output file to be created. The default is {}"
                                                     "".format(DEF_OUT_FILE_NAME), default=DEF_OUT_FILE_NAME)

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


def check_gausslog_fileset(file_set, hartree_loc):
    """
    checks include:
       using hartree to get info to check for:
           the correct number of imaginary freq
           the stoichiometry adds up   ##TODO
        using GoodVibes for its various checks (all same solvent, level of theory, etc.)
    :param file_set: list of reactant file(s) and TS file
    :param hartree_loc: location where hartree can be found to run
    :return: reaction_type: integer for molecularity of reaction
    """
    # Runs check that there is one terminal TS file, then calls GoodVibes for its checks
    if len(file_set) < 2 or len(file_set) > 4:
        raise InvalidDataError("Expected 2,3, or 4 files in a set, but found {}: {}".format(len(file_set), file_set))

    # Made the empty list below to make my IDE happy
    hartree_list = []
    total_react_charge = 0
    ts_charge = np.nan
    multiplicities = np.full([len(file_set)], np.nan)
    solvent = False
    for index, fname in enumerate(file_set):
        hartree_output = subprocess.check_output(["java", "-jar", hartree_loc,
                                                  "snap", "-f",  fname]).decode("utf-8").strip().split("\n")

        hartree_list = list(csv.DictReader(hartree_output, quoting=csv.QUOTE_NONNUMERIC))
        for row in hartree_list:
            multiplicities[index] = int(row[MULT])
            if index < len(file_set)-1:
                if float(row[FREQ1]) < 0 or float(row[FREQ2]) < 0:
                    raise InvalidDataError("Expected no imaginary frequencies in file: {}".format(fname))
                total_react_charge += int(row[CHARGE])
            else:
                if float(row[FREQ1]) > 0 or float(row[FREQ2]) < 0:
                    raise InvalidDataError("Expected one imaginary frequency in file: {}".format(fname))
                ts_charge = int(row[CHARGE])

    if hartree_list[0][SOLV] != 'N/A':
        file_set = file_set + ["-c", "1"]
        solvent = True

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


def get_thermochem(file_set, temp_range, solvent):
    """
    Calls GoodVibes to get thermochem at a range of temps
    :param file_set: list of reactant file(s) and TS file
    :param temp_range: string with range of temperatures at which to calculate thermochem
    :param solvent: boolean to decide whether to include
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
    for index in range(len(qh_gt)-1):
        if index == 0:
            gibbs_react = qh_gt[0]
        else:
            gibbs_react += qh_gt[index]
    gibbs_ts = qh_gt[-1]
    delta_gibbs = gibbs_ts - gibbs_react

    # rate coefficient from Eyring equation
    return KB / H * temps * np.exp(-delta_gibbs/RG/temps)  # [1/s]


def fit_arrhenius(temps, kt):
    """
    Fits a straight line of 1/temps vs. ln(kt) to get A (1/s), Ea (kcal/mol)
    :param temps: numpy array of temps in K
    :param kt: numpy array of rate coefficients in 1/s (or appropriate units)
    :return: calcs A (1/s), Ea (kcal/mol)
    """
    inv_temp = 1/temps
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

        # now the calculations and printing
        print_mode = 'w'
        for file_set in row_list:
            solvent = check_gausslog_fileset(file_set, args[0].hartree_location)
            temps, qh_gt = get_thermochem(file_set, args[0].temp_range, solvent)
            kt = get_kt(temps, qh_gt)
            a, ea = fit_arrhenius(temps, kt)
            print_results(a, ea, file_set, args[0].output_fname, print_mode)
            print_mode = 'a'

        # clean up GoodVibes detritus
        silent_remove("Goodvibes_output.dat")

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
