#!/usr/bin/env python
"""
Creates pdb data files from lammps data files, given a template pdb file.
"""

from __future__ import print_function
import os
import subprocess
import sys
import argparse
import tempfile
from nrel_tools.common import (InvalidDataError, warning, create_out_fname, list_to_file, ATOM_NUM_DICT,
                               NUM_ATOMS, GAU_COORD_PAT, GAU_SEP_PAT, GAU_E_PAT, read_csv,
                               GOOD_RET, INPUT_ERROR, IO_ERROR, INVALID_DATA, GAU_CHARGE_PAT)

try:
    # noinspection PyCompatibility
    from ConfigParser import ConfigParser, MissingSectionHeaderError
except ImportError:
    # noinspection PyCompatibility
    from configparser import ConfigParser, MissingSectionHeaderError

__author__ = 'hmayes'


# Constants #


# Config keys
DEF_LIST_FILE = 'log_list.txt'
GAUSSLOG_FILES_FILE = 'gausslog_list_file'
OUT_BASE_DIR = 'output_directory'
DEF_HARTREE_LOC = '/Users/hmayes/.local/bin/hartree-cli-1.2.4.jar'

# For log file processing
SEC_HEAD = 'head_section'
SEC_ATOMS = 'atoms_section'
SEC_TAIL = 'tail_section'
COM_TYPE = 'com_type'
CHARGE = 'charge'
MULTIPLICITY = 'multiplicity'
FREQ1 = 'Freq 1'
FREQ2 = 'Freq 2'


def str2bool(v):
    if isinstance(v, bool):
        return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')


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
    parser.add_argument("-hl", "--hartree_location", help="Optional: the call to invoke hartree. The default is: "
                                                          "{}".format(DEF_HARTREE_LOC), default=DEF_HARTREE_LOC)
    parser.add_argument("-l", "--list", help="The location of the list of Gaussian output files. "
                                             "The default file name.", default=False)
    args = None
    try:
        args = parser.parse_known_args(argv)
    except SystemExit as e:
        if hasattr(e, 'code') and e.code == 0:
            return args, GOOD_RET
        warning(e)
        parser.print_help()
        return args, INPUT_ERROR

    return args, GOOD_RET


def process_gausslog_file_set(file_set, hartree_loc):
    if len(file_set) == 2:
        react_type = 1  # for uni-molecular
    elif len(file_set) == 3:
        react_type = 2  # for bi-molecular
    elif len(file_set) == 4:
        react_type = 3  # for tri-molecular
    else:
        raise InvalidDataError("Expected 2,3, or 4 files in a set, but found {}: {}".format(len(file_set), file_set))

    hartree_csv = tempfile.NamedTemporaryFile()

    hartree_out_list = []
    for fname in file_set:
        hartree_output = subprocess.check_output(["java", "-jar", hartree_loc,
                                                  "snap", "-f",  fname]).decode("utf-8").strip().split("\n")
        list_to_file(hartree_output, hartree_csv, print_message=False)
        hartree_dict = read_csv(hartree_csv)

        # keys = [x.strip('"')for x in hartree_output[0].split(",")]
        # vals = [x.strip('"')for x in hartree_output[1].split(",")]
        # hartree_dict = dict(zip(keys, vals))
        # hartree_out_list.append(hartree_dict)

    for index in range(0, react_type):
        print(hartree_out_list[index][FREQ1], hartree_out_list[index][FREQ2])

    return react_type


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
        
        for file_set in row_list:
            process_gausslog_file_set(file_set, args[0].hartree_location)

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
