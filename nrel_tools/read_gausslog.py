#!/usr/bin/env python
"""
Creates pdb files from Gaussian log files
"""

from __future__ import print_function
import sys
import argparse
from nrel_tools.common import (InvalidDataError, warning, process_cfg, create_out_fname, list_to_file, process_pdb_file,
                               SEC_HEAD, SEC_ATOMS, SEC_TAIL, ATOM_NUM_DICT,
                               GAU_E_PAT,
                               GOOD_RET, INPUT_ERROR, IO_ERROR, INVALID_DATA, silent_remove)

try:
    # noinspection PyCompatibility
    from ConfigParser import ConfigParser, MissingSectionHeaderError
except ImportError:
    # noinspection PyCompatibility
    from configparser import ConfigParser, MissingSectionHeaderError

__author__ = 'hmayes'

# Constants #

# Config File Sections
MAIN_SEC = 'main'

# Config keys
GAUSSLOG_FILES = 'gausslog_files_list'
CP_JOB = 'cp'
STABILITY_JOB = 'stable'
JOB_TYPES = [CP_JOB, STABILITY_JOB]

# data file info
FILE_NAME = "File Name"
SOLVENT = "Solvent type"
STOICH = "Stoichiometry"
CHARGE = "Charge"
MULT = "Mult"
FUNC = "Functional"
BASIS = "Basis Set"
ENERGY = "Energy (Hartrees)"
DIPOLE = "Dipole (Debye)"
ZPE = "ZPE (Hartrees)"
CP = "Counterpoise Correction (Hartrees)"
H298 = "H298 (Hartrees)"
G298 = "G298 (Hartrees)"
FREQ1 = "Freq 1"
FREQ2 = "Freq 2"


def parse_cmdline(argv):
    """
    Returns the parsed argument list and return code.
    `argv` is a list of arguments, or `None` for ``sys.argv[1:]``.
    """
    if argv is None:
        argv = sys.argv[1:]

    # initialize the parser object:
    parser = argparse.ArgumentParser(description='Creates pdb files from Gaussian input files, given a template pdb '
                                                 'file.')
    parser.add_argument("-f", "--file", help="The location of the Gaussian log file to be read.",
                        default=None)
    parser.add_argument("-j", "--job_type", help="The type of calculation to be read. Options are {}".format(JOB_TYPES),
                        default=None)
    parser.add_argument("-l", "--list", help="The location of the file with the list of Gaussian log files to be read.",
                        default=None)
    args = None
    try:
        args = parser.parse_args(argv)
    except IOError as e:
        warning("Problems reading file:", e)
        parser.print_help()
        return args, IO_ERROR
    except (KeyError, InvalidDataError, MissingSectionHeaderError, SystemExit) as e:
        if hasattr(e, 'code') and e.code == 0:
            return args, GOOD_RET
        warning(e)
        parser.print_help()
        return args, INPUT_ERROR

    return args, GOOD_RET


def process_gausscom_files(gausslog_file, gausslog_file_list, job_type):
    log_file_list = []
    if gausslog_file:
        log_file_list.append(gausslog_file)
    if gausslog_file_list:
        with open(gausslog_file_list) as f:
            for data_file in f:
                log_file_list.append(data_file.strip())
    # for file in log_file_list:
    #     process_gausslog_file(file, job_type)


# def process_gausslog_file(cfg, gausslog_file, pdb_tpl_content, f_name):
#     with open(gausslog_file) as d:
#         section = SEC_HEAD
#         atom_id = 0
#         lines_after_coord = 2  # blank line, description, blank line, charge & multiplicity
#         message = True
#         coord_match = False
#         first_pass = True
#
#         for line in d:
#             line = line.strip()
#             if len(line) == 0:
#                 continue
#             if GAU_E_PAT.match(line):
#                 pdb_tpl_content[SEC_HEAD].append("REMARK    {}".format(line))
#                 if not cfg[ONLY_FINAL]:
#                     check_and_print(cfg, atom_id, pdb_tpl_content, gausslog_file, pdb_data_section,
#                                     f_name, mode, message)
#                     message = False
#                     mode = 'a'
#                 section = SEC_HEAD
#                 coord_match = False
#                 atom_id = 0
#                 lines_after_coord = 2
#
#     if cfg[ONLY_FINAL]:
#         check_and_print(cfg, atom_id, pdb_tpl_content, gausslog_file, pdb_data_section,
#                         f_name, mode, message)


def main(argv=None):
    # Read input
    args, ret = parse_cmdline(argv)
    if ret != GOOD_RET or args is None:
        return ret

    # Read template and data files
    try:
        process_gausscom_files(args.file, args.list, args.job_type)
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
