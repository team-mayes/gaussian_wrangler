#!/usr/bin/env python
"""
if atoms are in the same order, checks for duplicate conformers
"""

from __future__ import print_function
import os
import sys
import argparse
from common_wrangler.common import (InvalidDataError, warning,
                                    GOOD_RET, INPUT_ERROR, IO_ERROR, INVALID_DATA, DIHES)
from gaussian_wrangler.gw_common import (STOICH, CONVERG, ENERGY, ENTHALPY, CONVERG_ERR, process_gausslog_file)

try:
    # noinspection PyCompatibility
    from ConfigParser import ConfigParser, MissingSectionHeaderError
except ImportError:
    # noinspection PyCompatibility
    from configparser import ConfigParser, MissingSectionHeaderError

__author__ = 'hmayes'


# Constants #

# Config File Sections

# Config keys
LOG_LIST_FILE = 'input_log_list_file'
DIH_TOL = 'dihedral_tolerance'


# data file info

# Defaults
DEF_LIST_FILE = 'list.txt'
DEF_DIH_TOL = 5.0


# For file processing
MAX_BOND_DIST = 1.9  # same length units as in input and output file, here Angstroms
MAX_H_BOND_DIST = 1.5  # same length units as in input and output file, here Angstroms
MAX_M_BOND_DIST = 2.3  # same length units as in input and output file, here Angstroms
METALS = ['Ti', 'Sb', 'Ge']


def parse_cmdline(argv):
    """
    Returns the parsed argument list and return code.
    `argv` is a list of arguments, or `None` for ``sys.argv[1:]``.
    """
    if argv is None:
        argv = sys.argv[1:]

    # initialize the parser object:
    parser = argparse.ArgumentParser(description='Given a list of Gaussian output files, returns a list with '
                                                 'unique conformers, defined by having dihedral angles within the '
                                                 'specified tolerance.')
    parser.add_argument("-l", "--list", help="The file location of the list of Gaussian output files. There should "
                                             "be one output file listed per line. The default file name is '{}', "
                                             "located in the base directory where the program as "
                                             "run. This program assumes that all the given files have the same atom "
                                             "order.".format(DEF_LIST_FILE),
                        default=DEF_LIST_FILE)
    parser.add_argument("-t", "--tol", help="The tolerance, in degrees, for concluding that dihedral angles are "
                                            "equivalent. The default value is {}.".format(DEF_DIH_TOL),
                        default=DEF_DIH_TOL)
    parser.add_argument("-e", "--energy", help="Sort output by lowest electronic energy (not ZPE corrected)."
                                               "The default is False. This flag is superseded by the enthalpy flag.",
                        action='store_true')
    parser.add_argument("-n", "--enthalpy", help="Sort output by lowest enthalpy. If no enthalpy is found, it will "
                                                 "sort by the lowest electronic energy. The default is False.",
                        action='store_true')

    args = None
    try:
        args = parser.parse_args(argv)
    except (KeyError, InvalidDataError, MissingSectionHeaderError, SystemExit) as e:
        if hasattr(e, 'code') and e.code == 0:
            return args, GOOD_RET
        warning(e)
        parser.print_help()
        return args, INPUT_ERROR

    return args, GOOD_RET


def compare_gausslog_info(log_info, dih_tol):
    f_names = log_info.keys()
    conf_groups = []
    for f_index, fname in enumerate(f_names):
        if f_index == 0:
            conf_groups.append([fname])
        else:
            add_to_current_group = True  # not really necessary, but IDE won't complain this way
            for conf_list in conf_groups:
                if log_info[fname][STOICH] != log_info[conf_list[0]][STOICH]:
                    add_to_current_group = False
                    continue
                add_to_current_group = True
                check_dihes = log_info[conf_list[0]][DIHES]
                for dih_name, dih_val in log_info[fname][DIHES].items():
                    try:
                        dih_diff = abs(dih_val - check_dihes[dih_name])
                        if dih_diff > (360.0 - dih_tol):
                            dih_diff -= 360.0
                        if dih_diff > dih_tol:
                            add_to_current_group = False
                            break
                    except KeyError:
                        # probably an isomer
                        add_to_current_group = False
                        break
                if add_to_current_group:
                    conf_list.append(fname)
                    break
            if not add_to_current_group:
                conf_groups.append([fname])
    return conf_groups


def print_results(log_info, list_of_conf_lists, sort_by_enthalpy, sort_by_energy):
    winners = []
    warn_files_str = ''
    for conf_list in list_of_conf_lists:
        if len(conf_list) == 1:
            low_conv_log = conf_list[0]
        else:
            lowest_converg = 20000000.0
            low_conv_log = None  # here to make IDE happy
            for log_file in conf_list:
                if log_info[log_file][CONVERG] < lowest_converg:
                    lowest_converg = log_info[log_file][CONVERG]
                    low_conv_log = log_file
        winners.append((low_conv_log, log_info[low_conv_log][CONVERG],
                        log_info[low_conv_log][ENERGY], log_info[low_conv_log][ENTHALPY]))

    # sorting, if requested
    if sort_by_enthalpy:
        sort_by_energy = False
        for winner in winners:
            if winner[3] is None:
                sort_by_energy = True
                sort_by_enthalpy = False
                break
    if sort_by_enthalpy:
        sort_key = 3
    elif sort_by_energy:
        sort_key = 2
    else:
        sort_key = 0
    winners.sort(key=lambda tup: tup[sort_key])

    # now print results
    print(','.join(['File', CONVERG, ENERGY, ENTHALPY]))
    for winner, converg, energy, enthalpy in winners:
        try:
            print('{},{:.4f},{:.6f},{:.6f}'.format(winner, converg, energy, enthalpy))
        except TypeError:
            print('{},{:.4f},{:.6f},{}'.format(winner, converg, energy, enthalpy))
        if log_info[winner][CONVERG_ERR]:
            warn_files_str += '\n    {:}:  {:.4f}'.format(winner, converg, energy, enthalpy)
    if len(warn_files_str) > 0:
        warning("Check convergence of file(s):" + warn_files_str)


def main(argv=None):
    # Read input
    args, ret = parse_cmdline(argv)
    if ret != GOOD_RET or args is None:
        return ret

    # Read template and data files
    try:
        gausslog_files = []
        missing_files = []
        log_info = {}

        # check that we have files
        with open(args.list) as f:
            for line in f:
                fname = line.strip()
                if len(fname) == 0:
                    continue
                # check that each log file can be found
                if os.path.isfile(fname):
                    gausslog_files.append(fname)
                else:
                    missing_files.append(fname)
            if len(missing_files) > 0:
                raise IOError("Could not find the following file(s) listed in '{}':\n    "
                              "{}".format(args.list, '\n    '.join(sorted(set(missing_files)))))
            if len(gausslog_files) < 2:
                raise InvalidDataError("This program expects at least two files to compare to determine if they "
                                       "have the same conformation. Check input.")

        # get the data from the files
        for gausslog_file in gausslog_files:
            gausslog_content = process_gausslog_file(gausslog_file, find_dih=True, find_converg=True)
            log_info[os.path.basename(gausslog_file)] = gausslog_content

        # process data from files
        list_of_conf_lists = compare_gausslog_info(log_info, args.tol)
        print_results(log_info, list_of_conf_lists, args.enthalpy, args.energy)

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
