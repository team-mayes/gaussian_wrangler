#!/usr/bin/env python
"""
if atoms are in the same order, checks for duplicate conformers
"""

import os
import sys
import argparse
from numpy import isnan
from configparser import MissingSectionHeaderError
from common_wrangler.common import (InvalidDataError, warning,
                                    GOOD_RET, INPUT_ERROR, IO_ERROR, INVALID_DATA, DIHES, EHPART_TO_KCAL_MOL, quote,
                                    list_to_file)
from gaussian_wrangler.gw_common import (STOICH, CONVERG, ENERGY, ENTHALPY, GIBBS, CONVERG_ERR, process_gausslog_file)
from gaussian_wrangler import __version__

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
DEF_OUT_NAME = "within_cutoff.txt"


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
    parser.add_argument("-e", "--energy", help="Sort output by lowest electronic energy (not ZPE corrected)."
                                               "The default is False. This flag is superseded by the '-n'/'--enthalpy'"
                                               "and '-g'/'--gibbs' flags.",
                        action='store_true')
    parser.add_argument("-g", "--gibbs", help="Sort output by lowest Gibbs free energy. If not found, the script will "
                                              "sort output by the lowest electronic energy. The default is False.",
                        action='store_true')
    parser.add_argument("-l", "--list", help="The file location of the list of Gaussian output files. There should "
                                             "be one output file listed per line. The default file name is '{}', "
                                             "located in the base directory where the program as "
                                             "run. This program assumes that all the given files have the same atom "
                                             "order.".format(DEF_LIST_FILE),
                        default=DEF_LIST_FILE)
    parser.add_argument("-m", "--max_diff", help="If a numerical value is provided with this option, the output list "
                                                 "will be split between files within or not within this maximum "
                                                 "difference, in kcal/mol), from the lowest energy or enthalpy. "
                                                 "Additionally, the program will output a file with only the file "
                                                 "names of conformations within the cutoff; see the '-o'/'--out_name' "
                                                 "option to specify the name of this file.", default=None)
    parser.add_argument("-n", "--enthalpy", help="Sort output by lowest enthalpy. If no enthalpy is found, it will "
                                                 "sort by the lowest electronic energy. The default is False.",
                        action='store_true')
    parser.add_argument("-o", "--out_fname", help=f"When using the '-m'/'--max_diff' option, a file will be created "
                                                  f"with only the names of the files within the specified cutoff, one "
                                                  f"per line. This option allows the user to specify the output "
                                                  f"file name. By default, the name will be '{DEF_OUT_NAME}'.",
                        default=DEF_OUT_NAME)

    parser.add_argument("-t", "--tol", help="The tolerance, in degrees, for concluding that dihedral angles are "
                                            "equivalent. The default value is {}.".format(DEF_DIH_TOL),
                        default=DEF_DIH_TOL)

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
                if log_info[conf_list[0]][DIHES] is None:
                    add_to_current_group = False
                    continue
                if (log_info[fname][STOICH] != log_info[conf_list[0]][STOICH]) or log_info[fname][DIHES] is None:
                    add_to_current_group = False
                    continue
                add_to_current_group = True
                check_dihes = log_info[conf_list[0]][DIHES]
                for dih_name, dih_val in log_info[fname][DIHES].items():
                    try:
                        dih_diff = abs(dih_val - check_dihes[dih_name])
                        if dih_diff > (360.0 - dih_tol * 1.1):
                            dih_diff -= 360.0
                        if dih_diff > dih_tol:
                            add_to_current_group = False
                            break
                    except KeyError:
                        # may be a different molecule; possible not to reach but safer to leave
                        add_to_current_group = False
                        break
                if add_to_current_group:
                    conf_list.append(fname)
                    break
            if not add_to_current_group:
                conf_groups.append([fname])
    return conf_groups


def print_results(log_info, list_of_conf_lists, sort_by_enthalpy, sort_by_energy, max_diff=None, print_winners=True,
                  out_fname=DEF_OUT_NAME):
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
        winners.append((low_conv_log, log_info[low_conv_log][CONVERG], log_info[low_conv_log][ENERGY],
                        log_info[low_conv_log][ENTHALPY], log_info[low_conv_log][GIBBS]))

    # sorting, if requested
    sort_error = False
    if sort_by_enthalpy:
        sort_by_energy = False
        for winner in winners:
            if isnan(winner[3]):
                sort_by_energy = True
                sort_by_enthalpy = False
                break
    if sort_by_enthalpy:
        sort_key = 3
    elif sort_by_energy:
        sort_key = 2
    else:
        sort_key = 4
    winners.sort(key=lambda tup: tup[sort_key])
    winner_str = quote('","'.join(['File', CONVERG, ENERGY, ENTHALPY, GIBBS]))

    # now gather results
    cutoff_list = []
    if max_diff:
        winner_str += ',"Diff(kcal/mol)"\n'
        lowest_val = winners[0][sort_key]
        if sort_by_enthalpy:
            sort_type = "enthalpy"
        elif sort_by_energy:
            sort_type = "SCF energy"
        else:
            sort_type = "Gibbs free energy"
        winner_str += f'"Files within {sort_type} cutoff of {max_diff:.2f} kcal/mol"\n'
        within_cutoff = True
    else:
        winner_str += '\n'
        lowest_val = None  # to make IDE happy
        within_cutoff = False
    val_diff_str = ""
    val_diff = 0.
    for winner, converg, energy, enthalpy, gibbs, in winners:
        if not sort_error:
            if max_diff:
                if sort_by_enthalpy:
                    val_diff = (enthalpy - lowest_val) * EHPART_TO_KCAL_MOL
                elif sort_by_energy:
                    val_diff = (energy - lowest_val) * EHPART_TO_KCAL_MOL
                else:
                    val_diff = (gibbs - lowest_val) * EHPART_TO_KCAL_MOL
                val_diff_str = f",{val_diff:.2f}"

            if within_cutoff:
                if val_diff > max_diff:
                    winner_str += f'"Files outside of cutoff:"\n'
                    within_cutoff = False
                else:
                    cutoff_list.append(winner)

            winner_str += f'"{winner}",{converg:.4f},{energy:.6f},{enthalpy:.6f},{gibbs:.6f}{val_diff_str}\n'
        if log_info[winner][CONVERG_ERR]:
            warn_files_str += '\n    {:}:  {:.2f}'.format(winner, converg)
        elif log_info[winner][CONVERG_ERR] is None:
            warn_files_str += '\n    {:}:  Not found'.format(winner)
    if print_winners:
        print(winner_str)

    if cutoff_list:
        list_to_file(cutoff_list, out_fname)
    return winner_str, warn_files_str


def main(argv=None):
    print(f"Running GaussianWrangler script gausslog_unique version {__version__}")
    # Read input
    args, ret = parse_cmdline(argv)
    if ret != GOOD_RET or args is None:
        return ret

    # Read template and data files
    try:
        gausslog_files = []
        missing_files = []
        log_info = {}

        # check input
        if args.max_diff:
            args.max_diff = float(args.max_diff)
            if not args.energy and not args.gibbs:
                args.enthalpy = True

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
        winner_str, warn_files_str = print_results(log_info, list_of_conf_lists, args.enthalpy, args.energy,
                                                   args.max_diff, args.out_fname)
        if len(warn_files_str) > 0:
            warning("Check convergence of file(s):" + warn_files_str)

    except IOError as e:
        warning("Problems reading file:", e)
        return IO_ERROR
    except (InvalidDataError, UnicodeDecodeError) as e:
        warning("Problems reading data:", e)
        return INVALID_DATA
    except ValueError as e:
        warning(e.args[0])
        return INVALID_DATA
    return GOOD_RET  # success


if __name__ == '__main__':
    status = main()
    sys.exit(status)
