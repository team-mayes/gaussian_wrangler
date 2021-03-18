#!/usr/bin/env python
"""
Finds Gaussian files in a folder. If complete (Normal termination), moves them to a subdirectory.
If appears to have failed, add to a text file for investigation.
"""

import os
import re
import sys
import argparse
import matplotlib.pyplot as plt
import numpy as np
import warnings
from scipy.optimize import curve_fit, OptimizeWarning
from scipy.interpolate import make_interp_spline
from scipy.signal import argrelmin, argrelmax
from operator import itemgetter
from configparser import MissingSectionHeaderError
from common_wrangler.common import (GOOD_RET, INPUT_ERROR, IO_ERROR, INVALID_DATA, EHPART_TO_KCAL_MOL, COLOR_SEQUENCE,
                                    InvalidDataError, warning, create_out_fname, write_csv, check_for_files,
                                    assign_color)

from gaussian_wrangler.gw_common import (MAX_FORCE, RMS_FORCE, MAX_DISPL, RMS_DISPL, CONVERG, CONVERG_ERR,
                                         process_gausslog_file, CONVERG_STEP_DICT, ENERGY, SCAN_DICT)
from gaussian_wrangler import __version__


__author__ = 'hmayes'


# Constants #
NORM_TERM_PAT = re.compile(r"Normal termination of Gaussian*")
FAIL_OPEN_FILE = re.compile(r"open-new-file*")
FAIL_RDCARD = re.compile(r"In source file rdcard*")
FAIL_NTR = re.compile(r"NtrErr Called from FileIO*")
FAIL_LEN_PAT = re.compile(r"File lengths (MBytes): {2}RWF=.*")
FAIL_LEN = re.compile(r"File lengths (MBytes)*")
FAIL_PAT_LIST = [FAIL_OPEN_FILE, FAIL_LEN_PAT, FAIL_RDCARD, FAIL_NTR, FAIL_LEN]

DEF_COMPLETE_DIR = 'for_hartree'
DEF_EXT = '.log'
OUT_BASE_DIR = 'output_directory'

# For convergence check
STEP_NUM = 'step_number'
F_NAME = 'File Name'
STEP_CONVERG_HEADERS = [F_NAME, STEP_NUM, MAX_FORCE, RMS_FORCE, MAX_DISPL, RMS_DISPL, CONVERG, CONVERG_ERR]
FINAL_CONVERG_HEADERS = [F_NAME, CONVERG, CONVERG_ERR]

N_DIHE = np.asarray([1, 2, 3, 4, 6])


def parse_cmdline(argv):
    """
    Returns the parsed argument list and return code.
    `argv` is a list of arguments, or `None` for ``sys.argv[1:]``.
    """
    if argv is None:
        argv = sys.argv[1:]

    # initialize the parser object:
    parser = argparse.ArgumentParser(description='Checks for normal termination of Gaussian output files in a '
                                                 'specified directory, and moves them to a new location.')
    parser.add_argument("-a", "--all", help="Check convergence of all steps and print to standard out.",
                        action="store_true", default=False)
    parser.add_argument("-b", "--best", help="Check convergence of each step and list the convergence of the best 10 "
                                             "steps, sorted by convergence.", action="store_true", default=False)
    parser.add_argument("-d", "--directory", help="The directory where to look for Gaussian output files to check for "
                                                  "normal termination, without checking in subdirectories.",
                        metavar="path", default=None)
    parser.add_argument("-ds", "--dir_subdirs", help="The directory where to look for Gaussian output files to check "
                                                     "for normal termination, including checking in subdirectories.",
                        metavar="path", default=None)
    parser.add_argument("-e", "--extension", help="The extension of the Gaussian output file(s) to look for when "
                                                  "searching a directory for output files. The default is '{}'."
                                                  "".format(DEF_EXT), metavar="ext", default=DEF_EXT)
    parser.add_argument("-f", "--file_name", help="A file name (with path, if not the current directory) to check for "
                                                  "either normal termination or convergence. If used, this option "
                                                  "overrides the '-d' option, and no searching for files is "
                                                  "performed.", metavar="path", default=None)
    parser.add_argument("-l", "--file_list", help="A file name (with path, if not the current directory) with a "
                                                  "list of files (also with path, if not the current directory)  "
                                                  "overrides the '-d' option, and no searching for files is to check "
                                                  "for either normal termination or convergence. If used, this "
                                                  "option overrides the '-d' option, and no searching for files is "
                                                  "performed.", metavar="path", default=None)
    parser.add_argument("-o", "--output_directory", help="The directory where to put Gaussian output files that have "
                                                         "terminated normally. The default is '{}'."
                                                         "".format(DEF_COMPLETE_DIR), metavar="path",
                        default=DEF_COMPLETE_DIR)
    parser.add_argument("-s", "--step_converg", help="Report the convergence for each step value for the files in the "
                                                     "directory or those specified with the '-f' or '-l' options. When "
                                                     "this option is chosen, the check for normal termination is "
                                                     "skipped. The default is False.",
                        action="store_true", default=False)
    parser.add_argument("-t", "--to_step", help="Check convergence of each step only to provided step number, and "
                                                "before printing to standard out, sort by convergence.",
                        default=False)
    parser.add_argument("-z", "--final_converg", help="Report the final convergence value for the files in the "
                                                      "directory or those specified with the '-f' or '-l' options. "
                                                      "When this option is chosen, the check for normal termination "
                                                      "is skipped. The default is False.", action="store_true",
                        default=False)
    parser.add_argument("--scan", help="Read output file(s) from a scan and writes the converged energies from each "
                                       "point of the scan to a csv file and creates a plot saved as the given file "
                                       "name.", metavar="path", default=None)
    args = None
    try:
        args = parser.parse_args(argv)
        if args.to_step or args.best or args.all:
            args.step_converg = True
        if args.to_step:
            try:
                args.to_step = int(args.to_step)
            except ValueError:
                raise InvalidDataError("When the '-t' option is used, an integer must be provided.")
        if args.step_converg and args.final_converg:
            raise InvalidDataError("Choose either the '-a', '-b', '-s', '-t', or '-z' option.")
        # make the default output directory a subdirectory of the directory to search
        if args.output_directory == DEF_COMPLETE_DIR:
            if args.dir_subdirs:
                args.output_directory = os.path.relpath(os.path.join(args.dir_subdirs, DEF_COMPLETE_DIR))
            if args.directory:
                args.output_directory = os.path.relpath(os.path.join(args.directory, DEF_COMPLETE_DIR))

    except (KeyError, InvalidDataError, MissingSectionHeaderError, SystemExit) as e:
        if hasattr(e, 'code') and e.code == 0:
            return args, GOOD_RET
        warning(e)
        parser.print_help()
        return args, INPUT_ERROR

    return args, GOOD_RET


def check_file_termination(output_file, good_output_dir, completed_list, likely_failed_list, perhaps_running_list):
    try:
        with open(output_file, 'r') as fh:
            last_line = fh.readlines()[-1].strip()
    except IndexError:
        warning("Could not read the last line (may be blank) of file: {}".format(output_file))
        return
    if NORM_TERM_PAT.match(last_line):
        base_name = os.path.basename(output_file)
        completed_list.append(output_file)
        os.rename(output_file, os.path.join(good_output_dir, base_name))
        return
    for pattern in FAIL_PAT_LIST:
        if pattern.match(last_line):
            likely_failed_list.append(output_file)
            return
    perhaps_running_list.append(output_file)


def check_termination(args, check_file_list):
    completed_list = []
    perhaps_running_list = []
    likely_failed_list = []

    for fname in check_file_list:
        check_file_termination(fname, args.output_directory, completed_list, likely_failed_list, perhaps_running_list)
    # sort if list is at least 2 long:
    for file_list in [completed_list, likely_failed_list, perhaps_running_list]:
        if len(file_list) > 1:
            file_list.sort()
    if len(completed_list) > 0:
        print("The following files completed normally:")
        for fname in completed_list:
            print("    {}".format(os.path.relpath(fname)))
    else:
        print("No normally completed files found.")
    if len(likely_failed_list) > 0:
        print("The following files may have failed:")
        for fname in likely_failed_list:
            print("    {}".format(os.path.relpath(fname)))
    if len(perhaps_running_list) > 0:
        print("The following files may still be running:")
        for fname in perhaps_running_list:
            print("    {}".format(os.path.relpath(fname)))


def create_convergence_plots(out_fname, step_list):
    """
    To allow easy viewing of convergence
    :param out_fname: This is the name of the base csv file
    :param step_list: list of dicts with data re convergence
    :return: n/a, save file
    """
    png_out = create_out_fname(out_fname, prefix='', ext='.png')

    png_titles = [CONVERG, ENERGY, MAX_FORCE, RMS_FORCE, MAX_DISPL, RMS_DISPL]
    num_lists = len(png_titles)
    png_lists = [[] for _ in range(num_lists)]
    steps = []

    for s_dict in step_list:
        steps.append(s_dict[STEP_NUM])
        for list_id in range(num_lists):
            png_lists[list_id].append(s_dict[png_titles[list_id]])
    fig, axs = plt.subplots(num_lists, figsize=(7, 11.5))
    for list_id in range(num_lists):
        axs[list_id].plot(steps, png_lists[list_id])
        axs[list_id].set_title(png_titles[list_id])
    plt.subplots_adjust(hspace=0.4)
    plt.xlabel("Step number")
    plt.savefig(png_out, transparent=True, bbox_inches='tight',)
    plt.close()
    print(f"Wrote file: {os.path.relpath(png_out)}")


def check_convergence(check_file_list, step_converg, last_step, best_conv, all_steps_to_stdout):
    """
    Reads a Gaussian output file to check convergence
    :param all_steps_to_stdout: Boolean to print convergence to standard out
    :param check_file_list: list of file names
    :param step_converg: boolean; if True, capture convergence of each step. If false, only the final convergence.
    :param last_step: None or int; if int, the last step number to check for convergence
    :param best_conv: Boolean; if true, print ten steps with the best convergence
    :return: nothing: either saves a file or prints to stdout
    """
    fname_str_length = 36
    conv_str_length = 11
    for fname in check_file_list:
        if len(os.path.basename(fname)) > fname_str_length:
            fname_str_length = len(os.path.basename(fname))

    print(f"{F_NAME:{fname_str_length}} {CONVERG:{conv_str_length}} {CONVERG_ERR}")
    if step_converg:
        headers = STEP_CONVERG_HEADERS
    else:
        headers = FINAL_CONVERG_HEADERS
    for fname in check_file_list:
        log_content = process_gausslog_file(fname, find_converg=True, find_step_converg=step_converg,
                                            last_step_to_read=last_step)
        log_content[F_NAME] = os.path.basename(fname)
        if step_converg:
            # all_steps_to_stdout doesn't need an out_fname, but doesn't hurt either
            if last_step:
                out_fname = sys.stdout
            else:
                out_fname = create_out_fname(fname, prefix='', suffix='_conv_steps', ext='.csv')

            # create list of dicts for each step, for all step_converg options
            step_list = []
            for step_num in log_content[CONVERG_STEP_DICT].keys():
                # not sure necessary to make this new dict, but it is fast enough and clearer for next steps
                step_list.append({F_NAME: log_content[F_NAME], STEP_NUM: step_num,
                                  ENERGY: log_content[CONVERG_STEP_DICT][step_num][ENERGY],
                                  MAX_FORCE: log_content[CONVERG_STEP_DICT][step_num][MAX_FORCE],
                                  RMS_FORCE: log_content[CONVERG_STEP_DICT][step_num][RMS_FORCE],
                                  MAX_DISPL: log_content[CONVERG_STEP_DICT][step_num][MAX_DISPL],
                                  RMS_DISPL: log_content[CONVERG_STEP_DICT][step_num][RMS_DISPL],
                                  CONVERG: log_content[CONVERG_STEP_DICT][step_num][CONVERG],
                                  CONVERG_ERR: log_content[CONVERG_STEP_DICT][step_num][CONVERG_ERR],
                                  })

            # different output depending on which step_converg option
            if last_step or best_conv:
                if len(step_list) == 0:
                    print("No convergence data found for file: {}".format(log_content[F_NAME]))
                    continue
                sorted_by_converg = sorted(step_list, key=itemgetter(CONVERG))
                if last_step:
                    print("Steps sorted by convergence to step number {} for file: {}".format(last_step,
                                                                                              log_content[F_NAME]))
                    stop_step = last_step
                else:
                    print("Best (up to 10) steps sorted by convergence for file: {}".format(log_content[F_NAME]))
                    stop_step = 10
                print("    StepNum  Convergence")
                for print_num, step_dict in enumerate(sorted_by_converg):
                    if print_num == stop_step:
                        # break this for, and go to next file if there is one
                        break
                    print("    {:7} {:10.3f}".format(step_dict[STEP_NUM], step_dict[CONVERG]))
            elif all_steps_to_stdout:
                # print all steps to stdout, not sorted by convergence
                print("Convergence of all steps for file: {}".format(log_content[F_NAME]))
                print("    StepNum  Convergence")
                for step_dict in step_list:
                    print("    {:7} {:10.3f}".format(step_dict[STEP_NUM], step_dict[CONVERG]))
            else:
                # save all steps, not sorted by convergence
                print(f"{log_content[F_NAME]:{fname_str_length}} {step_list[-1][CONVERG]:{conv_str_length}.4f} "
                      f"{step_list[-1][CONVERG_ERR]}")
                write_csv(step_list, out_fname, headers, extrasaction="ignore", round_digits=6)
                # also make plots of step versus convergence
                create_convergence_plots(out_fname, step_list)
        else:
            # this is the printing for final termination step only (not step_converg)
            fname = log_content[headers[0]]
            print(f"{fname:{fname_str_length}} {log_content[headers[1]]:{conv_str_length}.4f} "
                  f"{log_content[headers[2]]}")


def process_scan_array(scan_array):
    """
    Script to find the first pair of values in 1st column of a numpy array and process scan arrays
    in degrees to remove jump in periodicity
    :param scan_array: a numpy array with at least two rows and 1 column
    :return: difference in values, accounting for possible angle periodicity in degrees
    """
    prev_val = scan_array[0][0]
    # next line to make ide happy (removes chance trying to return an undefined value)
    first_vals_diff = scan_array[1][0] - prev_val

    for idx in range(1, len(scan_array)):
        current_val = scan_array[idx][0]
        val_diff = current_val - prev_val
        if val_diff > 300.:
            val_diff -= 360.
            current_val -= 360.
        elif val_diff < -300:
            val_diff += 360.
            current_val += 360.
        if idx == 1:
            first_vals_diff = val_diff
        # will catch case of small negative number
        if abs(current_val) < 0.002:
            current_val = 0.0
        scan_array[idx][0] = current_val
        prev_val = current_val
    return first_vals_diff


def collect_output_scan_steps(check_file_list):
    """
    Looks for scan values in one or more files.
    Current functionality: returns one scan, or combines two scans if they search in opposite directions
    :param check_file_list:
    :return: a 2D numpy array with the scan values and energy differences in kcal/mol
    """
    scan_arrays = []
    for fname in check_file_list:
        log_content = process_gausslog_file(fname, collect_scan_steps=True)
        if len(log_content[SCAN_DICT]) > 0:
            scan_arrays.append(np.array(list(log_content[SCAN_DICT].items()), dtype=float))
    num_arrays = len(scan_arrays)
    # if only one scan file, return it
    if num_arrays == 1:
        return_array = scan_arrays[0]
    elif num_arrays == 0:
        raise InvalidDataError("No scan information found.")
    elif num_arrays == 2:
        first_array = scan_arrays[0]
        second_array = scan_arrays[1]
        first_diff = process_scan_array(first_array)
        second_diff = process_scan_array(second_array)
        # check if the first entry is in common, as for scan in two directions
        if abs(first_array[0][0] - second_array[0][0]) < 0.002:
            if first_diff < 0 < second_diff:
                first_array = np.flip(first_array, 0)
                return_array = np.vstack((first_array[:-1, :], second_array))
            elif first_diff > 0 > second_diff:
                second_array = np.flip(second_array, 0)
                return_array = np.vstack((second_array[:-1, :], first_array))
            else:
                raise InvalidDataError("Check how the scans are to be combined.")
        else:
            raise InvalidDataError("The program cannot currently handle these files. Check input, and if correct, "
                                   "please open an issue on github.")
    # convert dict to array
    else:
        raise InvalidDataError("The program can't yet handle this number of files. Please open an issue.")
    # find lowest energy and convert to differences in kcal/mol
    min_e = np.min(return_array[:, 1])
    return_array[:, 1] = (return_array[:, 1] - min_e) * EHPART_TO_KCAL_MOL
    return return_array


def charmm_dihedral(phi, k1, k2, k3, k4, k6, d1, d2, d3, d4, d6, n1, n2, n3, n4, n6):
    """
    Formula for representing a dihedral angle, where k and delta are vectors of the same dimension
    as the n vector below; see https://pubs.acs.org/doi/10.1021/ci500112w
    :param phi: float or np array, the dihedral angle(s) in degrees
    :param k1: float, param
    :param k2: float, param
    :param k3: float, param
    :param k4: float, param
    :param k6: float, param
    :param d1: float, param
    :param d2: float, param
    :param d3: float, param
    :param d4: float, param
    :param d6: float, param
    :param n1: int (0 or 1), indicates if n term is to be used
    :param n2: int (0 or 1), indicates if n term is to be used
    :param n3: int (0 or 1), indicates if n term is to be used
    :param n4: int (0 or 1), indicates if n term is to be used
    :param n6: int (0 or 1), indicates if n term is to be used
    :return: float, result
    """
    k_vec = np.asarray([k1, k2, k3, k4, k6])
    d_vec = np.asarray([d1, d2, d3, d4, d6])
    n_vec = np.asarray([n1, n2, n3, n4, n6])
    result_vec_list = []
    for n, k, d, n_use in zip(N_DIHE, k_vec, d_vec, n_vec):
        if n_use:
            result_vec_list.append(k * (1. + np.cos(n * np.deg2rad(phi) - d)))
    return np.add.reduce(result_vec_list)


def find_good_fit(x_vals, y_vals, x_fit, png_fname=None):
    """
    Find a good functional fit for scan data
    :param x_vals: np array, x values for fitting
    :param y_vals: np array, y values for fitting
    :param x_fit: np array, x values to use for creating curve
    :param png_fname: str, path to save plot, if desired
    :return:
    """
    smallest_resid = np.inf
    best_y_fit = None

    print("Residuals from curve fitting:")

    charmm_n_multipliers = [np.ones(5, dtype=int), np.asarray([0, 1, 1, 1, 1]),
                            np.asarray([1, 1, 1, 1, 0]), np.asarray([1, 1, 1, 0, 0]),
                            np.asarray([1, 1, 0, 0, 0]), np.asarray([1, 0, 0, 0, 0]),
                            np.asarray([0, 1, 0, 1, 1]), np.asarray([1, 0, 1, 0, 0])]
    if png_fname:
        plt.plot(x_vals, y_vals, '.', label='data')

    for idx, multipliers in enumerate(charmm_n_multipliers):
        n_vals = multipliers * N_DIHE
        # fit curve
        ini_vals = np.ones(len(N_DIHE) * 2)
        with warnings.catch_warnings():
            warnings.simplefilter("error", OptimizeWarning)
            try:
                popt, pcov = curve_fit(f=lambda x, *params: charmm_dihedral(x, *params, *multipliers),
                                       xdata=x_vals, ydata=y_vals, p0=ini_vals)
            except OptimizeWarning:
                pass

        y_fit = charmm_dihedral(x_fit, *popt, *multipliers)
        if png_fname:
            plt.plot(x_fit, y_fit, '-', color=assign_color(idx), label=f'fit: {multipliers}')

        y_from_fit = charmm_dihedral(x_vals, *popt, *multipliers)
        resid = np.sqrt(np.mean(np.square(y_from_fit - y_vals)))  # Root Mean Squared Error
        print(f'    CHARMM dihedral eq with n = {",".join([str(x) for x in n_vals[n_vals != 0]]) + ":":10} '
              f'{resid:5.2f}')
        if resid < smallest_resid:
            smallest_resid = resid
            best_y_fit = y_fit

    if png_fname:
        # plt.legend()
        charmm_fname = create_out_fname(png_fname, suffix="_charmm")
        plt.savefig(charmm_fname, transparent=True, bbox_inches='tight',)
        plt.close()
        print(f"Saved: {charmm_fname}")

    if png_fname:
        plt.plot(x_vals, y_vals, '.', label='data')

    for idx, order in enumerate(range(1, 12)):
        # noinspection PyTupleAssignmentBalance
        p, residuals, rank, singular_values, rcond = np.polyfit(x_vals, y_vals, order, full=True)
        y_fit = np.polyval(p, x_fit)

        if png_fname:
            plt.plot(x_fit, y_fit, '-', color=COLOR_SEQUENCE[idx], label=f'fit: poly order {order}')

        y_from_fit = np.polyval(p, x_vals)
        resid = np.sqrt(np.mean(np.square(y_from_fit - y_vals)))
        print(f'    Polynomial order {order:2}: {resid:5.2f}')

        if resid < smallest_resid:
            smallest_resid = resid
            best_y_fit = y_fit

    if png_fname:
        # plt.legend()
        poly_fname = create_out_fname(png_fname, suffix="_poly")
        plt.savefig(poly_fname, transparent=True, bbox_inches='tight',)
        plt.close()
        print(f"Saved: {poly_fname}")

    return best_y_fit


def plot_scan(scan_array, png_out_fname):
    """
    Given a 2D array of x and y values, plot with fitting for a dihedral scan
    :param scan_array: np array of shape [2, n]
    :param png_out_fname: str, path
    :return: n/a, saves plot
    """
    x_vals = scan_array[:, 0]
    y_vals = scan_array[:, 1]
    plt.plot(x_vals, y_vals, '.')

    b_spline = make_interp_spline(x_vals, y_vals)
    x_fit = np.linspace(min(x_vals), max(x_vals), 50)
    y_fit = b_spline(x_fit)
    plt.plot(x_fit, y_fit, '-', color=COLOR_SEQUENCE[0])

    axis_font_size = 14
    plt.xlabel("Dihedral Scan (degrees)", fontsize=axis_font_size)
    plt.ylabel("Change in Energy (kcal/mol)", fontsize=axis_font_size)

    plt.savefig(png_out_fname, transparent=True, bbox_inches='tight',)
    plt.close()
    print(f"Wrote file: {os.path.relpath(png_out_fname)}")

    return x_fit, y_fit


def find_stable_points(x_fit, y_fit):
    """
    Find x values for local min and max of y-values, and return differences in heights;
    Output formatting assumes the y-values correspond to energy in kcal/mol, and x-values are degrees
    :param x_fit: np array, e.g. values from dihedral scan in degrees
    :param y_fit: np array, e.g. energies from a dihedral scan in kcal/mol
    :return: n/a, prints results to stdout
    """
    # argrelmin, argrelmax, and argrelextrema return tuples, even for 1-D data
    local_max_idxs = argrelmin(y_fit)[0]
    local_min_idxs = argrelmax(y_fit)[0]
    stable_points = np.sort(np.hstack((local_max_idxs, local_min_idxs)))
    print("Barriers in kcal/mol:")
    previous_idx = None
    for idx in [0] + list(stable_points) + [-1]:
        if previous_idx is not None:
            print(f"{x_fit[previous_idx]:8.1f} to {x_fit[idx]:5.1f} degrees: "
                  f"{abs(y_fit[idx] - y_fit[previous_idx]):5.1f} kcal/mol")
        previous_idx = idx


def main(argv=None):
    print(f"Running GaussianWrangler script check_gauss version {__version__}")

    args, ret = parse_cmdline(argv)
    if ret != GOOD_RET or args is None:
        return ret

    try:
        # Find files to process, then process them
        check_sub_dirs = False
        search_dir = None
        if args.dir_subdirs:
            search_dir = args.dir_subdirs
            check_sub_dirs = True
        elif args.directory:
            search_dir = args.directory
        check_file_list = check_for_files(args.file_name, args.file_list, search_pattern=args.extension,
                                          search_dir=search_dir, search_sub_dir=check_sub_dirs)

        # now check either for convergence or termination
        if args.step_converg or args.final_converg:
            check_convergence(check_file_list, args.step_converg, args.to_step, args.best, args.all)
        else:
            # If output directory does not exist, make it:
            if not os.path.exists(args.output_directory):
                os.makedirs(args.output_directory)
            if args.scan:
                scan_array = collect_output_scan_steps(check_file_list)
                x_fit, y_fit = plot_scan(scan_array, args.scan)
                find_stable_points(x_fit, y_fit)
            else:
                check_termination(args, check_file_list)

    except IOError as e:
        warning("Problems reading file:", e)
        return IO_ERROR
    except InvalidDataError as e:
        warning("", e)
        return INVALID_DATA

    return GOOD_RET  # success


if __name__ == '__main__':
    status = main()
    sys.exit(status)
