#!/usr/bin/env python
"""
Finds Gaussian files in a folder. If complete (Normal termination), moves them to a subdirectory.
If appears to have failed, add to a text file for investigation.
"""

from __future__ import print_function

import os
import re
import sys
import argparse
from operator import itemgetter

from common_wrangler.common import (InvalidDataError, warning, GOOD_RET, INPUT_ERROR, IO_ERROR, INVALID_DATA,
                                    create_out_fname, write_csv, check_file_and_file_list)

from gaussian_wrangler.gw_common import (MAX_FORCE, RMS_FORCE, MAX_DISPL, RMS_DISPL, CONVERG, CONVERG_ERR,
                                         process_gausslog_file, CONVERG_STEP_DICT)

try:
    # noinspection PyCompatibility
    from ConfigParser import ConfigParser, MissingSectionHeaderError
except ImportError:
    # noinspection PyCompatibility
    from configparser import ConfigParser, MissingSectionHeaderError

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
F_NAME = 'File'
STEP_CONVERG_HEADERS = [F_NAME, STEP_NUM, MAX_FORCE, RMS_FORCE, MAX_DISPL, RMS_DISPL, CONVERG, CONVERG_ERR]
FINAL_CONVERG_HEADERS = [F_NAME, CONVERG, CONVERG_ERR]


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
                        default=None)
    parser.add_argument("-ds", "--dir_subdirs", help="The directory where to look for Gaussian output files to check "
                                                     "for normal termination, including checking in subdirectories.",
                        default=None)
    parser.add_argument("-e", "--extension", help="The extension of the Gaussian output file(s) to look for when "
                                                  "searching a directory for output files. The default is '{}'."
                                                  "".format(DEF_EXT), default=DEF_EXT)
    parser.add_argument("-f", "--file_name", help="A file name (with path, if not the current directory) to check for "
                                                  "either normal termination or convergence. If used, this option "
                                                  "overrides the '-d' option, and no searching for files is "
                                                  "performed.", default=None)
    parser.add_argument("-l", "--file_list", help="A file name (with path, if not the current directory) with a "
                                                  "list of files (also with path, if not the current directory)  "
                                                  "overrides the '-d' option, and no searching for files is to check "
                                                  "for either normal termination or convergence. If used, this "
                                                  "option overrides the '-d' option, and no searching for files is "
                                                  "performed.", default=None)
    parser.add_argument("-o", "--output_directory", help="The directory where to put Gaussian output files that have "
                                                         "terminated normally. The default is '{}'."
                                                         "".format(DEF_COMPLETE_DIR), default=DEF_COMPLETE_DIR)
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
    except (KeyError, InvalidDataError, MissingSectionHeaderError, SystemExit) as e:
        if hasattr(e, 'code') and e.code == 0:
            return args, GOOD_RET
        warning(e)
        parser.print_help()
        return args, INPUT_ERROR

    return args, GOOD_RET


def process_list_file(output_file, good_output_directory, completed_list, likely_failed_list, perhaps_running_list):
    try:
        with open(output_file, 'r') as fh:
            last_line = fh.readlines()[-1].strip()
    except IndexError:
        warning("Could not read the last line (may be blank) of file: {}".format(output_file))
        return
    if NORM_TERM_PAT.match(last_line):
        base_name = os.path.basename(output_file)
        completed_list.append(output_file)
        os.rename(output_file, os.path.join(good_output_directory, base_name))
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
        process_list_file(fname, args.output_directory, completed_list, likely_failed_list, perhaps_running_list)
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
    if step_converg:
        headers = STEP_CONVERG_HEADERS
    else:
        headers = FINAL_CONVERG_HEADERS
        print("{:36} {:11} {:}".format(*headers))
    for fname in check_file_list:
        log_content = process_gausslog_file(fname, find_converg=True, find_step_converg=step_converg,
                                            last_step_to_read=last_step)
        log_content[F_NAME] = os.path.basename(fname)
        if step_converg:
            # all_steps_to_stdout doesn't need an out_fname, but doesn't hurt either
            if last_step:
                out_loc = sys.stdout
            else:
                out_loc = create_out_fname(fname, prefix='', suffix='_conv_steps', ext='.csv')

            # create list of dicts for each step, for all step_converg options
            step_list = []
            for step_num in log_content[CONVERG_STEP_DICT].keys():
                # not sure necessary to make this new dict, but it is fast enough and clearer for next steps
                step_list.append({F_NAME: log_content[F_NAME], STEP_NUM: step_num,
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
                write_csv(step_list, out_loc, headers, extrasaction="ignore", round_digits=6)
        else:
            # this is the printing for final termination step only (not step_converg)
            fname = log_content[headers[0]]
            try:
                print("{:36} {:11.4f} {:}".format(fname, log_content[headers[1]], log_content[headers[2]]))
            except KeyError:
                print("{:36} {:>11} {:}".format(fname, "not found", "n/a"))


def search_dir_fnames(args, check_file_list):
    # set up so if some reason select args.dir_subdirs and args.directory, runs the former option (covers later option)
    if args.dir_subdirs:
        search_folder = args.dir_subdirs
        for root, dirs, files in os.walk(args.dir_subdirs):
            for fname in files:
                if fname.endswith(args.extension):
                    check_file_list.append(os.path.join(root, fname))
    else:
        if args.directory:
            if os.path.isdir(args.directory):
                search_folder = args.directory
            else:
                raise InvalidDataError("Could not find the specified input directory '{}'".format(args.directory))
        else:
            search_folder = os.getcwd()
            # search for output files
        for file in os.listdir(search_folder):
            if file.endswith(args.extension):
                check_file_list.append(os.path.join(search_folder, file))
    if len(check_file_list) == 0:
        raise InvalidDataError("Could not find files with extension '{}' in "
                               "directory '{}'".format(args.extension, search_folder))


def main(argv=None):
    # Read input
    args, ret = parse_cmdline(argv)
    if ret != GOOD_RET or args is None:
        return ret

    # Find files to process, then process them
    check_file_list = []
    try:
        # first make list of files to check, either by a directory search or from a specified file name and/or list,
        #    allowing multiple choices
        if args.file_name or args.file_list:
            check_file_list = check_file_and_file_list(args.file_name, args.file_list)
        # catches dir specified, dir+subdirs, or none of these options, which defaults to searching the current
        #    directory only
        if args.directory or args.dir_subdirs or (not args.file_name and not args.file_list):
            search_dir_fnames(args, check_file_list)
        check_file_list.sort()

        # now check either for convergence or termination
        if args.step_converg or args.final_converg:
            check_convergence(check_file_list, args.step_converg, args.to_step, args.best, args.all)
        else:
            # If output directory does not exist, make it:
            if not os.path.exists(args.output_directory):
                os.makedirs(args.output_directory)
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
