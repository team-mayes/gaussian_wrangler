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
from common_wrangler.common import (InvalidDataError, warning, GOOD_RET, INPUT_ERROR, IO_ERROR, INVALID_DATA)

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
FAIL_FILE_LEN = re.compile(r"File lengths (MBytes): RWF=*")
FAIL_RDCARD = re.compile(r"In source file rdcard*")
FAIL_NTR = re.compile(r"NtrErr Called from FileIO*")
FAIL_PAT_LIST = [FAIL_OPEN_FILE, FAIL_FILE_LEN, FAIL_RDCARD, FAIL_NTR]

DEF_COMPLETE_DIR = 'for_hartree'
DEF_EXT = '.log'
OUT_BASE_DIR = 'output_directory'


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
    parser.add_argument("-d", "--directory", help="The directory where to look for Gaussian output to check for "
                                                  "normal termination. The default is the current directory.",
                        default=None)
    parser.add_argument("-e", "--extension", help="The extension of the Gaussian output file. The default is '{}'."
                                                  "".format(DEF_EXT), default=DEF_EXT)
    parser.add_argument("-o", "--output_directory", help="The directory where to put Gaussian output files that have "
                                                         "terminated normally. The default is '{}'."
                                                         "".format(DEF_COMPLETE_DIR), default=DEF_COMPLETE_DIR)
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


def main(argv=None):
    # Read input
    args, ret = parse_cmdline(argv)
    if ret != GOOD_RET or args is None:
        return ret

    # Find files to process, then process them
    completed_list = []
    check_file_list = []
    perhaps_running_list = []
    likely_failed_list = []
    try:
        # check input for specified search directory (if specified)
        if args.directory:
            if os.path.isdir(args.directory):
                search_folder = args.directory
            else:
                raise InvalidDataError("Could not find the specified input directory '{}'".format(args.directory))
        else:
            search_folder = os.getcwd()
        # If output directory does not exist, make it:
        if not os.path.exists(args.output_directory):
            os.makedirs(args.output_directory)
        # search for output files
        for file in os.listdir(search_folder):
            if file.endswith(args.extension):
                check_file_list.append(os.path.join(search_folder, file))
        if len(check_file_list) == 0:
            raise InvalidDataError("Could not find files with extension '{}' in directory '{}'".format(args.extension,
                                                                                                       search_folder))
        for file in check_file_list:
            process_list_file(file, args.output_directory, completed_list, likely_failed_list, perhaps_running_list)
        # sort if list is at least 2 long:
        for file_list in [completed_list, likely_failed_list, perhaps_running_list]:
            if len(file_list) > 1:
                file_list.sort()
        if len(completed_list) > 0:
            print("The following files completed normally:")
            for file in completed_list:
                print("    {}".format(os.path.relpath(file)))
        else:
            print("No normally completed files found.")
        if len(likely_failed_list) > 0:
            print("The following files may have failed:")
            for file in likely_failed_list:
                print("    {}".format(os.path.relpath(file)))
        if len(perhaps_running_list) > 0:
            print("The following files may still be running:")
            for file in perhaps_running_list:
                print("    {}".format(os.path.relpath(file)))
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
