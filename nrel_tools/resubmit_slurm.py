#!/usr/bin/env python
"""
Given a list of Gaussian log file names, checks the last line of each to determine if it immediately failed and
needs to be restarted
"""

from __future__ import print_function
import re
import sys
import argparse
from nrel_tools.common import (InvalidDataError, warning, process_cfg, create_out_fname, list_to_file, process_pdb_tpl,
                               HEAD_CONTENT, ATOMS_CONTENT, TAIL_CONTENT, PDB_FORMAT, NUM_ATOMS,
                               GOOD_RET, INPUT_ERROR, IO_ERROR, INVALID_DATA)

try:
    # noinspection PyCompatibility
    from ConfigParser import ConfigParser, MissingSectionHeaderError
except ImportError:
    # noinspection PyCompatibility
    from configparser import ConfigParser, MissingSectionHeaderError

__author__ = 'hmayes'


# Constants #


FAIL_OPEN_FILE = re.compile(r"open-new-file*")
FAIL_FILE_LEN = re.compile(r"File lengths (MBytes): RWF=*")
FAIL_RDCARD = re.compile(r"In source file rdcard*")
FAIL_NTR = re.compile(r"NtrErr Called from FileIO*")
PAT_LIST = [FAIL_OPEN_FILE, FAIL_FILE_LEN, FAIL_RDCARD, FAIL_NTR]

DEF_LIST_FILE = 'list.txt'


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
    parser.add_argument("-l", "--list_file", help="The list of files to be checked."
                                                  "The default file name is {}, located in the "
                                                  "base directory where the program as run.".format(DEF_LIST_FILE),
                        default=DEF_LIST_FILE)
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


def process_list_file(list_file):
    with open(list_file) as d:
        for list_line in d:
            list_line = list_line.strip()
            if len(list_line) == 0:
                continue
            problem = False
            with open(list_line, 'r') as fh:
                last_line = fh.readlines()[-1].decode()
            for pattern in PAT_LIST:
                if pattern.match(last_line):
                    problem = True
                    continue
            if problem:
                print("Need to restart: {}; last line is: {}".format(list_file, last_line))

    # f_name = create_out_fname(gausscom_file, ext='.pdb', base_dir=cfg[OUT_BASE_DIR])
    # list_to_file(pdb_tpl_content[HEAD_CONTENT] + pdb_data_section + pdb_tpl_content[TAIL_CONTENT],
    #              f_name, list_format=PDB_FORMAT)


def main(argv=None):
    # Read input
    args, ret = parse_cmdline(argv)
    if ret != GOOD_RET or args is None:
        return ret

    # Read template and data files
    try:
        process_list_file(args.list_file)
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
