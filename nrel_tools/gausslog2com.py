#!/usr/bin/env python
"""
Creates pdb data files from lammps data files, given a template pdb file.
"""

from __future__ import print_function
import os
import re
import sys
import argparse
from nrel_tools.common import (InvalidDataError, warning, create_out_fname, list_to_file, ATOM_NUM_DICT,
                               NUM_ATOMS, GOOD_RET, INPUT_ERROR, IO_ERROR, INVALID_DATA)

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

GAU_COORD_PAT = re.compile(r"Center     Atomic      Atomic             Coordinates.*")
GAU_SEP_PAT = re.compile(r"---------------------------------------------------------------------.*")
GAU_E_PAT = re.compile(r"SCF Done:.*")

# For log file processing
SEC_HEAD = 'head_section'
SEC_ATOMS = 'atoms_section'
SEC_TAIL = 'tail_section'
COM_TYPE = 'com_type'


def parse_cmdline(argv):
    """
    Returns the parsed argument list and return code.
    `argv` is a list of arguments, or `None` for ``sys.argv[1:]``.
    """
    if argv is None:
        argv = sys.argv[1:]

    # initialize the parser object:
    parser = argparse.ArgumentParser(description='Creates Gaussian input files from Gaussian output files, given a '
                                                 'template file.')
    parser.add_argument("-f", "--file", help="The location of the Gaussian output file.")
    parser.add_argument("-l", "--list", help="The location of the list of Gaussian output files. "
                                             "The default file name is {}, located in the "
                                             "base directory where the program as run.".format(DEF_LIST_FILE),
                        default=DEF_LIST_FILE)
    parser.add_argument("-t", "--tpl", help="The location of the Gaussian input template file.")

    args = None
    try:
        args = parser.parse_args(argv)
    except SystemExit as e:
        if hasattr(e, 'code') and e.code == 0:
            return args, GOOD_RET
        warning(e)
        parser.print_help()
        return args, INPUT_ERROR

    return args, GOOD_RET


def process_gausscom_files(gausslog_files, com_tpl_content):
    for gausslog_file in gausslog_files:
        process_gausslog_file(gausslog_file, com_tpl_content)


def process_gausslog_file(gausslog_file, com_tpl_content):
    with open(gausslog_file) as d:
        section = SEC_HEAD
        atom_id = 0
        lines_after_coord = 2
        coord_match = False

        for line in d:
            line = line.strip()
            if len(line) == 0:
                continue
            # not currently keeping anything from the header
            if section == SEC_HEAD:
                if GAU_COORD_PAT.match(line):
                    coord_match = True
                    atoms_section = []
                    continue
                elif coord_match and lines_after_coord > 0:
                    lines_after_coord -= 1
                    if lines_after_coord == 0:
                        section = SEC_ATOMS
                    continue

            elif section == SEC_ATOMS:
                if GAU_SEP_PAT.match(line):
                    section = SEC_TAIL
                    continue

                split_line = line.split()

                try:
                    atom_type = ATOM_NUM_DICT[int(split_line[1])]
                except KeyError:
                    raise InvalidDataError("Currently, this code only expects atom numbers up to 36 (Kr), and the "
                                           "atomic number read was {}. Update the code to use this with your current "
                                           "output.".format(split_line[1]))
                if com_tpl_content[NUM_ATOMS]:
                    com_atom_type = com_tpl_content[SEC_ATOMS][atom_id].split('(')[0]
                    if com_atom_type != atom_type:
                        raise InvalidDataError("For atom number {}, {} has atom type {}, while the template has atom "
                                               "type {}".format(atom_id+1, gausslog_file, atom_type, com_atom_type))
                    atom_type = com_tpl_content[SEC_ATOMS][atom_id]  # This keeps the "fragment" number if there
                atom_type += '      '

                atom_xyz = ["{:>12}".format(x) for x in split_line[3:6]]
                atoms_section.append(atom_type + ''.join(atom_xyz))
                atom_id += 1
            elif section == SEC_TAIL:
                if com_tpl_content[NUM_ATOMS] and atom_id != com_tpl_content[NUM_ATOMS]:
                    raise InvalidDataError('In gausslog file: {}\n  found {} atoms, but the tpl expects '
                                           '{} atoms'.format(gausslog_file, atom_id, com_tpl_content[NUM_ATOMS]))
                if GAU_E_PAT.match(line):
                    section = SEC_HEAD
                    coord_match = False
                    atom_id = 0
                    lines_after_coord = 2

    f_name = create_out_fname(gausslog_file, ext='_' + com_tpl_content[COM_TYPE] + '.com')
    list_to_file(com_tpl_content[SEC_HEAD] + atoms_section + com_tpl_content[SEC_TAIL],
                 f_name)

    # Now that finished reading the file, first make sure didn't  exit before reaching the desired number of atoms


def process_gausscom_tpl(com_tpl_file):
    com_tpl_content = {SEC_HEAD: [], SEC_ATOMS: [], SEC_TAIL: ['', ]}
    section = SEC_HEAD
    lines_after_first_blank = 0
    with open(com_tpl_file) as d:
        file_name = os.path.basename(com_tpl_file)
        com_tpl_content[COM_TYPE] = os.path.splitext(file_name)[0]
        for line in d:
            line = line.strip()
            if section == SEC_HEAD:
                com_tpl_content[SEC_HEAD].append(line)
                if len(line) == 0:
                    lines_after_first_blank += 1
                elif lines_after_first_blank == 2:
                    section = SEC_ATOMS
            elif section == SEC_ATOMS:
                if len(line) == 0:
                    section = SEC_TAIL
                    continue
                line_split = line.split()
                com_tpl_content[SEC_ATOMS].append(line_split[0] + '    ')
            elif section == SEC_TAIL:
                com_tpl_content[SEC_TAIL].append(line)
    if len(com_tpl_content[SEC_ATOMS]) == 0:
        com_tpl_content[NUM_ATOMS] = None
    else:
        com_tpl_content[NUM_ATOMS] = len(com_tpl_content[SEC_ATOMS])

    return com_tpl_content


def main(argv=None):
    # Read input
    args, ret = parse_cmdline(argv)
    if ret != GOOD_RET or args is None:
        return ret

    # Make sure there are log files to process
    gausslog_files = []
    if os.path.isfile(args.list):
        with open(args.list) as f:
            for data_file in f:
                gausslog_files.append(data_file.strip())
    if args.file is not None:
        gausslog_files.append(args.file)
    if len(gausslog_files) == 0:
        raise InvalidDataError("No files to process: no single log file specified and "
                               "no list of files found")

    # Read template and data files
    try:
        com_tpl_content = process_gausscom_tpl(args.tpl)
        process_gausscom_files(gausslog_files, com_tpl_content)
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
