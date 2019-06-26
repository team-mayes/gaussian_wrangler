#!/usr/bin/env python
"""
Creates pdb data files from lammps data files, given a template pdb file.
"""

from __future__ import print_function
import os
import copy
import logging
import re
import sys
import argparse
from nrel_tools.common import (InvalidDataError, warning, process_cfg, create_out_fname, list_to_file, process_pdb_tpl,
                               HEAD_CONTENT, ATOMS_CONTENT, TAIL_CONTENT, PDB_FORMAT, NUM_ATOMS)

try:
    # noinspection PyCompatibility
    from ConfigParser import ConfigParser, MissingSectionHeaderError
except ImportError:
    # noinspection PyCompatibility
    from configparser import ConfigParser, MissingSectionHeaderError

__author__ = 'hmayes'


# Logging
logger = logging.getLogger('gausscom2pdb')
# logging.basicConfig(filename='gausscom2pdb.log', filemode='w', level=logging.DEBUG)
logging.basicConfig(level=logging.INFO)

# Error Codes
# The good status code
GOOD_RET = 0
INPUT_ERROR = 1
IO_ERROR = 2
INVALID_DATA = 3

# Constants #

# Config File Sections
MAIN_SEC = 'main'

# Config keys
PDB_TPL_FILE = 'pdb_tpl_file'
GAUSSCOM_FILES_FILE = 'gausscom_list_file'
GAUSSCOM_FILES = 'gausscom_files_list'
GAUSSCOM_FILE = 'gausscom_file'
CENTER_ATOM = 'center_to_atom_num'

OUT_BASE_DIR = 'output_directory'
MAKE_DICT_BOOL = 'make_dictionary_flag'

GAU_HEADER_PAT = re.compile(r"#.*")

# data file info

# Defaults
DEF_CFG_FILE = 'gausscom2pdb.ini'
# Set notation
DEF_CFG_VALS = {GAUSSCOM_FILES_FILE: 'gausscom_list.txt',
                OUT_BASE_DIR: None,
                MAKE_DICT_BOOL: False,
                GAUSSCOM_FILE: None,
                }
REQ_KEYS = {PDB_TPL_FILE: str,
            }

# For pdb template file processing
SEC_HEAD = 'head_section'
SEC_ATOMS = 'atoms_section'
SEC_TAIL = 'tail_section'


def read_cfg(f_loc, cfg_proc=process_cfg):
    """
    Reads the given configuration file, returning a dict with the converted values supplemented by default values.

    :param f_loc: The location of the file to read.
    :param cfg_proc: The processor to use for the raw configuration values.  Uses default values when the raw
        value is missing.
    :return: A dict of the processed configuration file's data.
    """
    config = ConfigParser()
    good_files = config.read(f_loc)

    if not good_files:
        raise IOError('Could not read file {}'.format(f_loc))
    main_proc = cfg_proc(dict(config.items(MAIN_SEC)), DEF_CFG_VALS, REQ_KEYS)

    # To fix; have this as default!
    main_proc[GAUSSCOM_FILES] = []
    if os.path.isfile(main_proc[GAUSSCOM_FILES_FILE]):
        with open(main_proc[GAUSSCOM_FILES_FILE]) as f:
            for data_file in f:
                main_proc[GAUSSCOM_FILES].append(data_file.strip())
    if main_proc[GAUSSCOM_FILE] is not None:
        main_proc[GAUSSCOM_FILES].append(main_proc[GAUSSCOM_FILE])
    if len(main_proc[GAUSSCOM_FILES]) == 0:
        raise InvalidDataError("No files to process: no '{}' specified and "
                               "no list of files found for: {}".format(GAUSSCOM_FILE, main_proc[GAUSSCOM_FILES_FILE]))

    return main_proc


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
    parser.add_argument("-c", "--config", help="The location of the configuration file in ini format. "
                                               "The default file name is {}, located in the "
                                               "base directory where the program as run.".format(DEF_CFG_FILE),
                        default=DEF_CFG_FILE, type=read_cfg)
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


def process_gausscom_files(cfg, pdb_tpl_content):
    # Don't want to change the original template data when preparing to print the new file:

    for gausscom_file in cfg[GAUSSCOM_FILES]:
        process_gausscom_file(cfg, gausscom_file, pdb_tpl_content)


def process_gausscom_file(cfg, gausscom_file, pdb_tpl_content):
    with open(gausscom_file) as d:
        pdb_data_section = copy.deepcopy(pdb_tpl_content[ATOMS_CONTENT])
        section = SEC_HEAD
        atom_id = 0
        lines_after_header = 4  # blank line, description, blank line, charge & multiplicity
        atom_types = []

        for line in d:
            line = line.strip()
            # not currently keeping anything from the header; just check num atoms
            if section == SEC_HEAD:
                if GAU_HEADER_PAT.match(line):
                    continue
                elif lines_after_header > 0:
                    lines_after_header -= 1
                    if lines_after_header == 0:
                        section = SEC_ATOMS
                    continue

            elif section == SEC_ATOMS:
                if len(line) == 0:
                    continue
                split_line = line.split()

                # check atom type
                atom_type = split_line[0]
                pdb_atom_type = pdb_data_section[atom_id][8].split(' ')[-1]
                if atom_type != pdb_atom_type:
                    warning("Atom types do not match for atom number {}; pdb atom type is {} while gausscom type is "
                            "{}".format(atom_id, pdb_atom_type, atom_type))
                # Keep as string; json save as string and this helps compare
                atom_types.append(atom_type)
                pdb_data_section[atom_id][5:8] = map(float, split_line[1:4])
                atom_id += 1
                # Check after increment because the counter started at 0
                if atom_id == pdb_tpl_content[NUM_ATOMS]:
                    # Since the tail will come only from the template, nothing more is needed.
                    break

    # Now that finished reading the file, first make sure didn't  exit before reaching the desired number of atoms
    if atom_id != pdb_tpl_content[NUM_ATOMS]:
        raise InvalidDataError('In gausscom file: {}\n'
                               '  found {} atoms, but pdb expects {} atoms'.format(gausscom_file, atom_id,
                                                                                   pdb_tpl_content[NUM_ATOMS]))
    f_name = create_out_fname(gausscom_file, ext='.pdb', base_dir=cfg[OUT_BASE_DIR])
    list_to_file(pdb_tpl_content[HEAD_CONTENT] + pdb_data_section + pdb_tpl_content[TAIL_CONTENT],
                 f_name, list_format=PDB_FORMAT)


def main(argv=None):
    # Read input
    args, ret = parse_cmdline(argv)
    if ret != GOOD_RET or args is None:
        return ret

    cfg = args.config

    # Read template and data files
    try:
        pdb_tpl_content = process_pdb_tpl(cfg[PDB_TPL_FILE])
        process_gausscom_files(cfg, pdb_tpl_content)
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
