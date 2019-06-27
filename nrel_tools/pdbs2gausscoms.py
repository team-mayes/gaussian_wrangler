#!/usr/bin/env python
"""
Creates Gaussian files from pdb files.

Provide template input file
List of pdb files to convert; there may be multiple structures in one file; each one must end with "END/n"
"""

from __future__ import print_function
import sys
import argparse
from nrel_tools.common import (list_to_file, InvalidDataError, create_out_fname, process_cfg, warning,
                               GOOD_RET, INPUT_ERROR, IO_ERROR, INVALID_DATA)
try:
    # noinspection PyCompatibility
    from ConfigParser import ConfigParser
except ImportError:
    # noinspection PyCompatibility
    from configparser import ConfigParser

__author__ = 'hmayes'

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
GAU_TPL_FILE = 'gau_tpl_file'
PDBS_FILE = 'pdb_list_file'
REMOVE_H = 'remove_final_h'
# PDB file info
PDB_SECTION_LAST_CHAR = 'pdb_section_last_char'
PDB_ATOM_NUM_LAST_CHAR = 'pdb_atom_num_last_char'
PDB_ATOM_INFO_LAST_CHAR = 'pdb_atom_info_last_char'
PDB_MOL_NUM_LAST_CHAR = 'pdb_mol_num_last_char'
PDB_X_LAST_CHAR = 'pdb_x_last_char'
PDB_Y_LAST_CHAR = 'pdb_y_last_char'
PDB_Z_LAST_CHAR = 'pdb_z_last_char'
PDB_BEFORE_ELE_LAST_CHAR = 'pdb_before_element_last_char'
PDB_ELE_LAST_CHAR = 'pdb_element_last_char'
PDB_FORMAT = 'pdb_print_format'

# Defaults
DEF_CFG_FILE = 'pdb2gau.ini'
DEF_CFG_VALS = {PDBS_FILE: 'pdb_list.txt',
                REMOVE_H: False,
                PDB_SECTION_LAST_CHAR: 6,
                PDB_ATOM_NUM_LAST_CHAR: 12,
                PDB_ATOM_INFO_LAST_CHAR: 22,
                PDB_MOL_NUM_LAST_CHAR: 28,
                PDB_X_LAST_CHAR: 38,
                PDB_Y_LAST_CHAR: 46,
                PDB_Z_LAST_CHAR: 54,
                PDB_BEFORE_ELE_LAST_CHAR: 76,
                PDB_ELE_LAST_CHAR: 78,
                }
REQ_KEYS = {GAU_TPL_FILE: str,
            }

# From data template file
NUM_ATOMS = 'num_atoms'
HEAD_CONTENT = 'head_content'
ATOMS_CONTENT = 'atoms_content'
TAIL_CONTENT = 'tail_content'
ATOM_TYPE_DICT = 'atom_type_dict'

# For data template file processing
SEC_HEAD = 'head_section'
SEC_ATOMS = 'atoms_section'
SEC_TAIL = 'tail_section'


def read_cfg(floc, cfg_proc=process_cfg):
    """
    Reads the given configuration file, returning a dict with the converted values supplemented by default values.

    :param floc: The location of the file to read.
    :param cfg_proc: The processor to use for the raw configuration values.  Uses default values when the raw
        value is missing.
    :return: A dict of the processed configuration file's data.
    """
    config = ConfigParser()
    good_files = config.read(floc)
    if not good_files:
        raise IOError('Could not read file {}'.format(floc))
    main_proc = cfg_proc(dict(config.items(MAIN_SEC)), def_cfg_vals=DEF_CFG_VALS, req_keys=REQ_KEYS)
    return main_proc


def parse_cmdline(argv):
    """
    Returns the parsed argument list and return code.
    `argv` is a list of arguments, or `None` for ``sys.argv[1:]``.
    """
    if argv is None:
        argv = sys.argv[1:]

    # initialize the parser object:
    parser = argparse.ArgumentParser(description='Creates Gaussian input files from pdb files, given a template input '
                                                 'file. The required input file provides the name/location of the '
                                                 'template file and a file with a list of pdb files to convert.')
    parser.add_argument("-c", "--config", help="The location of the configuration file. "
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
    except KeyError as e:
        warning("Input data missing:", e)
        parser.print_help()
        return args, INPUT_ERROR

    return args, GOOD_RET


def process_pdb_files(cfg, gau_tpl_content):
    with open(cfg[PDBS_FILE]) as f:
        for pdb_file in f.readlines():
            pdb_file = pdb_file.strip()
            with open(pdb_file) as d:
                mol_num = 0
                pdb_atom_line = []
                for line in d.readlines():
                    pdb_section = line[:cfg[PDB_SECTION_LAST_CHAR]]
                    if pdb_section == 'MODEL ':
                        mol_num += 1
                    elif pdb_section == 'ATOM  ' or pdb_section == 'HETATM':
                        element = line[cfg[PDB_BEFORE_ELE_LAST_CHAR]:cfg[PDB_ELE_LAST_CHAR]].strip()
                        pdb_xyz = line[cfg[PDB_MOL_NUM_LAST_CHAR]:cfg[PDB_Z_LAST_CHAR]]
                        pdb_atom_line.append(["{:6}".format(element), pdb_xyz])
                    elif pdb_section == 'END\n':
                        if mol_num == 0:
                            mol_id = ''
                        else:
                            mol_id = '_' + str(mol_num)
                        d_out = create_out_fname(pdb_file, suffix=mol_id, ext='.com')
                        if cfg[REMOVE_H]:
                            del pdb_atom_line[-1]
                        list_to_file(gau_tpl_content[HEAD_CONTENT] + pdb_atom_line + gau_tpl_content[TAIL_CONTENT],
                                     d_out)
                        pdb_atom_line = []


def process_gau_tpl(cfg):
    tpl_loc = cfg[GAU_TPL_FILE]
    tpl_data = {HEAD_CONTENT: [], TAIL_CONTENT: []}
    section = SEC_HEAD
    # atoms_sec_pat = re.compile(r"atoms")

    with open(tpl_loc) as f:
        for line in f.readlines():
            line = line.strip()
            # atoms_match = atoms_sec_pat.match(line)
            # if atoms_match:
            if line == '${atoms}':
                section = SEC_TAIL
                continue
            # head_content to contain Everything before 'Atoms' section
            elif section == SEC_HEAD:
                tpl_data[HEAD_CONTENT].append(line)
            elif section == SEC_TAIL:
                tpl_data[TAIL_CONTENT].append(line)
    return tpl_data


def main(argv=None):
    # Read input
    args, ret = parse_cmdline(argv)
    if ret != GOOD_RET:
        return ret

    # Read pdb files
    cfg = args.config
    try:
        gau_tpl_content = process_gau_tpl(cfg)
        process_pdb_files(cfg, gau_tpl_content)
    except IOError as e:
        warning("Problems reading file:", e)
        return IO_ERROR
    except InvalidDataError as e:
        warning("Problems reading data template:", e)
        return INVALID_DATA

    return GOOD_RET  # success


if __name__ == '__main__':
    status = main()
    sys.exit(status)
