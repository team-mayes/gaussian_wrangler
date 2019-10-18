#!/usr/bin/env python
"""
Creates Gaussian files from pdb files.

Provide template input file
List of pdb files to convert; there may be multiple structures in one file; each one must end with "END/n"
"""

from __future__ import print_function

import os
import sys
import argparse
from common_wrangler.common import (list_to_file, InvalidDataError, create_out_fname, process_cfg, warning,
                                    GOOD_RET, INPUT_ERROR, IO_ERROR, INVALID_DATA, PDB_LINE_TYPE_LAST_CHAR,
                                    PDB_MOL_NUM_LAST_CHAR, PDB_Z_LAST_CHAR, PDB_BEFORE_ELE_LAST_CHAR,
                                    PDB_ELE_LAST_CHAR, PDB_ATOM_NUM_LAST_CHAR, PDB_ATOM_TYPE_LAST_CHAR,
                                    MAIN_SEC, SEC_HEAD, SEC_TAIL)
from gaussian_wrangler.gw_common import (process_gausscom_file)

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
GAU_TPL_FILE = 'gau_tpl_file'
PDB_LIST_FILE = 'pdb_list_file'
PDB_FILE = 'pdb_file'
REMOVE_H = 'remove_final_h'
FIRST_ONLY = 'first_only'

# Defaults
DEF_CFG_FILE = 'pdb2gau.ini'
DEF_CFG_VALS = {PDB_LIST_FILE: 'pdb_list.txt',
                PDB_FILE: None,
                REMOVE_H: False,
                FIRST_ONLY: False,
                }
REQ_KEYS = {GAU_TPL_FILE: str,
            }

CONFIG_NAME = 'config_fname'
# From data template file
NUM_ATOMS = 'num_atoms'
ATOM_TYPE_DICT = 'atom_type_dict'


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
    if good_files:
        main_proc = cfg_proc(dict(config.items(MAIN_SEC)), def_cfg_vals=DEF_CFG_VALS, req_keys=REQ_KEYS)
    else:
        main_proc = {GAU_TPL_FILE: None, CONFIG_NAME: floc}
        for key, def_val in DEF_CFG_VALS.items():
            main_proc[key] = def_val
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
    parser.add_argument("-c", "--config", help="Optional: the location of the configuration file. The default file "
                                               "name is '{}', located in the base directory where the program as run. "
                                               "If a config file is not provided, use the command-line options to "
                                               "specify the '{}' (-t) and '{}' (-1) or '{}' (-p). The command lines "
                                               "for the '{}' flag (-r) or only the first entry in the pdb ('{}', -a) "
                                               "may also be specified.".format(DEF_CFG_FILE, GAU_TPL_FILE,
                                                                               PDB_LIST_FILE, PDB_FILE,
                                                                               REMOVE_H, FIRST_ONLY),
                        default=DEF_CFG_FILE, type=read_cfg)
    parser.add_argument("-t", "--tpl_file", help="Only if a config file is not provided, this command is required to "
                                                 "specify the '{}'".format(GAU_TPL_FILE),
                        default=None)
    parser.add_argument("-l", "--pdb_list_file", help="Only if a config file is not provided, this command can be used "
                                                      "to specify a file with a list of pdbs ('{}') to convert (one "
                                                      "file per line on the list).".format(PDB_LIST_FILE),
                        default=None)
    parser.add_argument("-p", "--pdb_file", help="Only if a config file is not provided, this command can be used to "
                                                 "specify a pdb file ('{}') to convert.".format(PDB_FILE),
                        default=None)
    parser.add_argument("-r", "--remove_final_h", help="Only if a config file is not provided, this command can be "
                                                       "used to specify removing the last H atom from the PDB file(s) "
                                                       "when creating the gausscom files. The default is False.",
                        action='store_true')
    parser.add_argument("-a", "--first_only", help="Only read if a config file is not provided. This command can be "
                                                   "used to specify only using the first set of coordinates in a pdb "
                                                   "file to create gausscom file(s). The default is False.",
                        action='store_true')

    args = None
    try:
        args = parser.parse_args(argv)
        if args.config[GAU_TPL_FILE] is None:
            if args.tpl_file is None:
                raise InvalidDataError("Could not read config file: {}\n    and did not specify a 'tpl_file' "
                                       "('-t' option). A tpl_file is needed to run this "
                                       "script.".format(args.config[CONFIG_NAME]))
            else:
                args.config[GAU_TPL_FILE] = args.tpl_file
                if args.first_only:
                    args.config[FIRST_ONLY] = True
                if args.remove_final_h:
                    args.config[REMOVE_H] = True
                if args.pdb_file:
                    args.config[PDB_FILE] = args.pdb_file
                if args.pdb_list_file:
                    args.config[PDB_LIST_FILE] = args.pdb_list_file
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


def process_pdb_files(cfg, gau_tpl_content):
    pdb_files = []
    if cfg[PDB_FILE]:
        if os.path.isfile(cfg[PDB_FILE]):
            pdb_files.append(cfg[PDB_FILE])
        else:
            raise IOError(cfg[PDB_FILE])
    if os.path.isfile(cfg[PDB_LIST_FILE]):
        with open(cfg[PDB_LIST_FILE]) as f:
            for pdb_file in f.readlines():
                pdb_file = pdb_file.strip()
                if len(pdb_file) > 0:
                    pdb_files.append(pdb_file)
    if len(pdb_files) == 0:
        raise InvalidDataError("No pdb files found to process.")
    for pdb_file in pdb_files:
        process_pdb_file(cfg, gau_tpl_content, pdb_file)


def process_pdb_file(cfg, gau_tpl_content, pdb_file):
    with open(pdb_file) as d:
        mol_num = 0
        pdb_atom_line = []
        for line in d.readlines():
            pdb_section = line[:PDB_LINE_TYPE_LAST_CHAR]
            if pdb_section == 'MODEL ':
                mol_num += 1
            elif pdb_section == 'ATOM  ' or pdb_section == 'HETATM':
                element = line[PDB_BEFORE_ELE_LAST_CHAR:PDB_ELE_LAST_CHAR].strip()
                if element == '':
                    element = line[PDB_ATOM_NUM_LAST_CHAR:PDB_ATOM_TYPE_LAST_CHAR].strip()
                pdb_xyz = line[PDB_MOL_NUM_LAST_CHAR:PDB_Z_LAST_CHAR]
                pdb_atom_line.append(["{:6}".format(element), pdb_xyz])
            elif pdb_section == 'END\n':
                if mol_num == 0:
                    mol_id = ''
                else:
                    mol_id = '_' + str(mol_num)
                d_out = create_out_fname(pdb_file, suffix=mol_id, ext='.com')
                if cfg[REMOVE_H]:
                    del pdb_atom_line[-1]
                list_to_file(gau_tpl_content[SEC_HEAD] + pdb_atom_line + gau_tpl_content[SEC_TAIL],
                             d_out)
                if cfg[FIRST_ONLY]:
                    return
                pdb_atom_line = []


def main(argv=None):
    # Read input
    args, ret = parse_cmdline(argv)
    if ret != GOOD_RET or args is None:
        return ret

    # Read pdb files
    cfg = args.config
    try:
        gau_tpl_content = process_gausscom_file(cfg[GAU_TPL_FILE])
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
