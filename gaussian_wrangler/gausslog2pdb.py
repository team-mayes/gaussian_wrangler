#!/usr/bin/env python
"""
Creates pdb files from Gaussian log files
"""

import os
import copy
import sys
import argparse
from configparser import ConfigParser, MissingSectionHeaderError
from common_wrangler.common import (MAIN_SEC, SEC_HEAD, SEC_ATOMS, SEC_TAIL, PDB_FORMAT, NUM_ATOMS, ATOM_NUM_DICT,
                                    GOOD_RET, INPUT_ERROR, IO_ERROR, INVALID_DATA, InvalidDataError, warning,
                                    process_cfg, create_out_fname, list_to_file, process_pdb_file, silent_remove)
from gaussian_wrangler.gw_common import (GAU_COORD_PAT, GAU_SEP_PAT, GAU_E_PAT)
from gaussian_wrangler import __version__

__author__ = 'hmayes'

# Constants #

# Config File Sections

# Config keys
PDB_TPL_FILE = 'pdb_tpl_file'
GAUSSLOG_FILES_FILE = 'gausslog_list_file'
GAUSSLOG_FILES = 'gausslog_files_list'
GAUSSLOG_FILE = 'gausslog_file'
ONLY_FIRST = 'only_first_coords'
ONLY_FINAL = 'only_final_coords'
OUT_BASE_DIR = 'output_directory'
OUTFILE_NAME = 'output_file_name'
COMBINE_LOGS = 'combine_logs'
ADD_NUM_TO_TYPE = 'add_nums_to_type'

# data file info

# Defaults
DEF_CFG_FILE = 'gausslog2pdb.ini'
DEF_LIST_FILE = 'log_list.txt'
# Set notation
DEF_CFG_VALS = {GAUSSLOG_FILES_FILE: DEF_LIST_FILE,
                OUT_BASE_DIR: None,
                GAUSSLOG_FILE: None,
                PDB_TPL_FILE: None,
                ONLY_FIRST: False,
                ONLY_FINAL: False,
                OUTFILE_NAME: None,
                COMBINE_LOGS: False,
                ADD_NUM_TO_TYPE: False,
                }
REQ_KEYS = {
}


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
        return DEF_CFG_VALS
        # raise IOError('Could not read file {}'.format(f_loc))

    main_proc = cfg_proc(dict(config.items(MAIN_SEC)), DEF_CFG_VALS, REQ_KEYS)

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
    parser.add_argument("-c", "--config", help="The location of the (optional) configuration file in ini format. The "
                                               "default file name is {}, located in the base directory where the " 
                                               "program as run. The program will run using either the specifications "
                                               "from the configuration file or from the command line. Command line "
                                               "specifications will override those in the configuration "
                                               "file.".format(DEF_CFG_FILE),
                        default=DEF_CFG_FILE, type=read_cfg)
    parser.add_argument("-d", "--out_dir", help="The directory where the output files will be placed. This will "
                                                "override any '{}' entry in the configuration file. The default is "
                                                "the same directory as the log file.".format(OUT_BASE_DIR),
                        default=None)
    parser.add_argument("-f", "--file", help="The location of a Gaussian output file. Will override any '{}' entry in "
                                             "the configuration file.".format(GAUSSLOG_FILE), default=None)
    parser.add_argument("-l", "--list", help="The location of the list of Gaussian output files. Will override any "
                                             "'{}' entry in a configuration file.".format(GAUSSLOG_FILES_FILE),
                        default=None)
    parser.add_argument("-o", "--out_fname", help="The name for the created pdb file. If none is provided, it will "
                                                  "take the basename from the provided Gaussian output file name, "
                                                  "with the '.pdb' extension.",
                        default=None)
    parser.add_argument("-t", "--tpl", help="The location of the pdb template file. Will override any '{}'entry in the "
                                            "config file.".format(PDB_TPL_FILE),
                        default=None)

    parser.add_argument("-a", "--only_first", help="Flag to have the program output a pdb only from the first "
                                                   "set of coordinates in the log file. Will override any '{}' entry "
                                                   "in the config file. The default is False.".format(ONLY_FIRST),
                        action="store_true", default=False)

    parser.add_argument("-z", "--only_final", help="Flag to have the program output a pdb only from the last "
                                                   "set of coordinates in the log file. Will override any '{}' entry "
                                                   "in the config file. The default is False.".format(ONLY_FINAL),
                        action="store_true", default=False)

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


def process_gausscom_files(cfg, pdb_tpl_content):
    f_name = ''
    if cfg[COMBINE_LOGS]:
        f_name = create_out_fname(cfg[OUTFILE_NAME], ext='.pdb', base_dir=cfg[OUT_BASE_DIR])
        silent_remove(f_name)
    for gausslog_file in cfg[GAUSSLOG_FILES]:
        if not cfg[PDB_TPL_FILE]:
            pdb_tpl_content[SEC_HEAD] = ["TITLE     {}".format(gausslog_file)]
            pdb_tpl_content[SEC_TAIL] = ["END"]
        if not cfg[COMBINE_LOGS]:
            if cfg[OUTFILE_NAME]:
                out_name = cfg[OUTFILE_NAME]
            else:
                out_name = gausslog_file
            f_name = create_out_fname(out_name, ext='.pdb', base_dir=cfg[OUT_BASE_DIR])
        process_gausslog_file(cfg, gausslog_file, pdb_tpl_content, f_name)


def check_and_print(cfg, atom_id, pdb_tpl_content, gausslog_file, pdb_data_section, f_name, mode, message):
    # Check Num atoms and print
    if cfg[PDB_TPL_FILE]:
        if atom_id != pdb_tpl_content[NUM_ATOMS]:
            raise InvalidDataError('In gausslog file: {}\nfound {} atoms, while the pdb template has {} atoms' 
                                   'atoms'.format(gausslog_file, atom_id, pdb_tpl_content[NUM_ATOMS]))
    list_to_file(pdb_tpl_content[SEC_HEAD] + pdb_data_section + pdb_tpl_content[SEC_TAIL],
                 f_name, list_format=PDB_FORMAT, mode=mode, print_message=message)


def process_gausslog_file(cfg, gausslog_file, pdb_tpl_content, f_name):
    with open(gausslog_file) as d:
        section = SEC_HEAD
        atom_id = 0
        lines_after_coord = 2  # blank line, description, blank line, charge & multiplicity
        if cfg[COMBINE_LOGS]:
            mode = 'a'
        else:
            mode = 'w'
        message = True
        coord_match = False
        first_pass = True

        for line in d:
            line = line.strip()
            if len(line) == 0:
                continue
            # not currently keeping anything from the header
            if section == SEC_HEAD:
                if GAU_COORD_PAT.match(line):
                    coord_match = True
                    if cfg[PDB_TPL_FILE]:
                        pdb_data_section = copy.deepcopy(pdb_tpl_content[SEC_ATOMS])
                    else:
                        pdb_data_section = []
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
                    element_type = ATOM_NUM_DICT[int(split_line[1])]
                except KeyError:
                    raise InvalidDataError("Currently, this code only expects atom numbers up to 36 (Kr), and the "
                                           "atomic number read was {}. Update the code to use this with your current "
                                           "output.".format(split_line[1]))
                # if working from a template, check atom type
                if cfg[PDB_TPL_FILE]:
                    try:
                        pdb_element_type = pdb_data_section[atom_id][8].split(' ')[-1]
                    except IndexError:
                        raise InvalidDataError('Gausslog file: {}\n   has more atoms than the expected {} atoms in '
                                               'the template file: {}'.
                                               format(gausslog_file, pdb_tpl_content[NUM_ATOMS], cfg[PDB_TPL_FILE]))
                    if element_type != pdb_element_type:
                        warning("Atom element types do not match for atom number {}; pdb atom type is {} while "
                                "gausscom type is {}".format(atom_id, pdb_element_type, element_type))
                else:
                    pdb_data_section.append(atom_id)
                    atom_type = element_type
                    if cfg[ADD_NUM_TO_TYPE]:
                        # catch too long atom_type after adding number...
                        max_atom_num_length = 4 - len(atom_type)
                        atom_type += str(atom_id + 1)[:max_atom_num_length]
                    pdb_data_section[atom_id] = ['HETATM', '{:5d}'.format(atom_id + 1), ' {:4} '.format(atom_type),
                                                 'UNL  ', 1, 0.0, 0.0, 0.0,
                                                 '  1.00  0.00          {:>2}'.format(element_type)]
                pdb_data_section[atom_id][5:8] = map(float, split_line[3:6])
                atom_id += 1
            elif section == SEC_TAIL:
                if GAU_E_PAT.match(line):
                    if first_pass:
                        first_pass = False
                    else:
                        del pdb_tpl_content[SEC_HEAD][-1]
                    pdb_tpl_content[SEC_HEAD].append("REMARK    {}".format(line))
                    if not cfg[ONLY_FINAL]:
                        check_and_print(cfg, atom_id, pdb_tpl_content, gausslog_file, pdb_data_section,
                                        f_name, mode, message)
                        if cfg[ONLY_FIRST]:
                            return
                        message = False
                        mode = 'a'
                    section = SEC_HEAD
                    coord_match = False
                    atom_id = 0
                    lines_after_coord = 2

    if len(pdb_tpl_content) == 0:
        raise InvalidDataError("Did not find Gaussian output coordinates in file {}".format(gausslog_file))

    if cfg[ONLY_FINAL]:
        check_and_print(cfg, atom_id, pdb_tpl_content, gausslog_file, pdb_data_section,
                        f_name, mode, message)


def check_input(args, cfg):
    # override config entries if command-line options used
    if args.file:
        cfg[GAUSSLOG_FILE] = args.file
    if args.list:
        cfg[GAUSSLOG_FILES_FILE] = args.list
    if args.tpl:
        cfg[PDB_TPL_FILE] = args.tpl
    if args.out_dir:
        cfg[OUT_BASE_DIR] = args.out_dir
    if args.only_first:
        cfg[ONLY_FIRST] = True
    if args.only_final:
        cfg[ONLY_FINAL] = True
    if args.out_fname:
        cfg[OUTFILE_NAME] = args.out_fname
    if args.out_dir:
        cfg[OUT_BASE_DIR] = args.out_dir

    # checking
    if cfg[COMBINE_LOGS] and not cfg[OUTFILE_NAME]:
        raise InvalidDataError("When combining outputs from multiple log files into one pdb, specify the output "
                               "file name")
    if cfg[COMBINE_LOGS] and not cfg[ONLY_FINAL]:
        warning("When combining outputs from multiple log files into one pdb, only the last coordinates of each "
                "log file will be kept.")
        cfg[ONLY_FINAL] = True

    if cfg[OUT_BASE_DIR]:
        if not os.path.exists(cfg[OUT_BASE_DIR]):
            os.makedirs(cfg[OUT_BASE_DIR])


def main(argv=None):
    print(f"Running GaussianWrangler script gausslog2pdb version {__version__}")

    # Read input
    args, ret = parse_cmdline(argv)
    if ret != GOOD_RET or args is None:
        return ret

    cfg = args.config

    # Read template and data files
    try:
        check_input(args, cfg)

        # set up list of files to process
        cfg[GAUSSLOG_FILES] = []
        if os.path.isfile(cfg[GAUSSLOG_FILES_FILE]):
            with open(cfg[GAUSSLOG_FILES_FILE]) as f:
                for data_file in f:
                    cfg[GAUSSLOG_FILES].append(data_file.strip())
        if cfg[GAUSSLOG_FILE] is not None:
            cfg[GAUSSLOG_FILES].append(cfg[GAUSSLOG_FILE])
        if len(cfg[GAUSSLOG_FILES]) == 0:
            raise InvalidDataError("No files to process: no '{}' specified and "
                                   "no list of files found for: {}".format(GAUSSLOG_FILE, cfg[GAUSSLOG_FILES_FILE]))
        if cfg[ONLY_FIRST] and cfg[ONLY_FINAL]:
            raise InvalidDataError("Cannot specify both '{}' and '{}'".format(ONLY_FIRST, ONLY_FINAL))

        # now start the actual work
        if cfg[PDB_TPL_FILE]:
            pdb_tpl_content = process_pdb_file(cfg[PDB_TPL_FILE])
        else:
            pdb_tpl_content = {}
        process_gausscom_files(cfg, pdb_tpl_content)
    except (IOError, UnicodeDecodeError) as e:
        warning("Problems reading file:", e)
        return IO_ERROR
    except InvalidDataError as e:
        warning("Problems reading data:", e)
        return INVALID_DATA

    return GOOD_RET  # success


if __name__ == '__main__':
    status = main()
    sys.exit(status)
