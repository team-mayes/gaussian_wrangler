#!/usr/bin/env python
"""
Homolytic fragmenter
"""

from __future__ import print_function
import re
import os
import sys
import argparse
import numpy as np
from nrel_tools.common import (InvalidDataError, warning, process_cfg, create_out_fname, list_to_file,
                               GOOD_RET, INPUT_ERROR, IO_ERROR, INVALID_DATA)

try:
    # noinspection PyCompatibility
    from ConfigParser import ConfigParser, MissingSectionHeaderError
except ImportError:
    # noinspection PyCompatibility
    from configparser import ConfigParser, MissingSectionHeaderError

__author__ = 'hmayes'


# Constants #

# Config File Sections
MAIN_SEC = 'main'

# Config keys
GAUSSCOM_FILE = 'gausscom_file'
OUT_BASE_DIR = 'output_directory'
CUT_ATOMS = 'cut_atoms'
GAUSS_COMMAND = 'gaussian_options_line'
GAUSS_CP_COMMAND = 'gaussian_cp_options_line'

# data file info

# Defaults
DEF_CFG_FILE = 'gausscom_fragment.ini'
DEF_GAUSS_COMMAND = '# m062x/Def2SVP nosymm scf=xqc opt freq'
DEF_GAUSS_CP_COMMAND = '# m062x/Def2TZVP nosymm Counterpoise=2'

# Set notation
DEF_CFG_VALS = {OUT_BASE_DIR: None,
                GAUSSCOM_FILE: None,
                GAUSS_COMMAND: DEF_GAUSS_COMMAND,
                GAUSS_CP_COMMAND: DEF_GAUSS_CP_COMMAND,
                }
REQ_KEYS = {CUT_ATOMS: str,
            }

# For file processing
CUT_PAIR_LIST = 'cut_pair_list'
GAU_HEADER_PAT = re.compile(r"#.*")
SEC_HEAD = 'head_section'
SEC_ATOMS = 'atoms_section'
SEC_TAIL = 'tail_section'
ATOM_TYPE = 'atom_type'
ATOM_COORDS = 'atom_coords'
FRAGMENT = 'fragment'
MAX_BOND_DIST = 2.0  # same length units as in input and output file, here Angstroms


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

    cut_pairs = main_proc[CUT_ATOMS].split(';')
    main_proc[CUT_PAIR_LIST] = []
    for pair in cut_pairs:
        if pair == '':
            continue
        atom_pair = [int(x) for x in pair.split(',')]
        if len(atom_pair) != 2:
            raise InvalidDataError("The '{}' values should be sets of two atoms, separated by commas, with each pair "
                                   "separated by ';'".format(CUT_PAIR_LIST))
        main_proc[CUT_PAIR_LIST].append(atom_pair)

    if not os.path.exists(main_proc[OUT_BASE_DIR]):
        os.makedirs(main_proc[OUT_BASE_DIR])

    return main_proc


def parse_cmdline(argv):
    """
    Returns the parsed argument list and return code.
    `argv` is a list of arguments, or `None` for ``sys.argv[1:]``.
    """
    if argv is None:
        argv = sys.argv[1:]

    # initialize the parser object:
    parser = argparse.ArgumentParser(description='Creates cp files from Gaussian input files, given a list of atom '
                                                 'numbers where to cut (list format: atom1, atom2; atom3, atom4.')
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


def process_gausscom_file(gausscom_file):
    with open(gausscom_file) as d:
        gausscom_content = {SEC_HEAD: [], SEC_ATOMS: {}, SEC_TAIL: []}
        section = SEC_HEAD
        atom_id = 1
        lines_after_header = 4  # blank line, description, blank line, charge & multiplicity

        for line in d:
            line = line.strip()

            if section == SEC_HEAD:
                gausscom_content[SEC_HEAD].append(line)
                if GAU_HEADER_PAT.match(line):
                    continue
                elif lines_after_header > 0:
                    lines_after_header -= 1
                    if lines_after_header == 0:
                        section = SEC_ATOMS
                    continue

            elif section == SEC_ATOMS:
                if len(line) == 0:
                    section = SEC_TAIL
                    gausscom_content[SEC_TAIL].append(line)
                    continue
                split_line = line.split()

                atom_type = split_line[0]
                atom_xyz = np.array(list(map(float, split_line[1:4])))
                gausscom_content[SEC_ATOMS][atom_id] = {ATOM_TYPE: atom_type, ATOM_COORDS: atom_xyz}
                atom_id += 1
            elif section == SEC_TAIL:
                gausscom_content[SEC_TAIL].append(line)

    return gausscom_content


def calc_dist(a, b):
    return np.linalg.norm(np.subtract(a, b))


def validate_atom_num(atom_pair, atoms_content, gausscom_file):
    # check that both atom numbers are not larger than the total number of atoms,
    # and that they are close enough to be bonded
    for atom_num in atom_pair:
        if atom_num not in atoms_content:
            raise InvalidDataError("Found atom id {} in '{}', but there are only {} atoms in the file {}"
                                   "".format(atom_num, CUT_ATOMS, len(atoms_content), gausscom_file))
    pair_dist = calc_dist(atoms_content[atom_pair[0]][ATOM_COORDS], atoms_content[atom_pair[1]][ATOM_COORDS])
    if pair_dist > MAX_BOND_DIST:
        raise InvalidDataError("Atom ids {} and {} are {:.2f} Angstroms apart, which is greater than the tested "
                               "maximum bond distance of {:.2f}".format(atom_pair[0], atom_pair[1], MAX_BOND_DIST))


def fragment_molecule(atom_pair, atoms_content):
    single_bond_atoms = ['H', 'Cl', ]
    atom_numbers = list(range(1,len(atoms_content)+1))
    frag1_list = []
    frag2_list = []
    for atom in atom_pair:
        # Check if fragment made up of just one atom
        lonely_frag = False
        if atoms_content[atom][ATOM_TYPE] == 'O':
            for other_atom in atom_numbers:
                if other_atom == atom:
                    continue
                pair_dist = calc_dist(atoms_content[atom][ATOM_COORDS], atoms_content[other_atom][ATOM_COORDS])
                if pair_dist < MAX_BOND_DIST:
                    lonely_frag = True
                    break
        elif atoms_content[atom][ATOM_TYPE] in single_bond_atoms:
            lonely_frag = True
        if lonely_frag:
            atoms_content[atom][FRAGMENT] = 1
            frag1_list.append(atom)
            atom_numbers.remove(atom)
            frag2_list = atom_numbers
            for other_atom in atom_numbers:
                atoms_content[other_atom][FRAGMENT] = 2
            return frag1_list, frag2_list
    return frag1_list, frag2_list


def write_com_file(cp_file_name, gauss_command, for_comment_line, atoms_content, frag_num = None, frag_list = []):
    # Don't bother making a separate file if just one atom; there would be lots of repeat calculations that way
    if len(frag_list) == 1:
        return
    if frag_num:
        comment_begin = 'radical calculation of fragment {} '.format(frag_num)
        charge_mult = '0 2'
    else:
        comment_begin = 'cp calculation '
        charge_mult = '0 1   0 2    0 2'
        frag_list = range(1,len(atoms_content)+1)
    print_list = [[gauss_command], [], [comment_begin + for_comment_line], [], [charge_mult]]
    for atom_num in frag_list:
        if frag_num:
            atom_line = atoms_content[atom_num][ATOM_TYPE] + \
                        '    {:11.6f} {:11.6f} {:11.6f}'.format(*atoms_content[atom_num][ATOM_COORDS])
        else:
            atom_line = atoms_content[atom_num][ATOM_TYPE] + \
                        '(Fragment={})    {:11.6f} {:11.6f} {:11.6f}'.format(atoms_content[atom_num][FRAGMENT],
                                                                             *atoms_content[atom_num][ATOM_COORDS])
        print_list.append(atom_line)
    print_list.append([])
    print_list.append([])
    list_to_file(print_list, cp_file_name)


def print_com_files(atom_pair, atoms_content, gausscom_file, cfg, frag1, frag2):
    for_comment_line = 'from fragment pair {} and {}'.format(atom_pair, gausscom_file)
    # First print template for CP calc (the coordinates should later be replaced by further optimized coordinates,
    # if desired)
    cp_file_name = create_out_fname(gausscom_file, suffix='_{}_{}_cp'.format(*atom_pair),
                                    ext='.com', base_dir=cfg[OUT_BASE_DIR])
    write_com_file(cp_file_name, cfg[GAUSS_CP_COMMAND], for_comment_line, atoms_content)
    frag1_file_name = create_out_fname(gausscom_file, suffix='_{}_{}_f1'.format(*atom_pair), base_dir=cfg[OUT_BASE_DIR])
    write_com_file(frag1_file_name, cfg[GAUSS_COMMAND], for_comment_line, atoms_content, 1, frag1)
    frag2_file_name = create_out_fname(gausscom_file, suffix='_{}_{}_f2'.format(*atom_pair), base_dir=cfg[OUT_BASE_DIR])
    write_com_file(frag2_file_name, cfg[GAUSS_COMMAND], for_comment_line, atoms_content, 2, frag2)


def main(argv=None):
    # Read input
    args, ret = parse_cmdline(argv)
    if ret != GOOD_RET or args is None:
        return ret

    cfg = args.config
    gausscom_file = cfg[GAUSSCOM_FILE]

    # Read template and data files
    try:
        gausscom_content = process_gausscom_file(gausscom_file)
        # Before making files, check that atom numbers are valid
        for atom_pair in cfg[CUT_PAIR_LIST]:
            validate_atom_num(atom_pair, gausscom_content[SEC_ATOMS], gausscom_file)
        for atom_pair in cfg[CUT_PAIR_LIST]:
            frag1, frag2 = fragment_molecule(atom_pair, gausscom_content[SEC_ATOMS])
            print_com_files(atom_pair, gausscom_content[SEC_ATOMS], gausscom_file, cfg, frag1, frag2)
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
