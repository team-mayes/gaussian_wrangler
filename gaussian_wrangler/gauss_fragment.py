#!/usr/bin/env python
"""
homolytic fragmenter
"""

from __future__ import print_function
import os
import sys
import argparse
from configparser import ConfigParser, MissingSectionHeaderError
from common_wrangler.common import (InvalidDataError, warning, process_cfg, create_out_fname, list_to_file,
                                    GOOD_RET, INPUT_ERROR, IO_ERROR, INVALID_DATA,
                                    ATOM_TYPE, ATOM_COORDS, MAIN_SEC, calc_dist, SEC_ATOMS)
from gaussian_wrangler.gw_common import process_gausscom_file, process_gausslog_file, CHARGE, MULT
from gaussian_wrangler import __version__


__author__ = 'hmayes'


# Constants #

# Config File Sections

# Config keys
GAUSSCOM_FILE = 'input_com_file'
GAUSSLOG_FILE = 'input_log_file'
OUT_BASE_DIR = 'output_directory'
CUT_ATOMS = 'cut_atoms'
GAUSS_COMMAND = 'gaussian_options_line'
GAUSS_END = 'gaussian_options_end'
GAUSS_CP_COMMAND = 'gaussian_cp_options_line'
GAUSS_CP_END = 'gaussian_cp_options_end'
TWO_MOLECULES = 'two_molecules'


# data file info

# Defaults
DEF_CFG_FILE = 'gausscom_fragment.ini'
DEF_GAUSS_COMMAND = '# m062x/Def2SVP nosymm scf=xqc opt guess=mix freq=noraman CPHF=Grid=Fine'
DEF_GAUSS_CP_COMMAND = '# m062x/Def2TZVP nosymm Counterpoise=2 CPHF=Grid=Fine'

# Set notation
DEF_CFG_VALS = {OUT_BASE_DIR: None,
                GAUSSCOM_FILE: None,
                GAUSSLOG_FILE: None,
                GAUSS_COMMAND: DEF_GAUSS_COMMAND,
                GAUSS_END: None,
                GAUSS_CP_COMMAND: DEF_GAUSS_CP_COMMAND,
                GAUSS_CP_END: None,
                TWO_MOLECULES: False,
                }
REQ_KEYS = {CUT_ATOMS: str,
            }

# For file processing
CUT_PAIR_LIST = 'cut_pair_list'
FRAGMENT = 'fragment'
MAX_BOND_DIST = 1.9  # same length units as in input and output file, here Angstroms
MAX_H_BOND_DIST = 1.5  # same length units as in input and output file, here Angstroms
MAX_M_BOND_DIST = 2.3  # same length units as in input and output file, here Angstroms
METALS = ['Ti', 'Sb', 'Ge']


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

    if main_proc[GAUSSCOM_FILE] and main_proc[GAUSSLOG_FILE]:
        raise InvalidDataError("Both an '{}' and an '{}' are specified. Choose only one file to read and fragment.".
                               format(GAUSSCOM_FILE, GAUSSLOG_FILE))

    if main_proc[OUT_BASE_DIR]:
        if not os.path.exists(main_proc[OUT_BASE_DIR]):
            os.makedirs(main_proc[OUT_BASE_DIR])

    for convert_to_line_breaks_key in [GAUSS_END, GAUSS_CP_END]:
        if main_proc[convert_to_line_breaks_key]:
            main_proc[convert_to_line_breaks_key] = main_proc[convert_to_line_breaks_key].replace(';', '\n')

    return main_proc


# noinspection DuplicatedCode
def parse_cmdline(argv):
    """
    Returns the parsed argument list and return code.
    `argv` is a list of arguments, or `None` for ``sys.argv[1:]``.
    """
    if argv is None:
        argv = sys.argv[1:]

    # initialize the parser object:
    parser = argparse.ArgumentParser(description='Creates cp files from Gaussian input files, given a list of atom '
                                                 'numbers where to cut (list format: atom1, atom2; atom3, atom4).')
    parser.add_argument("-c", "--config", help=f"The location of the configuration file in ini format. "
                                               f"The default file name is {DEF_CFG_FILE}, located in the "
                                               f"base directory where the program as run.",
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


def validate_atom_num(atom_pair, atoms_content, gauss_in_fname, ignore_max_dist):
    # check that both atom numbers are not larger than the total number of atoms,
    # and that they are close enough to be bonded
    for atom_num in atom_pair:
        if atom_num not in atoms_content:
            raise InvalidDataError("Found atom id {} in '{}', but there are only {} atoms in the file {}"
                                   "".format(atom_num, CUT_ATOMS, len(atoms_content), gauss_in_fname))
    pair_dist = calc_dist(atoms_content[atom_pair[0]][ATOM_COORDS], atoms_content[atom_pair[1]][ATOM_COORDS])
    if pair_dist > MAX_BOND_DIST and not ignore_max_dist:
        raise InvalidDataError("Atom ids {} and {} are {:.2f} Angstroms apart, which is greater than "
                               "maximum bond distance of {:.2f}".format(atom_pair[0], atom_pair[1], pair_dist,
                                                                        MAX_BOND_DIST))


def fragment_molecule(atom_pair, atoms_content, ignore_max_dist):
    broke_triple_bond = False
    broke_double_bond = False
    broke_double_bond_list = [False, False]
    single_bond_atoms = ['H', 'Cl', ]
    unassigned_atom_numbers = list(range(1, len(atoms_content)+1))
    frag1_list = []
    frag2_list = []
    if not ignore_max_dist:
        for atom in atom_pair:
            # Check if fragment made up of just one atom
            lonely_frag = False
            test_atom_type = atoms_content[atom][ATOM_TYPE]
            if test_atom_type in single_bond_atoms:
                lonely_frag = True
            elif test_atom_type == 'O' or test_atom_type == 'N':
                for other_atom in unassigned_atom_numbers:
                    lonely_frag = True
                    if test_atom_type == 'O':
                        broke_double_bond = True
                    else:
                        broke_triple_bond = True
                    if other_atom == atom_pair[0] or other_atom == atom_pair[1]:
                        continue
                    pair_dist = calc_dist(atoms_content[atom][ATOM_COORDS], atoms_content[other_atom][ATOM_COORDS])
                    if atoms_content[other_atom][ATOM_TYPE] == 'H':
                        max_dist = MAX_H_BOND_DIST
                    elif atoms_content[other_atom][ATOM_TYPE] in METALS:
                        max_dist = MAX_M_BOND_DIST
                    else:
                        max_dist = MAX_BOND_DIST
                    if pair_dist < max_dist:
                        lonely_frag = False
                        if test_atom_type == 'O':
                            broke_double_bond = False
                        else:
                            broke_triple_bond = False
                        break
            if lonely_frag:
                atoms_content[atom][FRAGMENT] = 1
                frag1_list.append(atom)
                unassigned_atom_numbers.remove(atom)
                frag2_list = unassigned_atom_numbers
                for other_atom in unassigned_atom_numbers:
                    atoms_content[other_atom][FRAGMENT] = 2
                return frag1_list, frag2_list, broke_double_bond, broke_triple_bond

    # Now the more difficult cases
    frag1_list.append(atom_pair[0])
    unassigned_atom_numbers.remove(atom_pair[0])
    atoms_content[atom_pair[0]][FRAGMENT] = 1
    frag2_list.append(atom_pair[1])
    unassigned_atom_numbers.remove(atom_pair[1])
    atoms_content[atom_pair[1]][FRAGMENT] = 2

    if not ignore_max_dist:
        # first check for C=C
        for counter, atom in enumerate(atom_pair):
            if atoms_content[atom][ATOM_TYPE] == 'C':
                bonded_to_c = []
                for other_atom in unassigned_atom_numbers:
                    pair_dist = calc_dist(atoms_content[atom][ATOM_COORDS], atoms_content[other_atom][ATOM_COORDS])
                    if atoms_content[other_atom][ATOM_TYPE] == 'H':
                        max_dist = MAX_H_BOND_DIST
                    elif atoms_content[other_atom][ATOM_TYPE] in METALS:
                        max_dist = MAX_M_BOND_DIST
                    else:
                        max_dist = MAX_BOND_DIST
                    if pair_dist < max_dist:
                        bonded_to_c.append(other_atom)
                if len(bonded_to_c) == 2:
                    type1 = atoms_content[bonded_to_c[0]][ATOM_TYPE]
                    type2 = atoms_content[bonded_to_c[0]][ATOM_TYPE]
                    if type1 == 'O' and type2 == 'O':
                        broke_double_bond_list[counter] = False
                    else:
                        broke_double_bond_list[counter] = True
        # if one atom has a double-bond but not the other, then a it is a single bond that is broken, so leave
        # broke_double_bond as false, otherwise:
        if broke_double_bond_list[0] == broke_double_bond_list[1]:
            broke_double_bond = broke_double_bond_list[0]

    # first add to frag 1
    atoms_to_check = [atom_pair[0]]
    add_atoms_to_fragment(unassigned_atom_numbers, atoms_content, atoms_to_check, frag1_list, 1, single_bond_atoms)
    # make sure no atoms in fragment 1 are within bonding distance of any atoms remaining in the
    # unassigned_atom_numbers list
    for f1_atom in frag1_list:
        for atom in unassigned_atom_numbers:
            pair_dist = calc_dist(atoms_content[atom][ATOM_COORDS], atoms_content[f1_atom][ATOM_COORDS])
            if atoms_content[atom][ATOM_TYPE] == 'H' or atoms_content[f1_atom][ATOM_TYPE] == 'H':
                max_dist = MAX_H_BOND_DIST
            elif atoms_content[atom][ATOM_TYPE] in METALS or atoms_content[f1_atom][ATOM_TYPE] in METALS:
                max_dist = MAX_M_BOND_DIST
            else:
                max_dist = MAX_BOND_DIST
            if pair_dist < max_dist:
                raise InvalidDataError("Found that atom {} assigned to fragment 1 is within {} Angstroms of atom {} "
                                       "which was not assigned to fragment 1".format(f1_atom, max_dist, atom))
    # check that all remaining atoms are bonded to each other
    atoms_to_check = [atom_pair[1]]
    add_atoms_to_fragment(unassigned_atom_numbers, atoms_content, atoms_to_check, frag2_list, 2, single_bond_atoms)
    if len(unassigned_atom_numbers) > 0:
        raise InvalidDataError("Atoms {} were not assigned to either fragment 1 or 2.".format(unassigned_atom_numbers))
    frag1_list.sort()
    frag2_list.sort()
    if len(frag1_list) > len(frag2_list):
        return frag2_list, frag1_list, broke_double_bond, broke_triple_bond
    else:
        return frag1_list, frag2_list, broke_double_bond, broke_triple_bond


def add_atoms_to_fragment(atom_numbers, atoms_content, atoms_to_check, frag_list, frag_num, single_bond_atoms):
    add_to_atoms_to_check = []
    while len(atoms_to_check) > 0:
        for check_atom in atoms_to_check:
            atoms_to_remove_from_atom_list = []
            for atom in atom_numbers:
                pair_dist = calc_dist(atoms_content[atom][ATOM_COORDS], atoms_content[check_atom][ATOM_COORDS])
                if atoms_content[atom][ATOM_TYPE] == 'H' or atoms_content[check_atom][ATOM_TYPE] == 'H':
                    max_dist = MAX_H_BOND_DIST
                elif atoms_content[atom][ATOM_TYPE] in METALS or atoms_content[check_atom][ATOM_TYPE] in METALS:
                    max_dist = MAX_M_BOND_DIST
                else:
                    max_dist = MAX_BOND_DIST
                if pair_dist < max_dist:
                    frag_list.append(atom)
                    atoms_content[atom][FRAGMENT] = frag_num
                    # avoid changing list while iterating
                    atoms_to_remove_from_atom_list.append(atom)
                    if atoms_content[atom][ATOM_TYPE] not in single_bond_atoms:
                        add_to_atoms_to_check.append(atom)
            for atom in atoms_to_remove_from_atom_list:
                atom_numbers.remove(atom)
        atoms_to_check = []
        for atom in add_to_atoms_to_check:
            atoms_to_check.append(atom)
        add_to_atoms_to_check = []


def write_com_file(com_file_name, gauss_command, gauss_end, for_comment_line, atoms_content, broke_double_bond,
                   broke_triple_bond, ignore_max_dist, charge, mult, frag_num=None, frag_list=None):
    """
    After figuring out the fragments, make Gaussian input files to calculate the counterpoint correction (if a non-zero
    list is passed to "frag_list". Otherwise, make a Gaussian input file to optimize any fragments with len > 1.
    :param com_file_name: str
    :param gauss_command: str
    :param gauss_end: str or None
    :param for_comment_line: str
    :param atoms_content: dictionary with atom type, atom coordinates, and fragment ID
    :param broke_double_bond: flag to change multiplicity
    :param broke_triple_bond: flag to change multiplicity
    :param ignore_max_dist: flag if not making radicals
    :param charge: int read from com file
    :param mult: int read from com file
    :param frag_num: optional integer that will be used to name the file
    :param frag_list: optional list that will be used to make a Gaussian input file with only the atoms in that fragment
    :return: nothing
    """
    main_charge = str(charge) + ' ' + str(mult)
    if broke_double_bond:
        frag_charge = str(charge) + ' ' + str(mult + 2)
    elif broke_triple_bond:
        frag_charge = str(charge) + ' ' + str(mult + 3)
    elif ignore_max_dist:
        frag_charge = str(charge) + ' ' + str(mult)
    else:
        frag_charge = str(charge) + ' ' + str(mult + 1)
    if frag_list is None:
        frag_list = []
    if frag_num:
        # Don't bother making a separate file if just a few atoms; it would re
        if len(frag_list) < 7:
            element_list = []
            for atom_num in frag_list:
                element_list.append(atoms_content[atom_num][ATOM_TYPE])
            print("Fragment {} is only {}".format(os.path.basename(com_file_name), element_list))
            return
        charge_mult = frag_charge
        if ignore_max_dist:
            comment_begin = 'optimization of fragment {} '.format(frag_num)
        else:
            comment_begin = 'radical calculation of fragment {} '.format(frag_num)
    else:
        comment_begin = 'cp calculation '
        charge_mult = main_charge + '   ' + frag_charge + '   ' + frag_charge
        frag_list = range(1, len(atoms_content)+1)
    print_list = [[gauss_command], [], [comment_begin + for_comment_line], [], [charge_mult]]
    for atom_num in frag_list:
        if frag_num:
            atom_line = "{:7}".format(atoms_content[atom_num][ATOM_TYPE]) + \
                        ' {:11.6f} {:11.6f} {:11.6f}'.format(*atoms_content[atom_num][ATOM_COORDS])
        else:
            atom_type_str = "{}(Fragment={})".format(atoms_content[atom_num][ATOM_TYPE],
                                                     atoms_content[atom_num][FRAGMENT])
            atom_line = '{:15}  {:11.6f} {:11.6f} {:11.6f}'.format(atom_type_str, *atoms_content[atom_num][ATOM_COORDS])
        print_list.append(atom_line)
    print_list.append([])
    if gauss_end:
        print_list.append([gauss_end])
        print_list.append([])
    print_list.append([])
    list_to_file(print_list, com_file_name)


def print_com_files(atom_pair, atoms_content, gauss_in_fname, cfg, frag1, frag2, broke_double_bond, broke_triple_bond,
                    ignore_max_dist, charge, mult):
    for_comment_line = 'from fragment pair {} and {}'.format(atom_pair, gauss_in_fname)
    # First print template for CP calc (the coordinates should later be replaced by further optimized coordinates,
    # if desired)
    cp_file_name = create_out_fname(gauss_in_fname, suffix='_{}_{}_cp'.format(*atom_pair),
                                    ext='.com', base_dir=cfg[OUT_BASE_DIR])
    write_com_file(cp_file_name, cfg[GAUSS_CP_COMMAND], cfg[GAUSS_CP_END], for_comment_line, atoms_content,
                   broke_double_bond, broke_triple_bond, ignore_max_dist, charge, mult)
    frag1_file_name = create_out_fname(gauss_in_fname, suffix='_{}_{}_f1'.format(*atom_pair), ext='.com',
                                       base_dir=cfg[OUT_BASE_DIR])
    write_com_file(frag1_file_name, cfg[GAUSS_COMMAND], cfg[GAUSS_END], for_comment_line, atoms_content,
                   broke_double_bond, broke_triple_bond, ignore_max_dist, charge, mult, 1, frag1)
    frag2_file_name = create_out_fname(gauss_in_fname, suffix='_{}_{}_f2'.format(*atom_pair), ext='.com',
                                       base_dir=cfg[OUT_BASE_DIR])
    write_com_file(frag2_file_name, cfg[GAUSS_COMMAND], cfg[GAUSS_END], for_comment_line, atoms_content,
                   broke_double_bond, broke_triple_bond, ignore_max_dist, charge, mult, 2, frag2)


def main(argv=None):
    print(f"Running GaussianWrangler script gauss_fragment version {__version__}")

    # Read input
    args, ret = parse_cmdline(argv)
    if ret != GOOD_RET or args is None:
        return ret

    cfg = args.config

    # Read template and data files
    try:
        if cfg[GAUSSCOM_FILE]:
            gauss_file = cfg[GAUSSCOM_FILE]
            gauss_in_content = process_gausscom_file(gauss_file)
            atom_data = gauss_in_content[SEC_ATOMS]
        elif cfg[GAUSSLOG_FILE]:
            gauss_file = cfg[GAUSSLOG_FILE]
            gauss_in_content = process_gausslog_file(gauss_file)
            atom_data = gauss_in_content[SEC_ATOMS]
        else:
            raise InvalidDataError("This program requires either a valid Gaussian input file ('{}') or Gaussian "
                                   "output file ('{}') from which to extract atoms with their coordinates.")
        # Before making files, check that atom numbers are valid
        for atom_pair in cfg[CUT_PAIR_LIST]:
            validate_atom_num(atom_pair, atom_data, gauss_file, cfg[TWO_MOLECULES])
        for atom_pair in cfg[CUT_PAIR_LIST]:
            frag1, frag2, broke_double_bond, broke_triple_bond = fragment_molecule(atom_pair, atom_data,
                                                                                   cfg[TWO_MOLECULES])
            print_com_files(atom_pair, atom_data, gauss_file, cfg, frag1, frag2, broke_double_bond, broke_triple_bond,
                            cfg[TWO_MOLECULES], gauss_in_content[CHARGE], gauss_in_content[MULT])
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
