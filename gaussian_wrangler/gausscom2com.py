#!/usr/bin/env python
"""
Creates Gaussian input files from other Gaussian input files.
"""

import argparse
import os
import sys
from configparser import MissingSectionHeaderError
from common_wrangler.common import (InvalidDataError, warning, create_out_fname, list_to_file,
                                    GOOD_RET, INPUT_ERROR, IO_ERROR, INVALID_DATA, check_for_files)
from gaussian_wrangler.gausslog2com import process_gausscom_tpl, ATOM_TYPES
from gaussian_wrangler.gw_common import GAU_HEADER_PAT
from gaussian_wrangler import __version__

__author__ = 'hmayes'


# Constants #

# Config File Sections

# Config keys
COM_TPL_FILE = 'com_tpl_file'
GAUSSCOM_FILES_FILE = 'gausscom_list_file'
GAUSSCOM_FILES = 'gausscom_files_list'
GAUSSCOM_FILE = 'gausscom_file'
OUT_BASE_DIR = 'output_directory'

# data file info

# For pdb template file processing
SEC_HEAD = 'head_section'
SEC_ATOMS = 'atoms_section'
SEC_TAIL = 'tail_section'


def parse_cmdline(argv):
    """
    Returns the parsed argument list and return code.
    `argv` is a list of arguments, or `None` for ``sys.argv[1:]``.
    """
    if argv is None:
        argv = sys.argv[1:]

    # initialize the parser object:
    parser = argparse.ArgumentParser(description="Creates Gaussian input files from Gaussian other files, given a "
                                                 "template input file. The default output name is the same as the "
                                                 "base name of the file with coordinates to use, with the '.com' "
                                                 "extension.")
    parser.add_argument("-c", "--charge_read_com", help="Flag to take the charge and multiplicity from the input file "
                                                        "to be read rather than from the template file. "
                                                        "The default is {}.".format(False),
                        action="store_true", default=False)
    parser.add_argument("-f", "--com_file", help="The location of the Gaussian input file with coordinates to be "
                                                 "input into the template file to create a new input file.",
                        default=None)
    parser.add_argument("-l", "--list_file", help="The location of a text file with a list of Gaussian input files "
                                                  "with coordinates to be input into the template file to create new "
                                                  "input files. Each file name should be on a separate line.",
                        default=None)
    parser.add_argument("-o", "--out_dir", help="A directory where output files should be saved. The default location "
                                                "is the current directory.", default=None)
    parser.add_argument("-t", "--tpl_file", help="The location of the Gaussian input file template (required). "
                                                 "The default file name is {}, located in the "
                                                 "base directory where the program as run.", default=None)

    args = None
    try:
        args = parser.parse_args(argv)
        if args.tpl_file is None:
            raise IOError("A template file (specified with '-t') must be provided.")
        if args.com_file is None and args.list_file is None:
            raise IOError("No files have been specified to be read. Use '-f' to specify a single file, or '-l' to "
                          "specify a file with a list of input files to be read and converted using the specified "
                          "template file.")
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


def process_gausscom_file(gausscom_file, tpl_com_content, read_new_charge, out_dir):
    # to make the later part easier to read
    tpl_atoms = tpl_com_content[SEC_ATOMS]
    tpl_atom_types = tpl_com_content[ATOM_TYPES]
    tpl_atom_num = len(tpl_atom_types)
    with open(gausscom_file) as d:
        section = SEC_HEAD
        atom_id = 0
        atom_content = []

        try:
            for line in d:
                line = line.strip()
                # not currently keeping anything from the header; just check num atoms
                if section == SEC_HEAD:
                    # there may be some instructions (which start with %, and can have some blank lines) before the
                    #    "route card lines" (which start with #)
                    while not GAU_HEADER_PAT.match(line):
                        line = next(d).strip()
                    # skip first line of route card
                    line = next(d).strip()
                    # for "route card" and then description, there may be more than one header line; look for blank line
                    for i in range(2):
                        while len(line) > 0:
                            line = next(d).strip()
                        # now move past the blank line, and get the content of the following line
                        line = next(d).strip()
                    # now on charge, multiplicity line, which we also skip unless we use its charge/mult
                    if read_new_charge:
                        # make sure reading a valid charge/mult line, with at least 2 integers
                        try:
                            charge_mult = line.split()
                            int(charge_mult[0])
                            int(charge_mult[1])
                            if len(charge_mult) % 2 != 0:
                                raise IndexError
                        except (IndexError, ValueError):
                            raise InvalidDataError("Problem while reading file: {}\nOption to read charge and "
                                                   "multiplicity from template not chosen, but found invalid data on "
                                                   "the expected line: {}".format(os.path.basename(gausscom_file),
                                                                                  line))
                        tpl_com_content[SEC_HEAD][-1] = line
                    section = SEC_ATOMS
                    continue

                elif section == SEC_ATOMS:
                    # stay in atom section until a blank line is reached
                    while len(line) > 0:
                        split_line = line.split()
                        # if there is a freeze/no freeze col, will be 5 columns (split by ' '); Keep atom info together
                        if len(split_line) == 5:
                            atom_info = "{:2}{:>8}".format(split_line[0], split_line[1])
                        else:
                            atom_info = split_line[0]

                        # if template has atoms, check atom type
                        if tpl_atom_num > 0:
                            atom_type = atom_info.split()[0].split('(')[0]
                            if atom_type != tpl_atom_types[atom_id]:
                                raise InvalidDataError("Problem while reading file: {}\nAtom types do not match for "
                                                       "atom number {}: file has type {} while tpl has type "
                                                       "{}".format(os.path.basename(gausscom_file), atom_id + 1,
                                                                   tpl_atom_types[atom_id], atom_type))
                            atom_info = tpl_atoms[atom_id]

                        atom_xyz = ["{:>12}".format(x) for x in split_line[-3:]]
                        atom_content.append('{:18}'.format(atom_info) + '  '.join(atom_xyz))
                        atom_id += 1
                        line = next(d).strip()
                    # Don't need to read the tail, because we won't use it
                    break
        except StopIteration:
            pass
        except UnicodeDecodeError:
            raise InvalidDataError(f"Error in reading file: {gausscom_file}\n           Exiting program.")

        # now loop is done; check atom number if atoms are in the tpl file
        check_num_atoms(atom_id, gausscom_file, tpl_atom_num)

        f_name = create_out_fname(gausscom_file, ext='.com', base_dir=out_dir)
        list_to_file(tpl_com_content[SEC_HEAD] + atom_content + tpl_com_content[SEC_TAIL], f_name)


def check_num_atoms(atom_id, gausscom_file, tpl_atom_num):
    if tpl_atom_num > 0:
        if tpl_atom_num != atom_id:
            raise InvalidDataError("Problem while reading file: {}\nFound {} atoms in this file,  "
                                   "while the tpl file has {} "
                                   "atoms".format(os.path.basename(gausscom_file), atom_id, tpl_atom_num))


def main(argv=None):
    print(f"Running GaussianWrangler script gausscom2com version {__version__}")
    args, ret = parse_cmdline(argv)
    if ret != GOOD_RET or args is None:
        return ret

    # Read template and data files
    try:
        # Make sure there are files to process
        gausscom_files = check_for_files(args.com_file, args.list_file)
        com_tpl_content = process_gausscom_tpl(args.tpl_file, not args.charge_read_com)
        for gausscom_file in gausscom_files:
            process_gausscom_file(gausscom_file, com_tpl_content, args.charge_read_com, args.out_dir)

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
