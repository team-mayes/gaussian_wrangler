#!/usr/bin/env python
"""
Creates pdb data files from lammps data files, given a template pdb file.
"""

import os
import re
import sys
import argparse
from common_wrangler.common import (ATOM_NUM_DICT, NUM_ATOMS, GOOD_RET, INPUT_ERROR, IO_ERROR, INVALID_DATA, BASE_NAME,
                                    SEC_HEAD, SEC_ATOMS, SEC_TAIL, InvalidDataError, warning, check_for_files,
                                    create_out_fname, list_to_file)
from gaussian_wrangler.gw_common import (GAU_COORD_PAT, GAU_SEP_PAT, GAU_E_PAT, GAU_CHARGE_PAT, GAU_STEP_PAT)
from gaussian_wrangler import __version__

__author__ = 'hmayes'


# Constants #

# todo: fix so template for smi2gausscom.py works here too (no extra space)

# Config keys
GAUSSLOG_FILES_FILE = 'gausslog_list_file'
OUT_BASE_DIR = 'output_directory'

# for reading com files
ATOM_TYPES = 'atom_type_dict'
ATOM_HOLDER_PAT = re.compile(r"{atoms}")

# For log file processing
SEC_INITIAL_COORDINATES = 'initial_coordinates_section'


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
    parser.add_argument("-l", "--list", help="The location of the list of Gaussian output files. ",
                        default=None)
    parser.add_argument("-t", "--tpl", help="The location of the Gaussian input template file.")
    parser.add_argument("-e", "--low_energy", help="Flag to take the lowest energy, rather than last, coordinates. "
                                                   "The default is {}.".format(False),
                        action="store_true", default=False)
    parser.add_argument("-c", "--charge_from_tpl", help="Flag to take the charge and multiplicity from the tpl file "
                                                        "rather than from the output file. "
                                                        "The default is {}.".format(False),
                        action="store_true", default=False)
    parser.add_argument("-d", "--out_dir", help="A directory where output files should be saved. The default location "
                                                "is the current directory.", default=None)
    parser.add_argument("-o", "--output_fname", help="The name of the output file to be created. The default is the "
                                                     "output file name with the template base name added to it, and the"
                                                     " '.com' extension.", default=None)
    parser.add_argument("-s", "--step_num", help="The step number in the output file optimization from which to take "
                                                 "the coordinates (instead of last or lowest energy conformation).",
                        default=None)

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


def process_gausslog_files(gausslog_files, com_tpl_content, charge_from_log_flag, find_low_energy, step_num,
                           base_dir, out_fname):
    for gausslog_file in gausslog_files:
        process_gausslog_file(gausslog_file, com_tpl_content, charge_from_log_flag, find_low_energy, step_num,
                              base_dir, out_fname)


def process_gausslog_file(gausslog_file, com_tpl_content, charge_from_log_flag, find_low_energy, step_num,
                          base_dir, out_fname):
    with open(gausslog_file) as d:
        rel_path_fname = os.path.relpath(gausslog_file)
        # The header may be more than 5 lines long--counting from end makes sure the comment goes in the correct line
        if find_low_energy:
            com_tpl_content[SEC_HEAD][-3] = "Low energy conformation from file {}".format(rel_path_fname)
        elif step_num:
            step_num = int(step_num)
            com_tpl_content[SEC_HEAD][-3] = "Conformation from step number {} in file {}".format(step_num,
                                                                                                 rel_path_fname)
        else:
            com_tpl_content[SEC_HEAD][-3] = "Last conformation from file {}".format(rel_path_fname)
        lowest_energy_found = 0.0
        current_step_num = None
        final_atoms_section = []
        atom_type_list = []
        section = SEC_HEAD
        atom_id = 0
        # so don't change the flag that is passed it, so if there is another log file it will also be checked
        if not charge_from_log_flag:
            find_charge = True
        else:
            find_charge = False

        for line in d:
            line = line.strip()
            if len(line) == 0:
                continue
            # not currently keeping anything from the header
            if section == SEC_HEAD:
                if find_charge:
                    if GAU_CHARGE_PAT.match(line):
                        charge_mult = []
                        while find_charge:
                            split_line = line.split('=')
                            charge_mult.append('{}  {}'.format(int(split_line[1].split()[0]),
                                                               int(split_line[2].split()[0])))
                            line = next(d).strip()
                            if not GAU_CHARGE_PAT.match(line):
                                if len(charge_mult) > 1:
                                    section = SEC_INITIAL_COORDINATES
                                    final_atoms_section = []
                                    # already reading the next section, so grab the needed info
                                    atom_type_list.append(line.split()[0])
                                com_tpl_content[SEC_HEAD][-1] = '   '.join(charge_mult)
                                find_charge = False
                        continue
                if step_num and GAU_STEP_PAT.match(line):
                    split_line = line.split()
                    current_step_num = int(split_line[2])
                    if current_step_num == step_num:
                        break

                if GAU_COORD_PAT.match(line):
                    atoms_section = []
                    next(d)
                    next(d)
                    section = SEC_ATOMS
                    continue

            elif section == SEC_INITIAL_COORDINATES:
                while len(line) > 0:
                    # originally just added whole line to final. Then found that this section prints fewer sig figs
                    #   than the coordinate section, so taking those instead
                    atom_type_list.append(line.split()[0])
                    line = next(d).strip()
                while not GAU_COORD_PAT.match(line):
                    line = next(d).strip()
                next(d)
                next(d)
                line = next(d).strip()
                while not GAU_SEP_PAT.match(line):
                    split_line = line.split()
                    atom_xyz = ["{:>12}".format(x) for x in split_line[3:6]]
                    final_atoms_section.append('{:16}'.format(atom_type_list[atom_id]) + ' '.join(atom_xyz))
                    atom_id += 1
                    line = next(d).strip()
                break
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
                    com_atom_type = re.split('[ (]', com_tpl_content[SEC_ATOMS][atom_id])[0].strip()
                    if com_atom_type != atom_type:
                        try:
                            if ATOM_NUM_DICT[int(com_atom_type)] != atom_type:
                                raise ValueError
                        except ValueError:
                            raise InvalidDataError("For atom number {}, {} has atom type '{}', while the template has "
                                                   "atom type '{}'".format(atom_id+1, gausslog_file, atom_type,
                                                                           com_atom_type))
                    atom_type = com_tpl_content[SEC_ATOMS][atom_id]  # This keeps the "fragment" number if there
                atom_type = '{:16}'.format(atom_type)

                atom_xyz = ["{:>12}".format(x) for x in split_line[3:6]]
                atoms_section.append(atom_type + ''.join(atom_xyz))
                atom_id += 1
            elif section == SEC_TAIL:
                if com_tpl_content[NUM_ATOMS] and atom_id != com_tpl_content[NUM_ATOMS]:
                    raise InvalidDataError('In gausslog file: {}\n  found {} atoms, but the tpl expects '
                                           '{} atoms'.format(gausslog_file, atom_id, com_tpl_content[NUM_ATOMS]))
                if GAU_E_PAT.match(line):
                    if find_low_energy:
                        split_line = line.split()
                        energy = float(split_line[4])
                        if energy < lowest_energy_found:
                            final_atoms_section = atoms_section[:]
                    else:
                        final_atoms_section = atoms_section[:]
                    section = SEC_HEAD
                    atom_id = 0

    if len(final_atoms_section) == 0:
        raise InvalidDataError("Check that the following log file has coordinates to use and/or specified step "
                               "number: {}".format(gausslog_file))
    if out_fname:
        f_name = create_out_fname(out_fname, base_dir=base_dir)

    else:
        f_name = create_out_fname(gausslog_file, suffix='_' + com_tpl_content[BASE_NAME], ext='.com', base_dir=base_dir)
    list_to_file(com_tpl_content[SEC_HEAD] + final_atoms_section + com_tpl_content[SEC_TAIL], f_name)

    # Now that finished reading the file, first make sure didn't  exit before reaching the desired number of atoms


def process_gausscom_tpl(com_tpl_file, check_for_charge_mult):
    com_tpl_content = {SEC_HEAD: [], SEC_ATOMS: [], SEC_TAIL: ['', ], ATOM_TYPES: []}
    section = SEC_HEAD
    num_blank_lines_read = 0
    with open(com_tpl_file) as d:
        file_name = os.path.basename(com_tpl_file)
        com_tpl_content[BASE_NAME] = os.path.splitext(file_name)[0]
        for line in d:
            line = line.strip()
            if section == SEC_HEAD:
                com_tpl_content[SEC_HEAD].append(line)
                if len(line) == 0:
                    num_blank_lines_read += 1
                    continue
                if num_blank_lines_read == 2:
                    section = SEC_ATOMS
            elif section == SEC_ATOMS:
                if len(line) == 0 or ATOM_HOLDER_PAT.match(line):
                    section = SEC_TAIL
                    continue
                line_split = line.split()
                # To just get the clean type, may have to remove (fragment=x)
                com_tpl_content[ATOM_TYPES].append(line_split[0].split('(')[0])
                # If there are 5 entries, that means that there is a freeze/no freeze column to keep
                if len(line_split) == 5:
                    com_tpl_content[SEC_ATOMS].append("{:2}{:>8}".format(line_split[0], line_split[1]))
                else:
                    com_tpl_content[SEC_ATOMS].append(line_split[0])
            elif section == SEC_TAIL:
                com_tpl_content[SEC_TAIL].append(line)
    if len(com_tpl_content[SEC_ATOMS]) == 0:
        com_tpl_content[NUM_ATOMS] = None
    else:
        com_tpl_content[NUM_ATOMS] = len(com_tpl_content[SEC_ATOMS])

    # Later code will assume there is a full header, which is five lines (minimum), so make it so now.
    # Blank lines will be filled in, unless the charge is from the template, which is checked next
    while len(com_tpl_content[SEC_HEAD]) < 5:
        com_tpl_content[SEC_HEAD].append('')
    # if getting the charge from the template, make sure they are read.
    if check_for_charge_mult:
        try:
            charge_mult = com_tpl_content[SEC_HEAD][-1].split()
            int(charge_mult[0])
            int(charge_mult[1])
        except (IndexError, ValueError):
            raise InvalidDataError("Option to read charge and multiplicity from template chosen, but not found "
                                   "in template file: {}".format(com_tpl_file))
    # Gaussian needs two blank lines at the end. If they are not already there, add them.
    while len(com_tpl_content[SEC_TAIL]) < 2 or \
            (not com_tpl_content[SEC_TAIL][-2] == '') or (not com_tpl_content[SEC_TAIL][-1] == ''):
        com_tpl_content[SEC_TAIL].append('')

    return com_tpl_content


def main(argv=None):
    print(f"Running GaussianWrangler script gausslog2com version {__version__}")
    args, ret = parse_cmdline(argv)
    if ret != GOOD_RET or args is None:
        return ret

    try:
        # Make sure there are files to process
        gausslog_files = check_for_files(args.file, args.list)

        # and a template file to process
        if not args.tpl:
            raise InvalidDataError("No template file ('-t' option) specified")
        if not os.path.isfile(args.tpl):
            raise IOError(args.tpl)

        # Read template and data files
        com_tpl_content = process_gausscom_tpl(args.tpl, args.charge_from_tpl)
        process_gausslog_files(gausslog_files, com_tpl_content, args.charge_from_tpl, args.low_energy, args.step_num,
                               args.out_dir, args.output_fname)
    except IOError as e:
        warning("Problems reading file:", e)
        return IO_ERROR
    except (InvalidDataError, UnicodeDecodeError) as e:
        warning("Problems reading data:", e)
        return INVALID_DATA

    return GOOD_RET  # success


if __name__ == '__main__':
    status = main()
    sys.exit(status)
