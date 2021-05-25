# coding=utf-8

"""
Common methods for this project.
"""
import re
import os
import collections
import numpy as np
from common_wrangler.common import (InvalidDataError, SEC_HEAD, SEC_ATOMS, SEC_TAIL, BASE_NAME,
                                    ATOM_TYPE, ATOM_COORDS, DIHES, ATOM_NUM_DICT, warning, get_fname_root,
                                    PDB_BEFORE_ELE_LAST_CHAR, PDB_ELE_LAST_CHAR, PDB_MOL_NUM_LAST_CHAR, PDB_Z_LAST_CHAR,
                                    list_to_file)


GAU_HEADER_PAT = re.compile(r"#.*")
GAU_COORD_PAT = re.compile(r"Center {5}Atomic {6}Atomic {13}Coordinates.*")
GAU_SEP_PAT = re.compile(r"---------------------------------------------------------------------.*")
GAU_E_PAT = re.compile(r"SCF Done:.*")
GAU_CHARGE_PAT = re.compile(r"Charge =.*")
GAU_STOICH_PAT = re.compile(r"Stoichiometry.*")
GAU_CONVERG_PAT = re.compile(r"Item {15}Value {5}Threshold {2}Converged?")
GAU_DIH_PAT = re.compile(r"! D.*")
GAU_DERIV_PAT = re.compile(r"! Name {2}Definition {14}Value {10}Derivative Info.*")
GAU_H_PAT = re.compile(r"Sum of electronic and thermal Enthalpies.*")
GAU_STEP_PAT = re.compile(r"Step number.*")
HARM_FREQ_PAT = re.compile(r"Harmonic frequencies.*")
FREQ_PAT = re.compile(r"Frequencies.*")
STOICH = 'Stoichiometry'
CHARGE = 'Charge'
MULT = 'Mult'
MULTIPLICITY = 'Multiplicity'
ENERGY = 'Energy'
CONVERG_STEP_DICT = 'converg_dict'
SCAN_DICT = 'scan_dict'
MAX_FORCE = 'Max Force'
RMS_FORCE = 'RMS Force'
MAX_DISPL = 'Max Displacement'
RMS_DISPL = 'RMS Displacement'
CONVERG = 'Convergence'
CONVERG_ERR = 'Convergence_Error'
ENTHALPY = 'Enthalpy'
GIBBS = 'Gibbs_Free_E'
TS = 'Transition_State'
SCAN_STR = "  Scan  "


def process_gausscom_file(gausscom_file):
    # Grabs and stores in gausscom_content as a dictionary with the keys:
    #    SEC_HEAD: header (route section, blank lines, comments, and full charge and multiplicity line)
    #    CHARGE: overall charge (only) as int
    #    MULT: overall multiplicity (only) as int
    #    SEC_ATOMS: atoms as a dict of dicts, with atom_id as key to dict with
    #        ATOM_TYPE: atom_type (str), ATOM_COORDS: (np array)
    #    SEC_TAIL: everything including and after the blank line following SEC_ATOMS
    with open(gausscom_file) as d:
        gausscom_content = {SEC_HEAD: [], SEC_ATOMS: {}, SEC_TAIL: [],
                            BASE_NAME: get_fname_root(gausscom_file)}
        section = SEC_HEAD
        atom_id = 1
        blank_header_lines = 0

        for line in d:
            line = line.strip()

            if section == SEC_HEAD:
                gausscom_content[SEC_HEAD].append(line)
                if GAU_HEADER_PAT.match(line):
                    continue
                elif len(line) == 0:
                    blank_header_lines += 1
                    if blank_header_lines == 2:
                        section = SEC_ATOMS
                        line = next(d).strip()
                        gausscom_content[SEC_HEAD].append(line)
                        split_line = line.split()
                        try:
                            gausscom_content[CHARGE] = int(split_line[0])
                            gausscom_content[MULT] = int(split_line[1])
                        except (IndexError, ValueError):
                            raise InvalidDataError("Error in reading file {}\n  as a Gaussian input file. On the line "
                                                   "where charge and multiplicity are expected, "
                                                   "found: '{}'".format(gausscom_file, line))

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


def process_gausslog_file(gausslog_file, find_dih=False, find_converg=False, find_step_converg=False,
                          last_step_to_read=None, collect_scan_steps=False):
    # Grabs and stores in gausslog_content as a dictionary with the keys:
    #    (fyi: unlike process_gausscom_file, no SEC_HEAD is collected)
    #    CHARGE: overall charge (only) as int
    #    MULT: overall multiplicity (only) as int
    #    SEC_ATOMS: atoms as a dict of dicts, with atom_id as key to dict with
    #        ATOM_TYPE: atom_type (str), ATOM_COORDS: (np array)
    #    SEC_TAIL: everything including and after the blank line following SEC_ATOMS

    # The mode argument is optional; 'r' will be assumed if it’s omitted. (so read-only below)
    base_name = os.path.basename(gausslog_file)
    with open(gausslog_file) as d:
        gausslog_content = {SEC_ATOMS: {}, BASE_NAME: base_name, STOICH: None, TS: None,
                            ENERGY: np.nan, ENTHALPY: np.nan, GIBBS: np.nan,
                            CONVERG_STEP_DICT: collections.OrderedDict(), SCAN_DICT: {}}
        section = SEC_HEAD
        atom_id = 1
        scan_parameter = None
        # using step convergence for collecting scan step energies, so set flag to true
        if collect_scan_steps:
            find_step_converg = True

        # add stop iteration catch because error can be thrown if EOF reached in one of the while loops
        try:
            for line in d:
                line = line.strip()

                if section == SEC_HEAD:
                    # first find charge & multiplicity; first thing encountered to keep
                    if GAU_CHARGE_PAT.match(line):
                        split_line = line.split('=')
                        gausslog_content[CHARGE] = int(split_line[1].split()[0])
                        gausslog_content[MULT] = int(split_line[2].split()[0])
                        # Dih is the next thing to encounter, which is found in "tail"
                        section = SEC_TAIL

                elif section == SEC_TAIL:
                    # there is not always a dih section in every step
                    while not (GAU_DERIV_PAT.match(line) or GAU_COORD_PAT.match(line) or GAU_STOICH_PAT.match(line)):
                        line = next(d).strip()
                    # sometimes this will execute
                    if (find_dih or collect_scan_steps) and GAU_DERIV_PAT.match(line):
                        if find_dih:
                            while not GAU_DIH_PAT.match(line):
                                line = next(d).strip()
                            dih_dict = {}
                            while GAU_DIH_PAT.match(line):
                                line_split = line.split()
                                dih_dict[line_split[2]] = float(line_split[3])
                                line = next(d).strip()
                            gausslog_content[DIHES] = dih_dict
                        else:
                            # the initial parameters indicate which parameter will be scanned, otherwise search for
                            #     parameter pattern
                            if scan_parameter:
                                search_term = scan_parameter
                            else:
                                search_term = SCAN_STR
                            # move past first separator
                            next(d)
                            line = next(d).strip()
                            while search_term not in line and not GAU_SEP_PAT.match(line):
                                line = next(d).strip()
                            if search_term in line:
                                line_split = line.split()

                                # Only if scan parameter already found will there be an energy collected to store
                                if search_term == SCAN_STR:
                                    scan_parameter = line_split[2]
                                else:
                                    gausslog_content[SCAN_DICT][line_split[3]] = gausslog_content[ENERGY]
                            else:
                                raise InvalidDataError(f"Did not find expected parameter scan info in file: "
                                                       f"{os.path.relpath(gausslog_file)}")
                    # may gave stopped by matching dih pat, so may need to continue
                    #   for some jobs, stoich before coordinates
                    while not GAU_COORD_PAT.match(line) and not GAU_STOICH_PAT.match(line):
                        line = next(d).strip()
                    # catch case when Stoich before coordinates (first step only); catch that case
                    if GAU_STOICH_PAT.match(line):
                        gausslog_content[STOICH] = line.split()[1]
                        while not GAU_COORD_PAT.match(line):
                            line = next(d).strip()
                    next(d)
                    next(d)
                    section = SEC_ATOMS

                elif section == SEC_ATOMS:
                    # will keep overwriting coordinates until it gets to the end
                    while not GAU_SEP_PAT.match(line):
                        split_line = line.split()
                        atom_type = ATOM_NUM_DICT[int(split_line[1])]
                        atom_xyz = np.array(list(map(float, split_line[3:6])))
                        # noinspection PyTypeChecker
                        gausslog_content[SEC_ATOMS][atom_id] = {ATOM_TYPE: atom_type, ATOM_COORDS: atom_xyz}
                        atom_id += 1
                        line = next(d).strip()
                    if not gausslog_content[STOICH]:
                        while not GAU_STOICH_PAT.match(line):
                            line = next(d).strip()
                        gausslog_content[STOICH] = line.split()[1]
                        line = next(d).strip()

                    # Sometimes there is energy & step number before hitting enthalpy, but not in CalcAll jobs
                    while not (GAU_E_PAT.match(line) or GAU_H_PAT.match(line) or HARM_FREQ_PAT.match(line)):
                        line = next(d).strip()
                    if GAU_E_PAT.match(line):
                        gausslog_content[ENERGY] = float(line.split('=')[1].split()[0])
                        line = next(d).strip()

                    # In CalcAll job reading from checkpoint, thermo after Dih then coordinates, then done reading:
                    #     first step: Charge, Coord, Dih, Stoich, SCF, Step, Converg
                    #     then, if not last step: Coord, SCF, Step, Converg
                    #     if last step: Coord, SCF, Step, Converg, Dih, Coord (repeat), Thermo
                    # In Opt then Freq job:
                    #     First and middle steps the same
                    #     Importantly, in Freq: SCF, **Thermo**, then step, conv dih
                    # Thus, Step always after SCF Done, but sometimes thermo in the middle
                    # In CalcAll job starting from coordinates:
                    #     first step: Charge, Dih, Stoich, **Coord**, SCF, Step, Converg
                    while not (HARM_FREQ_PAT.match(line) or GAU_H_PAT.match(line) or GAU_STEP_PAT.match(line)):
                        line = next(d).strip()
                    if HARM_FREQ_PAT.match(line):
                        # checking to hit GAU_H_PAT because one atom jobs won't have a FREQ_PAT
                        while not FREQ_PAT.match(line) and not GAU_H_PAT.match(line):
                            line = next(d).strip()
                        if FREQ_PAT.match(line):
                            split_line = line.split()
                            first_freq = float(split_line[2])
                            if len(split_line) > 3:
                                second_freq = float(split_line[3])
                            else:
                                second_freq = 0.
                            if first_freq < -10.:
                                gausslog_content[TS] = True
                                if first_freq > -40:
                                    warning("Low imaginary frequency for {}: {}".format(base_name, first_freq))
                                if second_freq < -10.:
                                    warning("At least two imaginary frequencies for {}: {:.2f}, {:.2f}"
                                            "".format(base_name, first_freq, second_freq))
                            else:
                                gausslog_content[TS] = False
                            while not GAU_H_PAT.match(line):
                                line = next(d).strip()

                    if GAU_H_PAT.match(line):
                        gausslog_content[ENTHALPY] = float(line.split('=')[1].strip())
                        line = next(d).strip()
                        gausslog_content[GIBBS] = float(line.split('=')[1].strip())
                        line = next(d).strip()

                    # Step num after SCF Done
                    if find_step_converg:
                        while not GAU_STEP_PAT.match(line):
                            line = next(d).strip()
                        split_line = line.split()
                        step_num = int(split_line[2])
                        line = next(d).strip()

                    # convergence after step number
                    if find_converg or find_step_converg:
                        converge_error = False
                        converg = 0.0
                        ind_converg = []
                        while not GAU_CONVERG_PAT.match(line):
                            line = next(d).strip()
                        # in all cases, use the loop below to advance lines, but do not get values if not needed
                        for i in range(4):
                            line = next(d).strip()
                            if find_converg or find_step_converg:
                                line_split = line.split()
                                try:
                                    # Catching case when convergence is so poor that Gaussian prints '********' instead
                                    # of a number (that won't fit) and assign an arbitrarily large convergence penalty
                                    ind_converg.append(float(line_split[2]))
                                    current_converge = ind_converg[i] / float(line_split[3])
                                    if current_converge > 1.0:
                                        converge_error = True
                                    converg += current_converge
                                except ValueError as e:
                                    if '********' in e.args[0]:
                                        ind_converg.append(9.999999)
                                        converg += 2000.00
                                        converge_error = True
                                    else:
                                        raise InvalidDataError(e)
                        # sometimes, Gaussian ignores the convergence error with the line below. Look for it.
                        next(d)
                        line = next(d).strip()
                        if "Optimization completed on the basis of negligible forces." in line:
                            converge_error = "False (negligible forces)"
                        if find_step_converg:
                            # the transformations with step_num_list may not be strictly necessary, but makes IDE happy
                            #    and (likely due to casting the list) removed an error that appeared on a unix cluster
                            #    but not on a local mac
                            step_num_list = list(gausslog_content[CONVERG_STEP_DICT].keys())
                            if step_num in step_num_list:
                                step_num = int(step_num_list[-1]) + 1
                            # noinspection PyTypeChecker
                            gausslog_content[CONVERG_STEP_DICT][step_num] = {ENERGY: gausslog_content[ENERGY],
                                                                             MAX_FORCE: ind_converg[0],
                                                                             RMS_FORCE: ind_converg[1],
                                                                             MAX_DISPL: ind_converg[2],
                                                                             RMS_DISPL: ind_converg[3],
                                                                             CONVERG: converg,
                                                                             CONVERG_ERR: converge_error}
                            # if want to stop after a particular step in a job (so we only care when we are keeping
                            #    track of steps), this is the place to do it,
                            #    as the last thing kept from a middle step is convergence
                            if last_step_to_read and last_step_to_read == step_num:
                                return gausslog_content
                        else:
                            gausslog_content[CONVERG] = converg
                            gausslog_content[CONVERG_ERR] = converge_error

                    section = SEC_TAIL
                    atom_id = 1
        except StopIteration:
            # continue to error checking
            pass

    # Sometimes we want data that is not in the output file. Put placeholders.
    if find_dih and DIHES not in gausslog_content and len(gausslog_content[SEC_ATOMS]) > 3:
        warning("Requested dihedral data not found for file:", os.path.basename(gausslog_file))
        gausslog_content[DIHES] = None
    if not collect_scan_steps:
        if (find_converg or find_step_converg) and CONVERG_ERR not in gausslog_content:
            warning("Did not find final convergence report for file:", os.path.basename(gausslog_file))
            gausslog_content[CONVERG] = np.nan
            gausslog_content[CONVERG_ERR] = None

    return gausslog_content


def get_pdb_coord_list(pdb_str):
    coord_list = []
    pdb_str_list = pdb_str.split("\n")
    for line in pdb_str_list:
        if line.startswith('ATOM') or line.startswith('HETATM'):
            element = line[PDB_BEFORE_ELE_LAST_CHAR:PDB_ELE_LAST_CHAR].strip()
            pdb_xyz = line[PDB_MOL_NUM_LAST_CHAR:PDB_Z_LAST_CHAR]
            coord_list.append(f"{element:6} {pdb_xyz}")
        elif line.startswith('CONECT') or line.startswith('END'):
            break
    return coord_list


def create_com_from_pdb_str(pdb_str, gau_tpl_content, com_fname):
    """
    Extracts one set of pdb coordinates from the "pdb_str" and combines with
    :param pdb_str: str in pdb format
    :param gau_tpl_content: dict with contents of the Gaussian template file
    :param com_fname: str, name of file to be created
    :return:
    """
    coord_list = []
    pdb_str_list = pdb_str.split("\n")
    for line in pdb_str_list:
        if line.startswith('ATOM') or line.startswith('HETATM'):
            element = line[PDB_BEFORE_ELE_LAST_CHAR:PDB_ELE_LAST_CHAR].strip()
            pdb_xyz = line[PDB_MOL_NUM_LAST_CHAR:PDB_Z_LAST_CHAR]
            coord_list.append(["{:6}".format(element), pdb_xyz])
        elif line.startswith('CONECT') or line.startswith('END'):
            break
    list_to_file(gau_tpl_content[SEC_HEAD] + coord_list + gau_tpl_content[SEC_TAIL], com_fname, print_message=False)
