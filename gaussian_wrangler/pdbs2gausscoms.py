#!/usr/bin/env python
"""
Creates Gaussian files from pdb files.

Provide template input file
List of pdb files to convert; there may be multiple structures in one file; each one must end with "END/n"
"""

import argparse
import os
import sys
import numpy as np
from configparser import ConfigParser, MissingSectionHeaderError
from operator import itemgetter
from rdkit import RDLogger
from rdkit.Chem.rdmolfiles import MolFromPDBFile, MolToPDBBlock
from rdkit.Chem.rdMolTransforms import GetDihedralDeg, SetDihedralDeg
from rdkit.Chem.AllChem import MMFFOptimizeMoleculeConfs
from common_wrangler.common import (GOOD_RET, INPUT_ERROR, IO_ERROR, INVALID_DATA,
                                    PDB_LINE_TYPE_LAST_CHAR, PDB_MOL_NUM_LAST_CHAR, PDB_Z_LAST_CHAR,
                                    PDB_BEFORE_ELE_LAST_CHAR, PDB_ELE_LAST_CHAR, PDB_ATOM_NUM_LAST_CHAR,
                                    PDB_ATOM_TYPE_LAST_CHAR, MAIN_SEC, SEC_HEAD, SEC_TAIL,
                                    InvalidDataError, warning, create_out_fname, list_to_file, process_cfg)
from gaussian_wrangler.gw_common import (process_gausscom_file, create_com_from_pdb_str)
from gaussian_wrangler import __version__

__author__ = 'hmayes'


# Constants #

# Config File Sections

# Config keys
GAU_TPL_FILE = 'gau_tpl_file'
PDB_LIST_FILE = 'pdb_list_file'
PDB_FILE = 'pdb_file'
REMOVE_H = 'remove_final_h'
NUM = 'num'
DIH_ROT = "dih_rot"  # enter as list: four atom ids (base 1), then rotation in degrees
DIH_DATA = "dih_data"
MAX_CONF = "max_conf"
ORIGINAL = "original"

# Defaults
DEF_CFG_FILE = 'pdb2gau.ini'
DEF_CFG_VALS = {PDB_LIST_FILE: 'pdb_list.txt',
                PDB_FILE: None,
                REMOVE_H: False,
                NUM: None,
                DIH_ROT: None,
                MAX_CONF: 500,
                ORIGINAL: False,
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
        if main_proc[NUM]:
            main_proc[NUM] = int(main_proc[NUM])
    else:
        main_proc = {GAU_TPL_FILE: None, CONFIG_NAME: floc}
        for key, def_val in DEF_CFG_VALS.items():
            main_proc[key] = def_val

    main_proc[DIH_DATA] = []
    if main_proc[DIH_ROT] is not None:
        try:
            dih_list = main_proc[DIH_ROT].split(";")
            for dih in dih_list:
                dih_data = dih.split(",")
                if len(dih_data) != 5:
                    raise IndexError
                # note: RDKit is zero-based with atom indices, thus subtracting one from each number
                dih_data[:4] = [int(x) - 1 for x in dih_data[:4]]
                # noinspection PyTypeChecker
                dih_data[4] = float(dih_data[4])
                main_proc[DIH_DATA].append(dih_data)
        except (ValueError, IndexError):
            raise InvalidDataError("Error in parsing dihedral entry. Enter multiple dihedrals by separating data "
                                   "with a semicolon (';'). Each dihedral should be specified with 5 values, were the "
                                   "first four are one-based integer atom ids, and the last value is the rotation "
                                   "increment in degrees. ")

        if main_proc[MAX_CONF]:
            main_proc[MAX_CONF] = int(main_proc[MAX_CONF])
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
                                               "specify the '{}' (-t) and '{}' (-1) or '{}' (-f). The command lines "
                                               "for the '{}' flag (-r) or only the first entry in the pdb ('{}', -a) "
                                               "may also be specified.".format(DEF_CFG_FILE, GAU_TPL_FILE,
                                                                               PDB_LIST_FILE, PDB_FILE, REMOVE_H, NUM),
                        default=DEF_CFG_FILE, type=read_cfg)
    parser.add_argument("-t", "--tpl_file", help="Specifies the '{}'".format(GAU_TPL_FILE),
                        default=None)
    parser.add_argument("-l", "--pdb_list_file", help="Option to specify a file with a list of pdbs ('{}') to convert "
                                                      "(one file per line on the list).".format(PDB_LIST_FILE),
                        default=None)
    parser.add_argument("-f", "--file", help="Option to specify a pdb file ('{}') to convert.".format(PDB_FILE),
                        default=None)
    parser.add_argument("-n", "--num", help="Only read if a config file is not provided. This command can be used to "
                                            "specify only using the first '-n'/'--num' set(s) of coordinates in a pdb "
                                            "file to create gausscom file(s). The default is to use all coordinates, "
                                            "making as many input files as there are molecules/conformations in the "
                                            "pdb.", default=None, type=int)
    parser.add_argument("-r", "--remove_final_h", help="Option to specify removing the last H atom from the PDB "
                                                       "file(s) when creating the gausscom files. The default is "
                                                       "False.", action='store_true')
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
                if args.num:
                    args.config[NUM] = args.num
                if args.remove_final_h:
                    args.config[REMOVE_H] = True
                if args.file:
                    args.config[PDB_FILE] = args.file
                if args.pdb_list_file:
                    args.config[PDB_LIST_FILE] = args.pdb_list_file
    except (IOError, KeyError, InvalidDataError, MissingSectionHeaderError, SystemExit) as e:
        if hasattr(e, 'code') and e.code == 0:
            return args, GOOD_RET
        warning(e)
        parser.print_help()
        return args, INPUT_ERROR

    return args, GOOD_RET


def add_atom_indices(mol):
    for i, a in enumerate(mol.GetAtoms()):
        a.SetAtomMapNum(i)


def rotate_dihes_pdb_files(cfg, gau_tpl_content, pdb_files):
    """
    If dihedral data is specified, use RDkit to rotate specified dihedrals
    :param cfg: dict of configuration input
    :param gau_tpl_content: dict of data to create gaussian input files
    :param pdb_files: list of pdb files with dihedral angles to be rotated
    :return: n/a, saves new Gaussian input files
    """
    for pdb_file in pdb_files:
        mol_orig = MolFromPDBFile(pdb_file, removeHs=False)
        all_confs = [mol_orig]
        num_atoms = mol_orig.GetNumAtoms()
        for dih in cfg[DIH_DATA]:
            max_id = max(dih[:4])
            rot_deg = dih[4]
            if max_id > num_atoms:
                raise InvalidDataError(f"Dihedral rotation specifies an atom id of {max_id}, while only {num_atoms} "
                                       f"atoms were found in the PBD file {os.path.relpath(pdb_file)}")
            new_confs = []
            for mol_id, current_mol in enumerate(all_confs):
                # atoms = [a for a in current_mol.GetAtoms()]
                # for a in atoms:
                #     print(a.GetIdx(), a.GetSymbol())
                dih_deg = GetDihedralDeg(current_mol.GetConformer(0), *dih[:4])
                for _ in range(int(round(360. / rot_deg, 0) - 1)):
                    dih_deg = dih_deg + rot_deg
                    # print(dih, dih_deg)
                    SetDihedralDeg(current_mol.GetConformer(0), *dih[:4], dih_deg)
                    new_confs.append(current_mol.__copy__())
            all_confs.extend(new_confs)
        create_coms_from_mol_list(all_confs, gau_tpl_content, pdb_file, cfg[MAX_CONF], cfg[ORIGINAL])


def create_coms_from_mol_list(conformer_list, gau_tpl_content, base_out_name, max_num_coms, print_original):
    """
    From a list of RDKit mol objects, create gaussian output files, optionally for only the specified number of
    objects
    :param conformer_list:
    :param gau_tpl_content:
    :param base_out_name:
    :param max_num_coms: int or infinity
    :param print_original: Boolean, whether to print the initial conformation
    :return:
    """
    energy_list = []
    if print_original:
        start_at = 0
    else:
        start_at = 1

    RDLogger.DisableLog('rdApp.*')
    for current_mol in conformer_list[start_at:]:
        opt_results = MMFFOptimizeMoleculeConfs(current_mol, maxIters=0)
        energy_list.append(opt_results[0][1])

    combined_lists = zip(energy_list, conformer_list)
    zipped_sorted = sorted(combined_lists, key=itemgetter(0))

    # for energy in sorted(energy_list):
    #     print(f"{energy:15.8f}")
    mol_num = 0
    last_energy = np.nan
    print_note = False
    com_fname = None
    for energy, current_mol in zipped_sorted:
        if mol_num >= max_num_coms:
            if np.isclose(energy, last_energy):
                print_note = True
            else:
                break
        mol_num += 1
        last_energy = energy
        com_fname = create_out_fname(base_out_name, suffix=f"_{mol_num}", ext=".com", rel_path=True)
        pdb_str = MolToPDBBlock(current_mol)
        create_com_from_pdb_str(pdb_str, gau_tpl_content, com_fname)
        print(f"{int(energy):12,} {com_fname}")

    if com_fname:
        print(f"Wrote {mol_num} files, ending with: {os.path.relpath(com_fname)}")
    else:
        print("No output created from rotating dihedrals.")
    if print_note:
        print(f"More than {max_num_coms} conformations were output to ties calculated energies.")


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
    if cfg[DIH_DATA]:
        rotate_dihes_pdb_files(cfg, gau_tpl_content, pdb_files)
    else:
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
                if cfg[NUM] and mol_num >= cfg[NUM]:
                    return
                pdb_atom_line = []


def main(argv=None):
    print(f"Running GaussianWrangler script pdbs2gausscoms version {__version__}")

    # Read input
    args, ret = parse_cmdline(argv)
    if ret != GOOD_RET or args is None:
        return ret

    # Read pdb files
    cfg = args.config
    try:
        gau_tpl_content = process_gausscom_file(cfg[GAU_TPL_FILE])
        process_pdb_files(cfg, gau_tpl_content)
    except (IOError, UnicodeDecodeError) as e:
        warning("Problems reading file:", e)
        return IO_ERROR
    except InvalidDataError as e:
        warning("Problems reading data template:", e)
        return INVALID_DATA

    return GOOD_RET  # success


if __name__ == '__main__':
    status = main()
    sys.exit(status)
