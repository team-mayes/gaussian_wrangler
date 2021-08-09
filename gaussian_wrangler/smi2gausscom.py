#!/usr/bin/env python
"""
Creates Gaussian files from SMILES strings.
"""

import argparse
import os
import sys
import pubchempy as pcp
from configparser import MissingSectionHeaderError
from rdkit import Chem
from rdkit.Chem.AllChem import (EmbedMultipleConfs)
from rdkit.Chem.rdmolfiles import MolToPDBBlock
from rdkit.Chem.rdmolops import AddHs
from common_wrangler.common import (GOOD_RET, INPUT_ERROR, IO_ERROR, INVALID_DATA,
                                    InvalidDataError, warning, create_out_fname, OUT_DIR,
                                    file_rows_to_list, make_dir, read_tpl)
from common_wrangler.fill_tpl import fill_save_tpl
from gaussian_wrangler.gw_common import (get_pdb_coord_list)
from gaussian_wrangler import __version__

__author__ = 'hmayes'


# Constants #

# Config File Sections

# Config keys
GAU_TPL_FILE = 'gau_tpl_file'
SMI = "smiles"
GATHER = "gather"
LIST_FNAME = "list_file"
MAX_CONFS = 'max_confs'
DEF_MAX_CONFS = 50
ATOMS = "atoms"
REQ_STR = f"{{{ATOMS}}}"
SPECIAL_REPLACE = {"(": "p", ")": "q", "[": "s", "]": "t",
                   "=": "e", "@  ": "a", "#": "h", "'": "i",
                   ",": "c", " ": "_", ".": "p", "+": "w",
                   "/": "d", "\\": "b"}

# # If add config file...
# DEF_CFG_VALS = {SMI: None,
#                 LIST_FNAME: None,
#                 OUT_DIR: None,
#                 MAX_CONFS: DEF_MAX_CONFS,
#                 }
# REQ_KEYS = {GAU_TPL_FILE: str,
#             }


def parse_cmdline(argv):
    """
    Returns the parsed argument list and return code.
    `argv` is a list of arguments, or `None` for ``sys.argv[1:]``.
    """
    if argv is None:
        argv = sys.argv[1:]

    # initialize the parser object:
    parser = argparse.ArgumentParser(description='Creates Gaussian input files from SMILES strings, given a template '
                                                 'input file. The created files will be named based on their CID.')
    parser.add_argument("-t", "--" + GAU_TPL_FILE, help=f"Required: the location of the Gaussian input template file. "
                                                        f"This file must contain the string '{REQ_STR}' in the "
                                                        f"location where the atom type and coordinates should be "
                                                        f"added.",
                        metavar='PATH', default=None)
    parser.add_argument("-l", "--" + LIST_FNAME, help=f"Option to specify a file with a list of SMILES strings "
                                                      f"(one file per line on the list).",
                        metavar='PATH', default=None)
    parser.add_argument("-m", "--" + MAX_CONFS, help=f"Option to specify the maximum number of conformations to be "
                                                     f"created for each SMILES string. The default value is "
                                                     f"{DEF_MAX_CONFS}.",
                        metavar='INT', default=DEF_MAX_CONFS, type=int)
    parser.add_argument("-o", "--" + OUT_DIR, help="Directory where created files should be saved. The default "
                                                   "is the working directory. If a directory is specified and does not "
                                                   "yet exist, it will be created.", metavar='PATH', default=None)
    parser.add_argument("-s", "--" + SMI, help="Option to specify a SMILES string. Multiple strings can be separated "
                                               "by ','; if '.' is present the two molecules/fragments will be "
                                               "considered together.", metavar='STR', default=None)
    args = None
    try:
        args = parser.parse_args(argv)
        if args.gau_tpl_file is None:
            raise InvalidDataError(f"A template Gaussian input file is required to be specified with the "
                                   f"'-t'/'--{GAU_TPL_FILE}' option.")
        if not os.path.isfile(args.gau_tpl_file):
            raise InvalidDataError(f"Could not locate the specified Gaussian input file.")

        if args.list_file is None:
            args.list_file = []
        else:
            args.list_file = file_rows_to_list(args.list_file)
        if args.smiles is not None:
            smi_list = [smi.strip() for smi in args.smiles.split(",")]
            args.list_file.extend(smi_list)
        if len(args.list_file) == 0:
            raise InvalidDataError(f"No SMILES input has been specified. Specify a single SMILES string with the "
                                   f"'-s'/'--{SMI}' option or a files with a list of SMILES strings (one per line) "
                                   f"with the '-l'/'--{LIST_FNAME}' option.")
        args.list_file = list(set(args.list_file))

        if args.out_dir is not None:
            make_dir(args.out_dir)
    except (IOError, KeyError, InvalidDataError, MissingSectionHeaderError, SystemExit) as e:
        if hasattr(e, 'code') and e.code == 0:
            return args, GOOD_RET
        warning(e)
        parser.print_help()
        return args, INPUT_ERROR

    return args, GOOD_RET


def gen_conformers(mol, num_confs=DEF_MAX_CONFS, max_attempts=1000,
                   prune_rms_thresh=0.025,
                   use_exp_torsion_angle_prefs=True, use_basic_knowledge=True, enforce_chirality=True
                   ):
    ids = EmbedMultipleConfs(mol, numConfs=num_confs, maxAttempts=max_attempts,
                             pruneRmsThresh=prune_rms_thresh,
                             useExpTorsionAnglePrefs=use_exp_torsion_angle_prefs,
                             useBasicKnowledge=use_basic_knowledge,
                             enforceChirality=enforce_chirality, numThreads=0, randomSeed=10)
    return list(ids)


def get_mol_name(smi):
    """
    Get the compound ID from pubchem if the molecule is in pubchem.
    Otherwise, return a SMILES string with all special characters replaced.
    :param smi: str, a SMILES string
    :return: mol_name, a str without special characters
    """
    try:
        # cid instead?
        results = pcp.get_compounds(smi, namespace='smiles')
        if results and results[0].cid:
            mol_name = f"cid_{results[0].cid}"
        else:
            raise pcp.BadRequestError
    except pcp.BadRequestError:
        mol_name = smi
        for spec_char, char in SPECIAL_REPLACE.items():
            mol_name = mol_name.replace(spec_char, char)

    return mol_name


def process_smiles(gau_tpl_fname, smi_list, max_num_confs, out_dir):
    """
    Creates Gaussian input files for each SMILES string provided
    https://www.rdkit.org/docs/GettingStartedInPython.html
    :param smi_list: list of SMILES strings
    :param gau_tpl_fname: str, the location of the template file to use to create input files
    :param max_num_confs: int, the maximum number of conformations to generate
    :param out_dir: str, directory where files are to be saved (if None, saves to working directory)
    :return: N/A, writes files and prints notes on files created
    """
    gau_tpl_str = read_tpl(gau_tpl_fname)
    if REQ_STR not in gau_tpl_str:
        raise InvalidDataError(f"Did not find the required string '{REQ_STR}' in the provided Gaussian input "
                               f"template file.")
    for smi in smi_list:
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            warning(f"Skipping SMILES input string '{smi}' due to error\n")
            continue
        Chem.Kekulize(mol)
        mol = AddHs(mol)
        confs = gen_conformers(mol, num_confs=max_num_confs)
        cid = get_mol_name(smi)
        base_fname = create_out_fname(cid, ext='com', base_dir=out_dir, rel_path=True)
        conf_id = -1  # make IDE happy
        for conf_id in confs:
            com_fname = create_out_fname(base_fname, suffix=f'_{conf_id}')
            pdb_str = MolToPDBBlock(mol, confId=conf_id)
            coord_list = get_pdb_coord_list(pdb_str)
            fill_save_tpl(gau_tpl_str, {ATOMS: "\n".join(coord_list)}, gau_tpl_fname, com_fname, print_info=False)
        print(f"Wrote {conf_id + 1} files with base name '{base_fname}'")


def main(argv=None):
    print(f"Running GaussianWrangler script smi2gausscom version {__version__}")

    # Read input
    args, ret = parse_cmdline(argv)
    if ret != GOOD_RET or args is None:
        return ret

    try:
        process_smiles(args.gau_tpl_file, args.list_file, args.max_confs, args.out_dir)
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
