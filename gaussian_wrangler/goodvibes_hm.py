# !/usr/bin/env python
"""
                            GoodVibes.py
Evaluation of quasi-harmonic thermochemistry from Gaussian output.
Partition functions are evaluated from vibrational frequencies and rotational temperatures from the standard output.

The rigid-rotor harmonic oscillator approximation is used as standard for all frequencies above a cut-off value. Below
this, two treatments can be applied to entropic values:
  (a) low frequencies are shifted to the cut-off value (as per Cramer-Truhlar)
  (b) a free-rotor approximation is applied below the cut-off (as per Grimme). In this approach, a damping function
  interpolates between the RRHO and free-rotor entropy treatment of Svib to avoid a discontinuity.
Both approaches avoid infinitely large values of Svib as wavenumbers tend to zero. With a cut-off set to 0, the results
will be identical to standard values output by the Gaussian program.

Enthalpy values below the cutoff value are treated similarly to Grimme's method (as per Head-Gordon) where below the
cutoff value, a damping function is applied as the value approaches a value of 0.5RT, appropriate for zeolitic systems

The free energy can be evaluated for variable temperature, concentration, vibrational scaling factor, and with a haptic
correction of the translational entropy in different solvents, according to the amount of free space available.

A potential energy surface may be evaluated for a given set of structures or conformers, in which case a correction to
the free energy due to multiple conformers is applied. Enantiomeric excess, diastereomeric ratios and ddG can also be
calculated to show preference of stereoisomers.

Careful checks may be applied to compare variables between multiple files such as Gaussian version, solvation models,
levels of theory, charge and multiplicity, potential duplicate structures errors in potential linear molecules,
correct or incorrect transition states, and empirical dispersion models.

Authors:     Rob Paton, Ignacio Funes-Ardoiz Guilian Luchini, Juan V. Alegre-Requena, Yanfei Guan
Last modified:  July 22, 2019
Further modified by Heather Mayes, Dec 2019
"""

import os
import sys
import argparse
import time
import numpy as np
from datetime import datetime, timedelta
from gaussian_wrangler.vib_scale_factors import (GetOutData, CalcBBE)
from gaussian_wrangler.goodvibes_functions import (ALPHABET, output_pes_temp_interval, create_plot, output_rel_e_data,
                                                   calc_enantio_excess, get_boltz, output_cosmos_rs_interval, all_same,
                                                   print_check_fails)
from common_wrangler.common import (InvalidDataError, warning,
                                    GAS_CONSTANT, ATM_TO_KPA, AU_TO_J,
                                    GOOD_RET, INPUT_ERROR, INVALID_DATA, file_rows_to_list)
from gaussian_wrangler import __version__

# # Below are the values originally used by GoodVibes; very close to current output
# GAS_CONSTANT = 8.3144621  # J / K / mol
# KB = 1.3806488e-23  # J / K, BOLTZMANN_CONSTANT
# H = 6.62606957e-34  # J * s, PLANCK_CONSTANT
# AVOGADRO_CONST = 6.0221415e23  # 1 / mol, AVOGADRO_CONSTANT
# AMU_to_KG = 1.66053886E-27  # UNIT CONVERSION
# ATM_TO_KPA = 101.325
# AU_TO_J = 4.184 * 627.509541 * 1000.0  # UNIT CONVERSION, J_TO_AU

# # To make the output exactly match Gaussian's, use the values below instead importing them from common_wrangler.common
# GAS_CONSTANT = 8.31441  # J / K / mol; in common, GAS_CONSTANT = 8.314462618
# ATM_TO_KPA = 101.325  # 1 atm in kPa (no change)
# EHPART_TO_KCAL_MOL = 627.5095  # kcal/mol/(Eh/part); in common, the value is 627.5094709
# AU_TO_J = 4.184 * EHPART_TO_KCAL_MOL * 1000.0  # This value changes based on which EHPART_TO_KCAL_MOL is used


# VERSION NUMBER
__version__ = "3.0.1.hmayes"

SUPPORTED_EXTENSIONS = {'.out', '.log'}


GOODVIBES_REF = ("Luchini, G.; Alegre-Requena J. V.; Guan, Y.; Funes-Ardoiz, I.; Paton, R. S. (2019)."
                 "\n         http://doi.org/10.5281/zenodo.595246")

# Some useful arrays

# Symmetry numbers for different point groups
PG_SM = {"C1": 1, "Cs": 1, "Ci": 1, "C2": 2, "C3": 3, "C4": 4, "C5": 5, "C6": 6, "C7": 7, "C8": 8, "D2": 4, "D3": 6,
         "D4": 8, "D5": 10, "D6": 12, "D7": 14, "D8": 16, "C2v": 2, "C3v": 3, "C4v": 4, "C5v": 5, "C6v": 6, "C7v": 7,
         "C8v": 8, "C2h": 2, "C3h": 3, "C4h": 4, "C5h": 5, "C6h": 6, "C7h": 7, "C8h": 8, "D2h": 4, "D3h": 6, "D4h": 8,
         "D5h": 10, "D6h": 12, "D7h": 14, "D8h": 16, "D2d": 4, "D3d": 6, "D4d": 8, "D5d": 10, "D6d": 12, "D7d": 14,
         "D8d": 16, "S4": 4, "S6": 6, "S8": 8, "T": 6, "Th": 12, "Td": 12, "O": 12, "Oh": 24, "Cinfv": 1, "Dinfh": 2,
         "I": 30, "Ih": 60, "Kh": 1}

# solvent name, molecular weight, density (at 20 C)
SOLVENTS = {"water": (18.02, 0.998), "oxidane": (18.02, 0.998), "methanol": (32.04, 0.791),
            "acetonitrile": (41.052, 0.7857), "ethanol": (46.07, 0.789), "acetone": (58.079, 0.7845),
            "acetic_acid": (60.052, 1.0446), "1-propanol": (60.10, 0.803), "propan-1-ol": (60.10, 0.803),
            "2-propanol": (60.10, 0.803), "propan-2-ol": (60.10, 0.785), "nitromethane": (61.04, 1.382),
            "ethylene_glycol": (62.07, 1.115), "ethane-1,2-diol": (62.07, 1.115),
            "tetrahydrofuran": (72.106, 0.8833), "oxolane": (72.106, 0.8833),
            "2-butanone": (72.11, 0.7999), "butan-2-one": (72.11, 0.7999), "pentane": (72.15, 0.626),
            "n,n-dimethyl_formamide": (73.09, 0.9445),
            "1-butanol": (74.12, 0.8095), "butan-1-ol": (74.12, 0.8095),
            "2-butanol": (74.12, 0.8063), "butan-2-ol": (74.12, 0.8063),
            "tert-butyl_alcohol": (74.12, 0.7887), "2-methylpropan-2-ol": (74.12, 0.7887),
            "diethyl_ether": (74.12, 0.713), "ethoxyethane": (74.12, 0.713),
            "benzene": (78.11, 0.8765), "dimethyl_sulfoxide": (78.13, 1.092), "methylsulfinylmethane": (78.13, 1.092),
            "hexane": (86.18, 0.659), "1,4-dioxane": (88.11, 1.033), "ethyl_acetate": (88.11, 0.895),
            "tert-butyl_methyl_ether": (88.15, 0.741), "methyl_tert-butyl_ether": (88.15, 0.741),
            "2-methoxy-2-methylpropane": (88.15, 0.741),
            "glycerol": (92.09, 1.261), "glycerin": (92.09, 1.261), "propane-1,2,3-triol": (92.09, 1.261),
            "1,2-dimethoxyethane": (90.12, 0.8637), "toluene": (92.14, 0.867), "1,2-dichloroethane": (98.96, 1.245),
            "n-methyl-2-pyrrolidone": (99.13, 1.033), "1-methylpyrrolidin-2-one": (99.13, 1.033),
            "triethylamine": (101.19, 0.728), "n,n-diethylethanamine": (101.19, 0.728),
            "diethylene_glycol": (106.12, 1.1197), "2-(2-hydroxyethoxy)ethanol": (106.12, 1.1197),
            "m-xylene": (106.17, 0.868), "1,3-xylene": (106.17, 0.868), "o-xylene": (106.17, 0.868),
            "1,2-xylene": (106.17, 0.868), "p-xylene": (106.17, 0.897), "1,4-xylene": (106.17, 0.861),
            "chlorobenzene": (112.56, 1.1058), "trichloromethane": (119.38, 1.4788), "chloroform": (119.38, 1.4788),
            "diglyme": (134.17, 0.943), "1-methoxy-2-(2-methoxyethoxy)ethane": (134.17, 0.943),
            "tetrachloromethane": (153.82, 1.594), "carbon_tetrachloride": (153.82, 1.594),
            "hexamethyl_phosphorous_triamide": (163.20, 0.898),
            "n-[bis(dimethylamino)phosphanyl]-n-methylmethanamine": (163.20, 0.898),
            "hexamethylphosphoramide": (179.20, 1.03),
            "n-[bis(dimethylamino)phosphoryl]-n-methylmethanamine": (179.20, 1.03), }


def find_level_of_theory(file):
    # Read output for the level of theory and basis set used
    repeated_theory = 0
    with open(file) as f:
        data = f.readlines()
    level, bs = 'none', 'none'

    for line in data:
        if line.strip().find('External calculation') > -1:
            level, bs = 'ext', 'ext'
            break
        bs, level = find_freq_sp_dlpno_cbs(bs, level, line, repeated_theory)
    level_of_theory = '/'.join([level, bs])
    return level_of_theory


def find_freq_sp_dlpno_cbs(bs, level, line, repeated_theory):
    if '\\Freq\\' in line.strip() and repeated_theory == 0:
        try:
            level, bs = (line.strip().split("\\")[4:6])
            repeated_theory = 1
        except IndexError:
            pass
    elif '|Freq|' in line.strip() and repeated_theory == 0:
        try:
            level, bs = (line.strip().split("|")[4:6])
            repeated_theory = 1
        except IndexError:
            pass
    if '\\SP\\' in line.strip() and repeated_theory == 0:
        try:
            level, bs = (line.strip().split("\\")[4:6])
        except IndexError:
            pass
    elif '|SP|' in line.strip() and repeated_theory == 0:
        try:
            level, bs = (line.strip().split("|")[4:6])
        except IndexError:
            pass
    if 'DLPNO BASED TRIPLES CORRECTION' in line.strip():
        level = 'DLPNO-CCSD(T)'
    if 'Estimated CBS total energy' in line.strip():
        try:
            bs = ("Extrapol." + line.strip().split()[4])
        except IndexError:
            pass
    # Remove the restricted R or unrestricted U label
    if level[0] in ('R', 'U'):
        level = level[1:]
    return bs, level


# At beginning of procedure, read level of theory, solvation model, and check for normal termination
def read_initial(file):
    with open(file) as f:
        data = f.readlines()
    level, bs, program, keyword_line = 'none', 'none', 'none', 'none'
    progress, orientation = 'Incomplete', 'Input'
    a, repeated_theory = 0, 0
    no_grid = True
    dft, dft_used, level, bs, scf_iradan, cphf_iradan = False, 'F', 'none', 'none', False, False

    for line in data:
        # Determine program to find solvation model used
        if "Gaussian" in line:
            program = "Gaussian"
        if "* O   R   C   A *" in line:
            program = "Orca"
        # Grab pertinent information from file
        if line.strip().find('External calculation') > -1:
            level, bs = 'ext', 'ext'
        if line.strip().find('Standard orientation:') > -1:
            orientation = 'Standard'
        if line.strip().find('IExCor=') > -1 and no_grid:
            dft_used = line.split('=')[2].split()[0]
            no_grid = False
        bs, level = find_freq_sp_dlpno_cbs(bs, level, line, repeated_theory)

    # print(file,level)
    # Grab solvation models - Gaussian files
    solvation_model = None
    if program == 'Gaussian':
        for i, line in enumerate(data):
            if '#' in line.strip() and a == 0:
                for j, d_line in enumerate(data[i:i + 10]):
                    if '--' in d_line.strip():
                        a = a + 1
                        break
                    if a != 0:
                        break
                    else:
                        for k in range(len(d_line.strip().split("\n"))):
                            keyword_line += d_line.strip().split("\n")[k]
            if 'Normal termination' in line:
                progress = 'Normal'
            elif 'Error termination' in line:
                progress = 'Error'
        keyword_line = keyword_line.lower()
        if 'scrf' not in keyword_line.strip():
            solvation_model = "gas phase"
        else:
            start_scrf = keyword_line.strip().find('scrf') + 5
            if keyword_line[start_scrf] == "(":
                end_scrf = keyword_line.find(")", start_scrf)
                solvation_model = "scrf=" + keyword_line[start_scrf:end_scrf]
                if solvation_model[-1] != ")":
                    solvation_model = solvation_model + ")"
            else:
                start_scrf2 = keyword_line.strip().find('scrf') + 4
                if keyword_line.find(" ", start_scrf) > -1:
                    end_scrf = keyword_line.find(" ", start_scrf)
                else:
                    end_scrf = len(keyword_line)
                if keyword_line[start_scrf2] == "(":
                    solvation_model = "scrf=(" + keyword_line[start_scrf:end_scrf]
                    if solvation_model[-1] != ")":
                        solvation_model = solvation_model + ")"
                else:
                    if keyword_line.find(" ", start_scrf) > -1:
                        end_scrf = keyword_line.find(" ", start_scrf)
                    else:
                        end_scrf = len(keyword_line)
                    solvation_model = "scrf=" + keyword_line[start_scrf:end_scrf]
    # ORCA parsing for solvation model
    elif program == 'Orca':
        keyword_line_1 = "gas phase"
        keyword_line_2 = ''
        keyword_line_3 = ''
        for i, line in enumerate(data):
            if 'CPCM SOLVATION MODEL' in line.strip():
                keyword_line_1 = "CPCM,"
            if 'SMD CDS free energy correction energy' in line.strip():
                keyword_line_2 = "SMD,"
            if "Solvent:              " in line.strip():
                keyword_line_3 = line.strip().split()[-1]
            if 'ORCA TERMINATED NORMALLY' in line:
                progress = 'Normal'
            elif 'error termination' in line:
                progress = 'Error'
        solvation_model = keyword_line_1 + keyword_line_2 + keyword_line_3
    level_of_theory = '/'.join([level, bs])

    return level_of_theory, solvation_model, progress, orientation, dft_used


# Calculate elapsed time
def add_time(tm, cpu):
    [days, hrs, mins, secs, msecs] = cpu
    full_date = datetime(100, 1, tm.day, tm.hour, tm.minute, tm.second, tm.microsecond)
    full_date = full_date + timedelta(days=days, hours=hrs, minutes=mins, seconds=secs, microseconds=msecs * 1000)
    return full_date


# Check for duplicate species from among all files based on energy, rotational constants and frequencies
# Energy cutoff = 1 microHartree; RMS Rotational Constant cutoff = 1kHz; RMS Freq cutoff = 10 wavenumbers
def check_dup(files, thermo_data):
    e_cutoff = 1e-4
    ro_cutoff = 1e-4
    mae_freq_cutoff = 10
    max_freq_cutoff = 10
    dup_list = []
    freq_diff, mae_freq_diff, max_freq_diff, e_diff, ro_diff = 100, 3, 10, 1, 1
    for i, file in enumerate(files):
        for j in range(0, i):
            bbe_i, bbe_j = thermo_data[files[i]], thermo_data[files[j]]
            if hasattr(bbe_i, "scf_energy") and hasattr(bbe_j, "scf_energy"):
                e_diff = bbe_i.scf_energy - bbe_j.scf_energy
            if hasattr(bbe_i, "roconst") and hasattr(bbe_j, "roconst"):
                if len(bbe_i.roconst) == len(bbe_j.roconst):
                    ro_diff = np.linalg.norm(np.array(bbe_i.roconst) - np.array(bbe_j.roconst))
            if hasattr(bbe_i, "frequency_wn") and hasattr(bbe_j, "frequency_wn"):
                if len(bbe_i.frequency_wn) == len(bbe_j.frequency_wn) and len(bbe_i.frequency_wn) > 0:
                    freq_diff = [np.linalg.norm(freq_i - freq_j) for freq_i, freq_j in
                                 zip(bbe_i.frequency_wn, bbe_j.frequency_wn)]
                    mae_freq_diff, max_freq_diff = np.mean(freq_diff), np.max(freq_diff)
                elif len(bbe_i.frequency_wn) == len(bbe_j.frequency_wn) and len(bbe_i.frequency_wn) == 0:
                    mae_freq_diff, max_freq_diff = 0., 0.
            if e_diff < e_cutoff and ro_diff < ro_cutoff and mae_freq_diff < mae_freq_cutoff and \
                    max_freq_diff < max_freq_cutoff:
                dup_list.append([files[i], files[j]])
    return dup_list


def check_files(files, thermo_data, options, delimiter_row, l_o_t):
    # Perform careful checks on calculation output files
    # Check for Gaussian version, solvation state/gas phase consistency, level of theory/basis set consistency,
    # charge and multiplicity consistency, standard concentration used, potential linear molecule error,
    # transition state verification, empirical dispersion models.
    print("Checks for thermochemistry calculations (frequency calculations):")
    print(delimiter_row)
    # Check program used and version
    version_check = [thermo_data[key].version_program for key in thermo_data]
    file_check = [thermo_data[key].file for key in thermo_data]
    if all_same(version_check):
        print("o  Using {} in all calculations.".format(version_check[0]))
    else:
        print_check_fails(version_check, file_check, "programs or versions")

    # Check level of theory
    if all_same(l_o_t):
        print("o  Using {} in all calculations.".format(l_o_t[0]))
    else:
        print_check_fails(l_o_t, file_check, "levels of theory")

    # Check for solvent models
    solvent_check = [thermo_data[key].solvation_model[0] for key in thermo_data]
    if all_same(solvent_check):
        solvent_check = [thermo_data[key].solvation_model[1] for key in thermo_data]
        print("o  Using {} in all calculations.".format(solvent_check[0]))
    else:
        solvent_check = [thermo_data[key].solvation_model[1] for key in thermo_data]
        print_check_fails(solvent_check, file_check, "solvation models")

    # Check for -c 1 when solvent is added
    if all_same(solvent_check):
        if solvent_check[0] == "gas phase" and str(round(options.conc, 4)) == str(round(0.0408740470708, 4)):
            print("o  Using a standard concentration of 1 atm for gas phase.")
        elif solvent_check[0] == "gas phase" and str(round(options.conc, 4)) != str(round(0.0408740470708, 4)):
            print("x  Caution! Standard concentration is not 1 atm for gas phase (using {} M).".
                  format(options.conc))
        elif solvent_check[0] != "gas phase" and str(round(options.conc, 4)) == str(round(0.0408740470708, 4)):
            print("x  Using a standard concentration of 1 atm for solvent phase (option -c 1 should be "
                  "included for 1 M).")
        elif solvent_check[0] != "gas phase" and str(options.conc) == str(1.0):
            print("o  Using a standard concentration of 1 M for solvent phase.")
        elif solvent_check[0] != "gas phase" and str(round(options.conc, 4)) != str(round(0.0408740470708, 4)) and str(
                options.conc) != str(1.0):
            print("x  Caution! Standard concentration is not 1 M for solvent phase (using {} M).".format(options.conc))
    if not all_same(solvent_check) and "gas phase" in solvent_check:
        print("x  Caution! The right standard concentration cannot be determined because the calculations use a "
              "combination of gas and solvent phases.")
    if (not all_same(solvent_check)) and ("gas phase" not in solvent_check):
        print("x  Caution! Different solvents used, fix this issue and use option -c 1 for a standard "
              "concentration of 1 M.")

    # Check charge and multiplicity
    charge_check = [thermo_data[key].charge for key in thermo_data]
    multiplicity_check = [thermo_data[key].multiplicity for key in thermo_data]
    if all_same(charge_check) and all_same(multiplicity_check):
        print("o  Using charge {} and multiplicity {} in all calculations.".format(charge_check[0],
                                                                                   multiplicity_check[0]))
    else:
        print_check_fails(charge_check, file_check, "charge and multiplicity", multiplicity_check)

    # Check for duplicate structures
    dup_list = check_dup(files, thermo_data)
    if len(dup_list) == 0:
        print("o  No duplicates or enantiomers found")
    else:
        print("x  Caution! Possible duplicates or enantiomers found:")
        for dup in dup_list:
            print('        {} and {}'.format(dup[0], dup[1]))

    # Check for linear molecules with incorrect number of vibrational modes
    linear_fails, linear_fails_atom, linear_fails_cart, linear_fails_files, linear_fails_list = [], [], [], [], []
    frequency_list = []
    for file in files:
        linear_fails = GetOutData(file)
        linear_fails_cart.append(linear_fails.cartesians)
        linear_fails_atom.append(linear_fails.atom_types)
        linear_fails_files.append(file)
        frequency_list.append(thermo_data[file].frequency_wn)

    linear_fails_list.append(linear_fails_atom)
    linear_fails_list.append(linear_fails_cart)
    linear_fails_list.append(frequency_list)
    linear_fails_list.append(linear_fails_files)

    linear_mol_correct, linear_mol_wrong = [], []
    for i in range(len(linear_fails_list[0])):
        count_linear = 0
        if len(linear_fails_list[0][i]) == 2:
            if len(linear_fails_list[2][i]) == 1:
                linear_mol_correct.append(linear_fails_list[3][i])
            else:
                linear_mol_wrong.append(linear_fails_list[3][i])
        if len(linear_fails_list[0][i]) == 3:
            if linear_fails_list[0][i] == ['I', 'I', 'I'] or linear_fails_list[0][i] == ['O', 'O', 'O'] or \
                    linear_fails_list[0][i] == ['N', 'N', 'N'] or linear_fails_list[0][i] == ['H', 'C', 'N'] or \
                    linear_fails_list[0][i] == ['H', 'N', 'C'] or linear_fails_list[0][i] == ['C', 'H', 'N'] or \
                    linear_fails_list[0][i] == ['C', 'N', 'H'] or linear_fails_list[0][i] == ['N', 'H', 'C'] or \
                    linear_fails_list[0][i] == ['N', 'C', 'H']:
                if len(linear_fails_list[2][i]) == 4:
                    linear_mol_correct.append(linear_fails_list[3][i])
                else:
                    linear_mol_wrong.append(linear_fails_list[3][i])
            else:
                for j in range(len(linear_fails_list[0][i])):
                    for k in range(len(linear_fails_list[0][i])):
                        if k > j:
                            for m in range(len(linear_fails_list[1][i][j])):
                                if linear_fails_list[0][i][j] == linear_fails_list[0][i][k]:
                                    if (-linear_fails_list[1][i][k][m] - 0.1) < linear_fails_list[1][i][j][m] < \
                                            (-linear_fails_list[1][i][k][m] + 0.1):
                                        count_linear = count_linear + 1
                                        if count_linear == 3:
                                            if len(linear_fails_list[2][i]) == 4:
                                                linear_mol_correct.append(linear_fails_list[3][i])
                                            else:
                                                linear_mol_wrong.append(linear_fails_list[3][i])
        if len(linear_fails_list[0][i]) == 4:
            if linear_fails_list[0][i] == ['C', 'C', 'H', 'H'] or linear_fails_list[0][i] == ['C', 'H', 'C', 'H'] or \
                    linear_fails_list[0][i] == ['C', 'H', 'H', 'C'] or \
                    linear_fails_list[0][i] == ['H', 'C', 'C', 'H'] or \
                    linear_fails_list[0][i] == ['H', 'C', 'H', 'C'] or \
                    linear_fails_list[0][i] == ['H', 'H', 'C', 'C']:
                if len(linear_fails_list[2][i]) == 7:
                    linear_mol_correct.append(linear_fails_list[3][i])
                else:
                    linear_mol_wrong.append(linear_fails_list[3][i])
    linear_correct_print, linear_wrong_print = "", ""
    for i in range(len(linear_mol_correct)):
        linear_correct_print += ', ' + linear_mol_correct[i]
    for i in range(len(linear_mol_wrong)):
        linear_wrong_print += ', ' + linear_mol_wrong[i]
    linear_correct_print = linear_correct_print[1:]
    linear_wrong_print = linear_wrong_print[1:]
    if len(linear_mol_correct) == 0:
        if len(linear_mol_wrong) == 0:
            print("-  No linear molecules found.")
        if len(linear_mol_wrong) >= 1:
            print("x  Caution! Potential linear molecules with wrong number of frequencies found "
                  "(correct number = 3N-5) -{}.".format(linear_wrong_print))
    elif len(linear_mol_correct) >= 1:
        if len(linear_mol_wrong) == 0:
            print("o  All the linear molecules have the correct number of frequencies -{}.".format(
                linear_correct_print))
        if len(linear_mol_wrong) >= 1:
            print("x  Caution! Potential linear molecules with wrong number of frequencies found -{}. Correct "
                  "number of frequencies (3N-5) found in other calculations -{}.".format(linear_wrong_print,
                                                                                         linear_correct_print))

    # Checks whether any TS have > 1 imaginary frequency and any GS have any imaginary frequencies
    for file in files:
        bbe = thermo_data[file]
        base_name = os.path.basename(file)
        if bbe.job_type.find('TS') > -1 and len(bbe.im_frequency_wn) != 1:
            print("x  Caution! TS {} does not have 1 imaginary frequency greater than -50 wavenumbers.".
                  format(base_name))
        if bbe.job_type.find('GS') > -1 and bbe.job_type.find('TS') == -1 and len(bbe.im_frequency_wn) != 0:
            print("x  Caution: GS {} has 1 or more imaginary frequencies greater than -50 wavenumbers.".
                  format(base_name))

    # Check for empirical dispersion
    dispersion_check = [thermo_data[key].empirical_dispersion for key in thermo_data]
    if all_same(dispersion_check):
        if dispersion_check[0] == 'No empirical dispersion detected':
            print("-  No empirical dispersion detected in any of the calculations.")
        else:
            print("o  Using " + dispersion_check[0] + " in all calculations.")
    else:
        print_check_fails(dispersion_check, file_check, "dispersion models")
    print(delimiter_row)

    # Check for single-point corrections
    if options.spc:
        print("    Checks for single-point corrections:")
        print(delimiter_row)
        names_spc, version_check_spc = [], []
        for file in files:
            name, ext = os.path.splitext(file)
            if os.path.exists(name + '_' + options.spc + '.log'):
                names_spc.append(name + '_' + options.spc + '.log')
            elif os.path.exists(name + '_' + options.spc + '.out'):
                names_spc.append(name + '_' + options.spc + '.out')

        # Check SPC program versions
        version_check_spc = [thermo_data[key].sp_version_program for key in thermo_data]
        if all_same(version_check_spc):
            print("o  Using {} in all the single-point corrections.".format(version_check_spc[0]))
        else:
            print_check_fails(version_check_spc, file_check, "programs or versions")

        # Check SPC solvation
        solvent_check_spc = [thermo_data[key].sp_solvation_model for key in thermo_data]
        if all_same(solvent_check_spc):
            print("o  Using " + solvent_check_spc[0] + " in all the single-point corrections.")
        else:
            print_check_fails(solvent_check_spc, file_check, "solvation models")

        # Check SPC level of theory
        l_o_t_spc = [find_level_of_theory(name) for name in names_spc]
        if all_same(l_o_t_spc):
            print("o  Using {} in all the single-point corrections.".format(l_o_t_spc[0]))
        else:
            print_check_fails(l_o_t_spc, file_check, "levels of theory")

        # Check SPC charge and multiplicity
        charge_spc_check = [thermo_data[key].sp_charge for key in thermo_data]
        multiplicity_spc_check = [thermo_data[key].sp_multiplicity for key in thermo_data]
        if all_same(charge_spc_check) and all_same(multiplicity_spc_check):
            print("o  Using charge and multiplicity {} {} in all the single-point corrections.".
                  format(charge_spc_check[0], multiplicity_spc_check[0]))
        else:
            print_check_fails(charge_spc_check, file_check, "charge and multiplicity", multiplicity_spc_check)

        # Check if the geometries of freq calculations match their corresponding structures in single-point calculations
        geom_duplic_list, geom_duplic_list_spc, geom_duplic_cart, geom_duplic_files, geom_duplic_cart_spc, \
            geom_duplic_files_spc = [], [], [], [], [], []
        for file in files:
            geom_duplic = GetOutData(file)
            geom_duplic_cart.append(geom_duplic.cartesians)
            geom_duplic_files.append(file)
        geom_duplic_list.append(geom_duplic_cart)
        geom_duplic_list.append(geom_duplic_files)

        for name in names_spc:
            geom_duplic_spc = GetOutData(name)
            geom_duplic_cart_spc.append(geom_duplic_spc.cartesians)
            geom_duplic_files_spc.append(name)
        geom_duplic_list_spc.append(geom_duplic_cart_spc)
        geom_duplic_list_spc.append(geom_duplic_files_spc)
        spc_mismatching = "Caution! Potential differences found between frequency and single-point geometries -"
        if len(geom_duplic_list[0]) == len(geom_duplic_list_spc[0]):
            for i in range(len(files)):
                count = 1
                for j in range(len(geom_duplic_list[0][i])):
                    if count == 1:
                        if geom_duplic_list[0][i][j] == geom_duplic_list_spc[0][i][j]:
                            count = count
                        elif '{0:.3f}'.format(geom_duplic_list[0][i][j][0]) == '{0:.3f}'. \
                                format(geom_duplic_list_spc[0][i][j][0] * (-1)) or '{0:.3f}'. \
                                format(geom_duplic_list[0][i][j][0]) == '{0:.3f}'. \
                                format(geom_duplic_list_spc[0][i][j][0]):
                            if '{0:.3f}'.format(geom_duplic_list[0][i][j][1]) == '{0:.3f}'. \
                                    format(geom_duplic_list_spc[0][i][j][1] * (-1)) or '{0:.3f}'. \
                                    format(geom_duplic_list[0][i][j][1]) == '{0:.3f}'. \
                                    format(geom_duplic_list_spc[0][i][j][1] * (-1)):
                                count = count
                            if '{0:.3f}'.format(geom_duplic_list[0][i][j][2]) == '{0:.3f}'. \
                                    format(geom_duplic_list_spc[0][i][j][2] * (-1)) or '{0:.3f}'. \
                                    format(geom_duplic_list[0][i][j][2]) == '{0:.3f}'. \
                                    format(geom_duplic_list_spc[0][i][j][2] * (-1)):
                                count = count
                        else:
                            spc_mismatching += ", " + geom_duplic_list[1][i]
                            count = count + 1
            if spc_mismatching == "Caution! Potential differences found between frequency and single-point " \
                                  "geometries -":
                print("o  No potential differences found between frequency and single-point geometries (based "
                      "on input coordinates).")
            else:
                spc_mismatching_1 = spc_mismatching[:84]
                spc_mismatching_2 = spc_mismatching[85:]
                print("x  " + spc_mismatching_1 + spc_mismatching_2 + '.')
        else:
            print("x  One or more geometries from single-point corrections are missing.")

        # Check for SPC dispersion models
        dispersion_check_spc = [thermo_data[key].sp_empirical_dispersion for key in thermo_data]
        if all_same(dispersion_check_spc):
            if dispersion_check_spc[0] == 'No empirical dispersion detected':
                print("-  No empirical dispersion detected in any of the calculations.")
            else:
                print("o  Using " + dispersion_check_spc[0] + " in all the singe-point calculations.")
        else:
            print_check_fails(dispersion_check_spc, file_check, "dispersion models")
        print(delimiter_row)


def parse_cmdline(argv):
    """
    Returns the parsed argument list and return code.
    `argv` is a list of arguments, or `None` for ``sys.argv[1:]``.
    """
    if argv is None:
        argv = sys.argv[1:]

    # initialize the parser object:
    parser = argparse.ArgumentParser(description='This script is based on https://github.com/bobbypaton/GoodVibes, '
                                                 'copied here to allow tailored input and output.')
    parser.add_argument("-c", dest="conc", default=False, type=float, metavar="CONC",
                        help="Concentration (mol/l) (default 1 atm)")
    parser.add_argument("-f", dest="freq_cutoff", default=100, type=float, metavar="FREQ_CUTOFF",
                        help="Cut-off frequency for both entropy and enthalpy (wavenumbers) (default = 100)", )
    parser.add_argument("--fs", dest="S_freq_cutoff", default=100.0, type=float, metavar="S_FREQ_CUTOFF",
                        help="Cut-off frequency for entropy (wavenumbers) (default = 100)")
    parser.add_argument("--fh", dest="h_freq_cutoff", default=100.0, type=float, metavar="h_freq_cutoff",
                        help="Cut-off frequency for enthalpy (wavenumbers) (default = 100)")
    parser.add_argument("-l", dest="file_list", default=None, type=str, metavar="LIST",
                        help="List of file names to process, with one file per line.",)
    parser.add_argument("-q", dest="Q", action="store_true", default=False,
                        help="Quasi-harmonic entropy correction and enthalpy correction applied (default S=Grimme, "
                             "H=Head-Gordon)")
    parser.add_argument("--qs", dest="qs", default="grimme", type=str.lower, metavar="qs",
                        choices=('grimme', 'truhlar'),
                        help="Type of quasi-harmonic entropy correction (Grimme or Truhlar) (default Grimme)", )
    parser.add_argument("--qh", dest="qh", action="store_true", default=False,
                        help="Type of quasi-harmonic enthalpy correction (Head-Gordon)")
    parser.add_argument("-t", dest="temperature", default=298.15, type=float, metavar="TEMP",
                        help="Temperature (K) (default 298.15)")
    parser.add_argument("--ti", dest="temperature_interval", default=False, metavar="TI",
                        help="Initial temp, final temp, step size (K)")
    parser.add_argument("-v", dest="freq_scale_factor", default=False, type=float, metavar="SCALE_FACTOR",
                        help="Frequency scaling factor. If not set, the program will look for a suitable value in "
                             "its database. If not found, the program will use 1.0.")
    parser.add_argument("-z", dest="zpe_scale_factor", default=False, type=float, metavar="ZPE_SCALE_FACTOR",
                        help="Frequency scaling factor for calculation of the zero-point energy correction. "
                             "If not set, the same value will be used as for the frequency scaling factor "
                             "('-v' option) used for estimating entropy.")
    parser.add_argument("--vmm", dest="mm_freq_scale_factor", default=False, type=float, metavar="MM_SCALE_FACTOR",
                        help="Additional frequency scaling factor used in ONIOM calculations")
    parser.add_argument("--spc", dest="spc",
                        help="Indicates single point corrections (default False)", action="store_true", default=False)
    parser.add_argument("--boltz", dest="boltz", action="store_true", default=False,
                        help="Show Boltzmann factors")
    parser.add_argument("--cpu", dest="cputime", action="store_true", default=False,
                        help="Total CPU time")
    parser.add_argument("--d3", dest="D3", action="store_true", default=False,
                        help="Zero-damped DFTD3 correction will be computed")
    parser.add_argument("--d3bj", dest="D3BJ", action="store_true", default=False,
                        help="Becke-Johnson damped DFTD3 correction will be computed")
    parser.add_argument("--atm", dest="ATM", action="store_true", default=False,
                        help="Axilrod-Teller-Muto 3-body dispersion correction will be computed")
    parser.add_argument("--imag", dest="imag_freq", action="store_true", default=False,
                        help="Print imaginary frequencies (default False)")
    parser.add_argument("--invertifreq", dest="invert", nargs='?', const=True, default=False,
                        help="Make low lying imaginary frequencies positive (cutoff > -50.0 wavenumbers)")
    parser.add_argument("--freespace", dest="freespace", default="none", type=str, metavar="FREESPACE",
                        help="Solvent (H2O, toluene, DMF, AcOH, chloroform) (default none)")
    parser.add_argument("--dup", dest="duplicate", action="store_true", default=False,
                        help="Remove possible duplicates from thermochemical analysis")
    parser.add_argument("--cosmo", dest="cosmo", default=False, metavar="COSMO-RS",
                        help="Filename of a COSMO-RS .tab output file")
    parser.add_argument("--cosmo_int", dest="cosmo_int", default=False, metavar="COSMO-RS",
                        help="Filename of a COSMO-RS .tab output file along with a temperature range (K): "
                             "file.tab,'Initial_T, Final_T'")
    parser.add_argument("--pes", dest="pes", default=False, metavar="PES",
                        help="Tabulate relative values")
    parser.add_argument("--nogconf", dest="gconf", action="store_false", default=True,
                        help="Calculate a free-energy correction related to multi-configurational space (default "
                             "calculate Gconf)")
    parser.add_argument("--ee", dest="ee", default=False, type=str,
                        help="Tabulate selectivity values (excess, ratio) from a mixture, provide pattern for two "
                             "types such as *_R*,*_S*")
    parser.add_argument("--check", dest="check", action="store_true", default=False,
                        help="Checks if calculations were done with the same program, level of theory and solvent, "
                             "as well as detects potential duplicates")
    parser.add_argument("--media", dest="media", default=False, metavar="MEDIA",
                        help="Entropy correction for standard concentration of solvents")
    parser.add_argument("--custom_ext", type=str, default='',
                        help="List of additional file extensions to support, beyond .log or .out, use separated by "
                             "commas (ie, '.qfi, .gaussian'). It can also be specified with environment variable "
                             "GOODVIBES_CUSTOM_EXT")
    parser.add_argument("--graph", dest='graph', default=False, metavar="GRAPH",
                        help="Graph a reaction profile based on free energies calculated. ")
    parser.add_argument("--ssymm", dest='ssymm', action="store_true", default=False,
                        help="Turn on the symmetry correction.")

    args = None
    try:
        args = parser.parse_known_args(argv)

    except (SystemExit, ValueError) as e:
        if hasattr(e, 'code') and e.code == 0:
            return args, GOOD_RET
        warning(e)
        parser.print_help()
        return args, INPUT_ERROR

    return args, GOOD_RET


def output_pes_data(options, thermo_data, delimiter_row, interval, interval_bbe_data, interval_thermo_data, file_list):
    if options.gconf:
        print('\n   Gconf correction requested to be applied to below relative values using quasi-harmonic '
              'Boltzmann factors\n')
    for key in thermo_data:
        if not hasattr(thermo_data[key], "qh_gibbs_free_energy"):
            raise InvalidDataError("\nWarning! Could not find thermodynamic data for " + key + "\n")
        if (not hasattr(thermo_data[key], "sp_energy")) and options.spc:
            raise InvalidDataError("\nWarning! Could not find thermodynamic data for " + key + "\n")

    if options.temperature_interval:
        output_pes_temp_interval(options, delimiter_row, interval, interval_bbe_data, interval_thermo_data,
                                 file_list)
    else:
        # Output the relative energy data
        output_rel_e_data(options, delimiter_row, thermo_data)


def variable_temp_analysis(options, delimiter_row, files, t_interval, interval_bbe_data, gas_phase):
    print("Variable-Temperature analysis of the enthalpy, entropy and the entropy at a constant "
          "pressure between")
    if options.cosmo_int:
        interval = t_interval
        print("    T init:  {:.2f},   T final: {:.2f}\n".format(interval[0], interval[-1]))
    else:
        temperature_interval = [float(temp) for temp in options.temperature_interval.split(',')]
        # If no temperature step was defined, divide the region into 10
        if len(temperature_interval) == 2:
            temperature_interval.append((temperature_interval[1] - temperature_interval[0]) / 10.0)
        # below assumes that the interval is great than 1; no big deal if it isn't
        interval = np.arange(float(temperature_interval[0]), float(temperature_interval[1] + 1),
                             float(temperature_interval[2]))
        print("    T init:  {:.2f},  T final:  {:.2f},  T interval: {:.2f}\n".
              format(temperature_interval[0], temperature_interval[1], temperature_interval[2]))

    if options.qh:
        qh_print_format = "{:<39} {:>13} {:>24} {:>13} {:>10} {:>10} {:>13} {:>13}"
        if options.spc and options.cosmo_int:
            print(qh_print_format.format("Structure", "Temp/K", "H_SPC", "qh-H_SPC", "T.S", "T.qh-S",
                                         "G(T)_SPC", "COSMO-RS-qh-G(T)_SPC"))
        elif options.cosmo_int:
            print(qh_print_format.format("Structure", "Temp/K", "H", "qh-H", "T.S", "T.qh-S", "G(T)",
                                         "qh-G(T)", "COSMO-RS-qh-G(T)"))
        elif options.spc:
            print(qh_print_format.format("Structure", "Temp/K", "H_SPC", "qh-H_SPC", "T.S", "T.qh-S",
                                         "G(T)_SPC", "qh-G(T)_SPC"))
        else:
            print(qh_print_format.format("Structure", "Temp/K", "H", "qh-H", "T.S", "T.qh-S", "G(T)", "qh-G(T)"))
    else:
        print_format_3 = '{:<39} {:>13} {:>24} {:>10} {:>10} {:>13} {:>13}'
        if options.spc and options.cosmo_int:
            print(print_format_3.format("Structure", "Temp/K", "H_SPC", "T.S", "T.qh-S", "G(T)_SPC",
                                        "COSMO-RS-qh-G(T)_SPC"))
        elif options.cosmo_int:
            print(print_format_3.format("Structure", "Temp/K", "H", "T.S", "T.qh-S", "G(T)", "qh-G(T)",
                                        "COSMO-RS-qh-G(T)"))
        elif options.spc:
            print(print_format_3.format("Structure", "Temp/K", "H_SPC", "T.S", "T.qh-S", "G(T)_SPC", "qh-G(T)_SPC"))
        else:
            print(print_format_3.format("Structure", "Temp/K", "H", "T.S", "T.qh-S", "G(T)", "qh-G(T)"))

    for h, file in enumerate(files):  # Temperature interval
        bbe = None  # Add because it is possible for this not to be defined
        print(delimiter_row)
        base_name = os.path.basename(file)
        name_str = '{:<39}'.format(base_name)
        interval_bbe_data.append([])
        for i in range(len(interval)):  # Iterate through the temperature range
            temp = interval[i]
            if gas_phase:
                conc = ATM_TO_KPA / GAS_CONSTANT / temp
            else:
                conc = options.conc
            linear_warning = []
            if options.cosmo_int:
                # haven't implemented D3 for this option
                pass
            else:
                bbe = CalcBBE(file, options.qs, options.qh, options.S_freq_cutoff, options.h_freq_cutoff, temp,
                              conc, options.freq_scale_factor, options.zpe_scale_factor, options.freespace,
                              options.spc, options.invert, 0.0, cosmo=False)
            interval_bbe_data[h].append(bbe)
            linear_warning.append(bbe.linear_warning)
            if linear_warning == [['Warning! Potential invalid calculation of linear molecule from Gaussian.']]:
                print("x  {}".format(name_str))
                print('          Warning! Potential invalid calculation of linear molecule from Gaussian ...')
            else:
                # Gaussian spc files
                if hasattr(bbe, "scf_energy") and not hasattr(bbe, "gibbs_free_energy"):
                    print("x  {}".format(name_str))
                # ORCA spc files
                elif not hasattr(bbe, "scf_energy") and not hasattr(bbe, "gibbs_free_energy"):
                    print("x  {}".format(name_str))
                if not hasattr(bbe, "gibbs_free_energy"):
                    print("Warning! Couldn't find frequency information ...")
                else:
                    name_temp = '{:<39} {:13.2f}'.format(base_name, temp)
                    if not options.media:
                        if all(getattr(bbe, attrib) for attrib in
                               ["enthalpy", "entropy", "qh_entropy", "gibbs_free_energy",
                                "qh_gibbs_free_energy"]):
                            if options.qh:
                                if options.cosmo_int:
                                    print('{} {:24.6f} {:13.6f} {:10.6f} {:10.6f} {:13.6f} {:13.6f}'.format(
                                        name_temp, bbe.enthalpy, bbe.qh_enthalpy, (temp * bbe.entropy),
                                        (temp * bbe.qh_entropy), bbe.gibbs_free_energy, bbe.cosmo_qhg))
                                else:
                                    print('{} {:24.6f} {:13.6f} {:10.6f} {:10.6f} {:13.6f} {:13.6f}'.format(
                                        name_temp, bbe.enthalpy, bbe.qh_enthalpy, (temp * bbe.entropy),
                                        (temp * bbe.qh_entropy), bbe.gibbs_free_energy,
                                        bbe.qh_gibbs_free_energy))
                            else:
                                if options.cosmo_int:
                                    print('{} {:24.6f} {:10.6f} {:10.6f} {:13.6f} {:13.6f}'.
                                          format(name_temp, bbe.enthalpy, (temp * bbe.entropy),
                                                 (temp * bbe.qh_entropy), bbe.gibbs_free_energy, bbe.cosmo_qhg))
                                else:
                                    print('{} {:24.6f} {:10.6f} {:10.6f} {:13.6f} {:13.6f}'.
                                          format(name_temp, bbe.enthalpy, (temp * bbe.entropy),
                                                 (temp * bbe.qh_entropy), bbe.gibbs_free_energy,
                                                 bbe.qh_gibbs_free_energy))
                    else:
                        if options.media.lower() in SOLVENTS and options.media.lower() == \
                                os.path.splitext(os.path.basename(file))[0].lower():
                            mw_solvent = SOLVENTS[options.media.lower()][0]
                            density_solvent = SOLVENTS[options.media.lower()][1]
                            concentration_solvent = (density_solvent * 1000) / mw_solvent
                            media_correction = -(GAS_CONSTANT / AU_TO_J) * np.log(concentration_solvent)
                            if all(getattr(bbe, attrib) for attrib in
                                   ["enthalpy", "entropy", "qh_entropy", "gibbs_free_energy",
                                    "qh_gibbs_free_energy"]):
                                if options.qh:
                                    print('{} {:10.6f} {:13.6f} {:13.6f} {:10.6f} {:10.6f} {:13.6f} '
                                          '{:13.6f}'.format(name_temp, bbe.zpe, bbe.enthalpy, bbe.qh_enthalpy,
                                                            (temp * (bbe.entropy + media_correction)),
                                                            (temp * (bbe.qh_entropy + media_correction)),
                                                            bbe.gibbs_free_energy + (temp * (-media_correction)),
                                                            bbe.qh_gibbs_free_energy + (temp * (-media_correction))))
                                    print("  Solvent")
                            else:
                                print('{} {:10.6f} {:13.6f} {:10.6f} {:10.6f} {:13.6f} '
                                      '{:13.6f}'.format(name_temp, bbe.zpe, bbe.enthalpy,
                                                        (temp * (bbe.entropy + media_correction)),
                                                        (temp * (bbe.qh_entropy + media_correction)),
                                                        bbe.gibbs_free_energy + (temp * (-media_correction)),
                                                        bbe.qh_gibbs_free_energy + (temp * (-media_correction))))
                                print("  Solvent")
                        else:
                            if all(getattr(bbe, attrib) for attrib in
                                   ["enthalpy", "entropy", "qh_entropy", "gibbs_free_energy", "qh_gibbs_free_energy"]):
                                if options.qh:
                                    print('{} {:10.6f} {:13.6f} {:13.6f} {:10.6f} {:10.6f} {:13.6f} '
                                          '{:13.6f}'.format(name_temp, bbe.zpe, bbe.enthalpy, bbe.qh_enthalpy,
                                                            (temp * bbe.entropy), (temp * bbe.qh_entropy),
                                                            bbe.gibbs_free_energy, bbe.qh_gibbs_free_energy))
                                else:
                                    print('{} {:10.6f} {:13.6f} {:10.6f} {:10.6f} {:13.6f} '
                                          '{:13.6f}'.format(name_temp, bbe.zpe, bbe.enthalpy,
                                                            (temp * bbe.entropy), (temp * bbe.qh_entropy),
                                                            bbe.gibbs_free_energy, bbe.qh_gibbs_free_energy))
        print(delimiter_row)


def main(argv=None):
    # Read input
    args, ret = parse_cmdline(argv)
    if ret != GOOD_RET or args is None:
        return ret

    (options, args) = args
    files = []
    clusters = []
    command = 'Requested: '
    clustering = False

    try:
        # Start printing results
        start_time = time.strftime("%Y/%m/%d %H:%M:%S", time.localtime())
        print("GoodVibes v{} {}\n    REF: {}\n".format(__version__, start_time, GOODVIBES_REF))

        # If requested, turn on head-gordon enthalpy correction
        if options.Q:
            options.qh = True
        if options.qh:
            delimiter_row = "-" * 142
        else:
            delimiter_row = "-" * 128
        # If user has specified different file extensions
        if options.custom_ext or os.environ.get('GOODVIBES_CUSTOM_EXT', ''):
            custom_extensions = options.custom_ext.split(',') + os.environ.get('GOODVIBES_CUSTOM_EXT', '').split(',')
            for ext in custom_extensions:
                SUPPORTED_EXTENSIONS.add(ext.strip())

        # Default value for inverting imaginary frequencies
        if options.invert:
            options.invert = -50.0
        elif options.invert > 0:
            options.invert = -1 * options.invert

        # Initialize the total CPU time
        total_cpu_time = datetime(100, 1, 1, 00, 00, 00, 00)
        add_days = 0
        if len(args) > 1:
            for elem in args:
                if elem == 'clust:':
                    clustering = True
                    options.boltz = True
        # Get the filenames from the command line prompt
        # add any that come from the list command
        if options.file_list:
            file_list = file_rows_to_list(options.file_list)
            args += file_list
        missing_files = set()
        for elem in args:
            if clustering:
                if elem == 'clust:':
                    clusters.append([])
            try:
                # Look for file names
                if os.path.splitext(elem)[1].lower() in SUPPORTED_EXTENSIONS:
                    if os.path.isfile(elem):
                        # skip repeats
                        if elem not in files:
                            if options.spc:
                                name, ext = os.path.splitext(os.path.relpath(elem))
                                if os.path.exists(name + '_spc.log') or os.path.exists(name + '_spc.out'):
                                    files.append(elem)
                                else:
                                    raise InvalidDataError("SPC calculation file '{}.{}' not found!\n    Make sure "
                                                           "files are named with the convention: "
                                                           "'filename_spc.{}'.\n    For help, use option '-h'\n"
                                                           .format(name + '_spc', ext, ext))
                            else:
                                files.append(elem)
                    else:
                        missing_files.add(elem)
                elif elem != 'clust:':  # Look for requested options
                    command += elem + ' '
            except IndexError:
                pass

        if len(missing_files) > 0:
            raise IOError("Could not find the following file(s):\n    {}".format("\n    ".join(missing_files)))

        # Check if user has specified any files, if not quit now
        if len(files) == 0:
            raise InvalidDataError("No file names provided. Please provide GoodVibes with calculation output files on "
                                   "the command line.\n    For help, use option '-h'\n")
        if clustering:
            command += '(clustering active)'
        print(command)
        if not options.temperature_interval:
            print("    Temperature = " + str(options.temperature) + " Kelvin")
        # If not at standard temp, need to correct the molarity of 1 atmosphere (assuming pressure is still 1 atm)
        if options.conc:
            gas_phase = False
            print("    Concentration = " + str(options.conc) + " mol/L")
        else:
            gas_phase = True
            options.conc = ATM_TO_KPA / (GAS_CONSTANT * options.temperature)
            print("    Pressure = 1 atm")
        print('\nAll energetic values below shown in Hartree unless otherwise specified.')
        # Initial read of files,
        # Grab level of theory, solvation model, check for Normal Termination
        l_o_t, s_m, progress, orientation, grid = [], [], {}, {}, {}
        for file in files:
            lot_sm_prog = read_initial(file)
            l_o_t.append(lot_sm_prog[0])
            s_m.append(lot_sm_prog[1])
            progress[file] = lot_sm_prog[2]
            orientation[file] = lot_sm_prog[3]
            grid[file] = lot_sm_prog[4]
        remove_key = []
        # Remove problem files and print errors
        for i, key in enumerate(files):
            if progress[key] == 'Error':
                warning("Error termination found in file: {}\n    This file will be omitted from further "
                        "calculations.".format(key))
                remove_key.append([i, key])
            elif progress[key] == 'Incomplete':
                warning("File {} may not have terminated normally or the calculation may still be "
                        "running. This file will be omitted from further calculations.".format(key))
                remove_key.append([i, key])
        for [i, key] in list(reversed(remove_key)):
            files.remove(key)
            del l_o_t[i]
            del s_m[i]
            del orientation[key]
            del grid[key]
        if len(files) == 0:
            raise InvalidDataError("Please try again with normally terminated output files.\n    For help, "
                                   "use option '-h'\n")

        cosmo_solv, gsolv_dicts, t_interval = output_cosmos_rs_interval(files, options, s_m, l_o_t)

        # Check for special options
        if options.ssymm:
            ssymm_option = options.ssymm
        else:
            ssymm_option = False
        if options.mm_freq_scale_factor:
            vmm_option = options.mm_freq_scale_factor
        else:
            vmm_option = False

        # Loop over all specified output files and compute thermochemistry
        thermo_data = compute_thermochem(files, options,
                                         cosmo_solv=cosmo_solv, ssymm_option=ssymm_option, vmm_option=vmm_option)
        interval_bbe_data, interval_thermo_data = [], []

        inverted_freqs, inverted_files = [], []
        for file in files:
            if len(thermo_data[file].inverted_freqs) > 0:
                inverted_freqs.append(thermo_data[file].inverted_freqs)
                inverted_files.append(file)

        # Check if user has chosen to make any low lying imaginary frequencies positive
        if options.invert:
            for i, file in enumerate(inverted_files):
                if len(inverted_freqs[i]) == 1:
                    print("\n   The following frequency was made positive and used in calculations: " +
                          str(inverted_freqs[i][0]) + " from " + file)
                elif len(inverted_freqs[i]) > 1:
                    print("\n   The following frequencies were made positive and used in calculations: " +
                          str(inverted_freqs[i]) + " from " + file)

        # Adjust printing according to options requested
        if options.spc:
            delimiter_row += '-' * 14
        if options.cosmo:
            delimiter_row += '-' * 30
        if options.imag_freq:
            delimiter_row += '-' * 9
        if options.boltz:
            delimiter_row += '-' * 7
        if options.ssymm:
            delimiter_row += '-' * 13

        # Perform checks for consistent options provided in calculation files (level of theory)
        if options.check:
            check_files(files, thermo_data, options, delimiter_row, l_o_t)

        # Standard mode: tabulate thermochemistry output from file(s) at a single temperature and concentration
        interval = None
        dup_list = []
        # Running a variable temperature analysis of the enthalpy, entropy and the free energy
        if options.temperature_interval:
            variable_temp_analysis(options, delimiter_row, files, t_interval, interval_bbe_data, gas_phase)
        else:
            if options.spc:
                print("\n")
                if options.qh:
                    print('{:<39} {:>13} {:>13} {:>10} {:>13} {:>13} {:>10} {:>10} {:>13} {:>13}'.
                          format("Structure", "E_SPC", "E", "ZPE", "H_SPC", "qh-H_SPC", "T.S", "T.qh-S", "G(T)_SPC",
                                 "qh-G(T)_SPC"))
                else:
                    print('{:<39} {:>13} {:>13} {:>10} {:>13} {:>10} {:>10} {:>13} {:>13}'.
                          format("Structure", "E_SPC", "E", "ZPE", "H_SPC", "T.S", "T.qh-S", "G(T)_SPC", "qh-G(T)_SPC"))
            else:
                if options.qh:
                    print('{:<39} {:>13} {:>10} {:>13} {:>13} {:>10} {:>10} {:>13} '
                          '{:>13}'.format("Structure", "E", "ZPE", "H", "qh-H", "T.S", "T.qh-S", "G(T)", "qh-G(T)"))
                else:
                    print('{:<39} {:>13} {:>10} {:>13} {:>10} {:>10} {:>13} {:>13}'.
                          format("Structure", "E", "ZPE", "H", "T.S", "T.qh-S", "G(T)", "qh-G(T)"))

            if options.cosmo:
                print('{:>13} {:>16}'.format("COSMO-RS", "COSMO-qh-G(T)"))
            if options.boltz:
                print('{:>7}'.format("Boltz"))
            if options.imag_freq:
                print('{:>9}'.format("im freq"))
            if options.ssymm:
                print('{:>13}'.format("Point Group"))
            print(delimiter_row)
            # Look for duplicates or enantiomers
            if options.duplicate:
                dup_list = check_dup(files, thermo_data)

            # Boltzmann factors and averaging over clusters
            boltz_facts, boltz_sum, weighted_free_energy = None, None, None  # make IDE happy
            if options.boltz:
                boltz_facts, weighted_free_energy, boltz_sum = get_boltz(files, thermo_data, clustering, clusters,
                                                                         options.temperature, dup_list)

            for file in files:  # Loop over the output files and compute thermochemistry
                duplicate = False
                base_name = os.path.basename(file)
                name_str = '{:<39}'.format(base_name)
                if len(dup_list) != 0:
                    for dup in dup_list:
                        if dup[0] == file:
                            duplicate = True
                            print('\nx  {} is a duplicate or enantiomer of {}'.format(dup[0].rsplit('.', 1)[0],
                                                                                      dup[1].rsplit('.', 1)[0]))
                            break
                if not duplicate:
                    bbe = thermo_data[file]
                    if options.cputime:  # Add up CPU times
                        if hasattr(bbe, "cpu"):
                            if bbe.cpu is not None:
                                total_cpu_time = add_time(total_cpu_time, bbe.cpu)
                        if hasattr(bbe, "sp_cpu"):
                            if bbe.sp_cpu is not None:
                                total_cpu_time = add_time(total_cpu_time, bbe.sp_cpu)
                    if total_cpu_time.month > 1:
                        add_days += 31

                    # Check for possible error in Gaussian calculation of linear molecules which can return 2
                    # rotational constants instead of 3
                    if bbe.linear_warning:
                        print("\nx  " + name_str)
                        print('          ----   Caution! Potential invalid calculation of linear molecule '
                              'from Gaussian')
                    else:
                        e_str = ''
                        if hasattr(bbe, "gibbs_free_energy"):
                            if options.spc:
                                if bbe.sp_energy != '!':
                                    print('{} {:13.6f}'.format(name_str, bbe.sp_energy))
                                if bbe.sp_energy == '!':
                                    print("\nx  {:} {:>13}".format(name_str, '----'))
                        # Gaussian SPC file handling
                        if hasattr(bbe, "scf_energy") and not hasattr(bbe, "gibbs_free_energy"):
                            name_str = "x  " + name_str
                        # ORCA spc files
                        elif not hasattr(bbe, "scf_energy") and not hasattr(bbe, "gibbs_free_energy"):
                            name_str = "x  " + name_str
                        if hasattr(bbe, "scf_energy"):
                            e_str = '{:13.6f}'.format(bbe.scf_energy)
                        # No freqs found
                        if not hasattr(bbe, "gibbs_free_energy"):
                            warning("Could not find frequency information for: {}".format(os.path.relpath(file)))
                        else:
                            if not options.media:
                                if all(getattr(bbe, attrib) for attrib in ["enthalpy", "entropy", "qh_entropy",
                                                                           "gibbs_free_energy",
                                                                           "qh_gibbs_free_energy"]):
                                    if options.qh:
                                        print('{} {} {:10.6f} {:13.6f} {:13.6f} {:10.6f} {:10.6f} {:13.6f} {:13.6f}'.
                                              format(name_str, e_str, bbe.zpe, bbe.enthalpy, bbe.qh_enthalpy,
                                                     (options.temperature * bbe.entropy),
                                                     (options.temperature * bbe.qh_entropy), bbe.gibbs_free_energy,
                                                     bbe.qh_gibbs_free_energy))
                                    else:
                                        print('{} {} {:10.6f} {:13.6f} {:10.6f} {:10.6f} {:13.6f} {:13.6f}'.
                                              format(name_str, e_str, bbe.zpe, bbe.enthalpy,
                                                     (options.temperature * bbe.entropy),
                                                     (options.temperature * bbe.qh_entropy),
                                                     bbe.gibbs_free_energy, bbe.qh_gibbs_free_energy))
                            else:
                                # Media correction based on standard concentration of solvent
                                if options.media.lower() in SOLVENTS and options.media.lower() == \
                                        os.path.splitext(os.path.basename(file))[0].lower():
                                    mw_solvent = SOLVENTS[options.media.lower()][0]
                                    density_solvent = SOLVENTS[options.media.lower()][1]
                                    concentration_solvent = (density_solvent * 1000) / mw_solvent
                                    media_correction = -(GAS_CONSTANT / AU_TO_J) * np.log(concentration_solvent)
                                    if all(getattr(bbe, attrib) for attrib in ["enthalpy", "entropy", "qh_entropy",
                                                                               "gibbs_free_energy",
                                                                               "qh_gibbs_free_energy"]):
                                        if options.qh:
                                            print('{} {:10.6f} {:13.6f} {:13.6f} {:10.6f} {:10.6f} {:13.6f} {:13.6f}'
                                                  .format(name_str, bbe.zpe, bbe.enthalpy, bbe.qh_enthalpy,
                                                          (options.temperature * (bbe.entropy + media_correction)),
                                                          (options.temperature * (bbe.qh_entropy +
                                                                                  media_correction)),
                                                          bbe.gibbs_free_energy + (options.temperature *
                                                                                   (-media_correction)),
                                                          bbe.qh_gibbs_free_energy + (options.temperature *
                                                                                      (-media_correction))))
                                            print("  Solvent")
                                        else:
                                            print('{} {:10.6f} {:13.6f} {:10.6f} {:10.6f} {:13.6f} {:13.6f}'
                                                  ''.format(name_str, bbe.zpe, bbe.enthalpy,
                                                            (options.temperature * (bbe.entropy +
                                                                                    media_correction)),
                                                            (options.temperature * (
                                                                    bbe.qh_entropy + media_correction)),
                                                            bbe.gibbs_free_energy + (options.temperature * (
                                                                -media_correction)),
                                                            bbe.qh_gibbs_free_energy + (options.temperature * (
                                                                -media_correction))))
                                            print("  Solvent")
                                else:
                                    if all(getattr(bbe, attrib) for attrib in ["enthalpy", "entropy", "qh_entropy",
                                                                               "gibbs_free_energy",
                                                                               "qh_gibbs_free_energy"]):
                                        if options.qh:
                                            print('{} {:10.6f} {:13.6f} {:13.6f} {:10.6f} {:10.6f} {:13.6f} '
                                                  '{:13.6f}'.format(name_str, bbe.zpe, bbe.enthalpy, bbe.qh_enthalpy,
                                                                    (options.temperature * bbe.entropy),
                                                                    (options.temperature * bbe.qh_entropy),
                                                                    bbe.gibbs_free_energy,
                                                                    bbe.qh_gibbs_free_energy))
                                        else:
                                            print('{} {:10.6f} {:13.6f} {:10.6f} {:10.6f} {:13.6f} '
                                                  '{:13.6f}'.format(name_str, bbe.zpe, bbe.enthalpy,
                                                                    (options.temperature * bbe.entropy),
                                                                    (options.temperature * bbe.qh_entropy),
                                                                    bbe.gibbs_free_energy,
                                                                    bbe.qh_gibbs_free_energy))
                    # Append requested options to end of output
                    if options.cosmo and cosmo_solv is not None:
                        print('{:13.6f} {:16.6f}'.format(cosmo_solv[file], bbe.qh_gibbs_free_energy +
                                                         cosmo_solv[file]))
                    if options.boltz:
                        print('{:7.3f}'.format(boltz_facts[file] / boltz_sum))
                    if options.imag_freq and hasattr(bbe, "im_frequency_wn"):
                        for freq in bbe.im_frequency_wn:
                            print('{:9.2f}'.format(freq))
                    if options.ssymm:
                        if hasattr(bbe, "qh_gibbs_free_energy"):
                            print('{:>13}'.format(bbe.point_group))
                        else:
                            print('{:>37}'.format('---'))
                # Cluster files if requested
                if clustering:
                    dashes = "-" * (len(delimiter_row) - 3)
                    for n, cluster in enumerate(clusters):
                        for struct_id, structure in enumerate(cluster):
                            if structure == file:
                                if struct_id == len(cluster) - 1:
                                    print("\n   " + dashes)
                                    print("\n   " + '{name:<{var_width}} {gval:13.6f} {weight:6.2f}'.format(
                                        name='Boltzmann-weighted Cluster ' + ALPHABET[n].upper(),
                                        var_width=len(delimiter_row) - 24,
                                        gval=weighted_free_energy['cluster-' + ALPHABET[n].upper()] / boltz_facts[
                                            'cluster-' + ALPHABET[n].upper()],
                                        weight=100 * boltz_facts['cluster-' + ALPHABET[n].upper()] / boltz_sum))
                                    print("\n   " + dashes)
            print(delimiter_row)

        # Print CPU usage if requested
        if options.cputime:
            print('   {:<13} {:>2} {:>4} {:>2} {:>3} {:>2} {:>4} {:>2} {:>4}\n'
                  .format('TOTAL CPU', total_cpu_time.day + add_days - 1, 'days', total_cpu_time.hour, 'hrs',
                          total_cpu_time.minute, 'mins', total_cpu_time.second, 'secs'))

        # Tabulate relative values
        if options.pes:
            output_pes_data(options, thermo_data, delimiter_row, interval, interval_bbe_data, interval_thermo_data,
                            files)

        if options.ee:
            # Compute enantiomeric excess
            calc_enantio_excess(clustering, clusters, dup_list, files, options, thermo_data)

        if options.graph:
            # Graph reaction profiles
            create_plot(options, thermo_data)

    except (InvalidDataError, IOError) as e:
        warning(e)
        return INVALID_DATA

    return GOOD_RET  # success


def compute_thermochem(files, options, cosmo_solv=None, ssymm_option=False, vmm_option=False):
    bbe_vals = []
    for file in files:
        if options.cosmo:
            cosmo_option = cosmo_solv[file]
        else:
            cosmo_option = None

        # computes D3 term if requested, which is then sent to calc_bbe as a correction
        d3_energy = 0.0
        # The following is commented out because the called repo/code is no longer available
        # if options.D3 or options.D3BJ:
        #     verbose, intermolecular, pairwise, abc_term = False, False, False, False
        #     s6, rs6, s8, bj_a1, bj_a2 = 0.0, 0.0, 0.0, 0.0, 0.0
        #     functional = find_level_of_theory(file).split('/')[0]
        #     if options.D3:
        #         damp = 'zero'
        #     elif options.D3BJ:
        #         damp = 'bj'
        #     if options.ATM: abc_term = True
        #     try:
        #         file_data = GetOutData(file)
        #         d3_calc = dftd3.calcD3(file_data, functional, s6, rs6, s8, bj_a1, bj_a2, damp, abc_term,
        #                                intermolecular, pairwise, verbose)
        #         d3_energy = (d3_calc.attractive_r6_vdw + d3_calc.attractive_r8_vdw + d3_calc.repulsive_abc) / \
        #                     EHPART_TO_KCAL_MOL
        #     except:
        #         print('\n   ! Dispersion Correction Failed')
        #         d3_energy = 0.0
        bbe = CalcBBE(file, options.qs, options.qh, options.S_freq_cutoff, options.h_freq_cutoff,
                      options.temperature, options.conc, options.freq_scale_factor, options.zpe_scale_factor,
                      options.freespace, options.spc, options.invert, d3_energy=d3_energy,
                      cosmo=cosmo_option, ssymm=ssymm_option, mm_freq_scale_factor=vmm_option)

        # Populate bbe_vals with individual bbe entries for each file
        bbe_vals.append(bbe)
    # Creates a new dictionary object thermo_data, which attaches the bbe data to each file-name
    thermo_data = dict(zip(files, bbe_vals))  # The collected thermochemical data for all files
    return thermo_data


if __name__ == '__main__':
    status = main()
    sys.exit(status)
