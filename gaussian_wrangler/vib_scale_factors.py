#!/usr/bin/env python
"""
Comments and/or additions are welcome (send e-mail to: robert.paton@chem.ox.ac.uk

vib_scale_factors.py
Written by:  Rob Paton and Guilian Luchini
Last modified:  2019
Further modified by Heather Mayes, Dec 2019

Frequency scaling factors, taken from version 4 of The Truhlar group database
(https://t1.chem.umn.edu/freqscale/index.html
I. M. Alecu, J. Zheng, Y. Zhao, and D. G. Truhlar, J. Chem. Theory Comput. 6, 2872-2887 (2010).

The array is ordered as:
[level/basis set, zpe_fac, zpe_ref, zpe_meth, harm_fac, harm_ref, harm_meth, fund_fac, fund_ref, fund_meth]
where zpe_fac, harm_fac and fund_fac are the scaling factors for ZPEs, harmonic frequencies, and fundamentals,
respectively.

The ref and meth elements refer to original references and method of determination by the Truhlar group.
All information taken from https://comp.chem.umn.edu/freqscale/190107_Database_of_Freq_Scale_Factors_v4.pdf

Methods
D: The scale factor was directly obtained from the ZPVE15/10 or F38/10 databases given in Ref. 1.
C: The scale factor was obtained by applying a small systematic correction of -0.0025 to preexisting scale factor.
   The references for the preexisting (uncorrected) scale factors are given in Supporting Information of Ref. 1 and in
   Version 1 of this database
R: The scale factor was obtained via the Reduced Scale Factor Optimization Model described in Ref. 1. Briefly, this
   entails using the ZPE6 database for determining ZPE scale factors, and/or using the universal scale factor ratios
   of aF/ZPE = 0.974 and aH/ZPE = 1.014 to obtain the respective values for the scale factors for fundamental and
   harmonic frequencies.
"""

import os
import sys
import numpy as np
from common_wrangler.common import (InvalidDataError,
                                    SPEED_OF_LIGHT, GAS_CONSTANT, KB, PLANCK_CONST_JS, AVOGADRO_CONST, AMU_TO_KG,
                                    AU_TO_J)

# # To make the output exactly match Gaussian's, use the values below instead importing them from common_wrangler.common
# SPEED_OF_LIGHT = 2.99792458e10  # cm / s (same as in common)
# GAS_CONSTANT = 8.31441  # J / K / mol; in common, GAS_CONSTANT = 8.314462618 (in Goodvibes: 8.3144621)
# KB = 1.380662e-23  # Boltzmann's constant in J/K; in common, the value is 1.380649e-23 (in Goodvibes: 1.3806488e-23)
# H = 6.626176e-34  # Planck's constant in J/s; in common, the value is 6.62607015e-34 (in Goodvibes: 6.62606957e-34)
# AVOGADRO_CONST = 6.0221415e23  # 1 / mol; in common, the value is 6.02214076e23 (Goodvibes == Gaussian)
# AMU_TO_KG = 1.66053886E-27  # UNIT CONVERSION per GoodVibes; 1.66053906660e-27 in common
# EHPART_TO_KCAL_MOL = 627.5095  # kcal/mol/(Eh/part); in common, the value is 627.5094709  (in Goodvibes: 627.509541)
# AU_TO_J = 4.184 * EHPART_TO_KCAL_MOL * 1000.0  # This value changes based on which EHPART_TO_KCAL_MOL is used

Str_char = "U%d"

SCALING_REFS = ["none",
                "I. M. Alecu, J. Zheng, Y. Zhao, and D. G. Truhlar, J. Chem. Theory Comput. 6, 2872-2887 (2010).",
                "Y. Zhao and D. G. Truhlar, unpublished (2003), modified by systematic correction of -0.0025 by I. M. "
                "Alecu (2010).",
                "I. M. Alecu, unpublished (2011).",
                "J. Zheng, R. J. Rocha, M. Pelegrini, L. F. A. Ferrao, E. F. V. Carvalho, O. Roberto-Neto, F. B. C. "
                "Machado, and D. G. Truhlar, J. Chem. Phys. 136, 184310/1-10 (2012).",
                "J. Zheng and D. G. Truhlar, unpublished (2014).", "J. Bao and D. G. Truhlar, unpublished (2014).",
                "H. Yu, J. Zheng, and D. G. Truhlar, unpublished (2015)",
                "S. Kanchanakungwankul, J. L. Bao, J. Zheng, I. M. Alecu, B. J. Lynch, Y. Zhao, and D. G. Truhlar, "
                "unpublished (2018)"]

SCALING_DATA = np.array([('AM1', 0.948, 1, 'R', 0.961, 1, 'R', 0.923, 1, 'R'),
                         ('B1B95/6-31+G(d,p)', 0.971, 1, 'C', 0.985, 1, 'R', 0.946, 1, 'R'),
                         ('B1B95/MG3S', 0.973, 1, 'C', 0.987, 1, 'R', 0.948, 1, 'R'),
                         ('B1LYP/MG3S', 0.978, 1, 'D', 0.994, 1, 'D', 0.955, 1, 'D'),
                         ('B3LYP/6-31G(2df,2p)', 0.981, 1, 'C', 0.995, 1, 'R', 0.955, 1, 'R'),
                         ('B3LYP/6-31G(d)', 0.977, 1, 'R', 0.991, 1, 'R', 0.952, 1, 'R'),
                         ('B3LYP/aug-cc-pVTZ', 0.985, 3, 'R', 0.999, 3, 'R', 0.959, 3, 'R'),
                         ('B3LYP/def2TZVP', 0.985, 3, 'R', 0.999, 3, 'R', 0.959, 3, 'R'),
                         ('B3LYP/ma-TZVP', 0.986, 1, 'R', 1.0, 1, 'R', 0.96, 1, 'R'),
                         ('B3LYP/MG3S', 0.983, 1, 'D', 0.998, 1, 'D', 0.96, 1, 'D'),
                         ('B3P86/6-31G(d)', 0.971, 1, 'R', 0.985, 1, 'R', 0.946, 1, 'R'),
                         ('B3PW91/6-31G(d)', 0.972, 1, 'R', 0.986, 1, 'R', 0.947, 1, 'R'),
                         ('B973/def2TZVP', 0.974, 8, 'D', 0.988, 8, 'R', 0.949, 8, 'R'),
                         ('B973/ma-TZVP', 0.975, 1, 'R', 0.989, 1, 'R', 0.95, 1, 'R'),
                         ('B973/MG3S', 0.972, 1, 'D', 0.986, 1, 'D', 0.947, 1, 'D'),
                         ('B98/def2TZVP', 0.984, 1, 'R', 0.998, 1, 'R', 0.958, 1, 'R'),
                         ('B98/ma-TZVP', 0.985, 1, 'R', 0.999, 1, 'R', 0.959, 1, 'R'),
                         ('B98/MG3S', 0.982, 1, 'D', 0.995, 1, 'D', 0.956, 1, 'D'),
                         ('BB1K/6-31+G(d,p)', 0.954, 1, 'C', 0.967, 1, 'R', 0.929, 1, 'R'),
                         ('BB1K/MG3S', 0.957, 1, 'C', 0.97, 1, 'R', 0.932, 1, 'R'),
                         ('BB95/6-31+G(d,p)', 1.011, 1, 'C', 1.025, 1, 'R', 0.985, 1, 'R'),
                         ('BB95/MG3S', 1.012, 1, 'C', 1.026, 1, 'R', 0.986, 1, 'R'),
                         ('BLYP/6-311G(df,p)', 1.013, 1, 'R', 1.027, 1, 'R', 0.987, 1, 'R'),
                         ('BLYP/6-31G(d)', 1.009, 1, 'R', 1.023, 1, 'R', 0.983, 1, 'R'),
                         ('BLYP/MG3S', 1.013, 1, 'D', 1.031, 1, 'D', 0.991, 1, 'D'),
                         ('BMC-CCSD', 0.985, 1, 'D', 1.001, 1, 'D', 0.962, 1, 'D'),
                         ('BMK/ma-TZVP', 0.972, 1, 'R', 0.986, 1, 'R', 0.947, 1, 'R'),
                         ('BMK/MG3S', 0.971, 1, 'D', 0.984, 1, 'D', 0.945, 1, 'D'),
                         ('BP86/6-31G(d)', 1.007, 1, 'R', 1.021, 1, 'R', 0.981, 1, 'R'),
                         ('BP86/ma-TZVP', 1.014, 1, 'R', 1.028, 1, 'R', 0.988, 1, 'R'),
                         ('BPW60/6-311+G(d,p)', 0.934, 2, 'C', 0.91, 2, 'R', 0.947, 2, 'R'),
                         ('BPW63/MG3S', 0.923, 2, 'C', 0.899, 2, 'R', 0.936, 2, 'R'),
                         ('CAM-B3LYP/ma-TZVP', 0.976, 1, 'R', 0.99, 1, 'R', 0.951, 1, 'R'),
                         ('CCSD(T)/jul-cc-pVTZ', 0.984, 1, 'R', 0.998, 1, 'R', 0.958, 1, 'R'),
                         ('CCSD(T)/aug-cc-pVTZ', 0.987, 1, 'R', 1.001, 1, 'R', 0.961, 1, 'R'),
                         ('CCSD(T)-F12/jul-cc-pVTZ', 0.981, 1, 'R', 0.995, 1, 'R', 0.955, 1, 'R'),
                         ('CCSD(T)-F12a/cc-pVDZ-F12', 0.983, 11, 'R', 0.997, 11, 'R', 0.957, 11, 'R'),
                         ('CCSD(T)-F12a/cc-pVTZ-F12', 0.984, 1, 'R', 0.998, 1, 'R', 0.958, 1, 'R'),
                         ('CCSD(T)-F12b/VQZF12//CCSD(T)-F12a/TZF', 0.984, 13, 'R', 0.998, 13, 'R', 0.958, 13, 'R'),
                         ('CCSD(T)-F12b/VQZF12//CCSD(T)-F12a/DZF', 0.983, 13, 'R', 0.997, 13, 'R', 0.957, 13, 'R'),
                         ('CCSD/jul-cc-pVTZ', 0.973, 1, 'R', 0.987, 1, 'R', 0.948, 1, 'R'),
                         ('CCSD-F12/jul-cc-pVTZ', 0.971, 1, 'R', 0.985, 1, 'R', 0.946, 1, 'R'),
                         ('G96LYP80/6-311+G(d,p)', 0.911, 2, 'C', 0.887, 2, 'R', 0.924, 2, 'R'),
                         ('G96LYP82/MG3S', 0.907, 2, 'C', 0.883, 2, 'R', 0.92, 2, 'R'),
                         ('GAM/def2TZVP', 0.98, 7, 'D', 0.994, 7, 'D', 0.955, 7, 'D'),
                         ('GAM/ma-TZVP', 0.981, 7, 'D', 0.995, 7, 'D', 0.956, 7, 'D'),
                         ('HF/3-21G', 0.919, 1, 'R', 0.932, 1, 'R', 0.895, 1, 'R'),
                         ('HF/6-31+G(d)', 0.911, 1, 'R', 0.924, 1, 'R', 0.887, 1, 'R'),
                         ('HF/6-31+G(d,p)', 0.915, 1, 'C', 0.928, 1, 'R', 0.891, 1, 'R'),
                         ('HF/6-311G(d,p)', 0.92, 1, 'R', 0.933, 1, 'R', 0.896, 1, 'R'),
                         ('HF/6-311G(df,p)', 0.92, 1, 'R', 0.933, 1, 'R', 0.896, 1, 'R'),
                         ('HF/6-31G(d)', 0.909, 1, 'R', 0.922, 1, 'R', 0.885, 1, 'R'),
                         ('HF/6-31G(d,p)', 0.913, 1, 'R', 0.926, 1, 'R', 0.889, 1, 'R'),
                         ('HF/MG3S', 0.919, 1, 'D', 0.932, 1, 'D', 0.895, 1, 'D'),
                         ('HFLYP/MG3S', 0.899, 1, 'D', 0.912, 1, 'D', 0.876, 1, 'D'),
                         ('HSEh1PBE/ma-TZVP', 0.979, 1, 'R', 0.993, 1, 'R', 0.954, 1, 'R'),
                         ('M05/aug-cc-pVTZ', 0.978, 1, 'R', 0.992, 1, 'R', 0.953, 1, 'R'),
                         ('M05/def2TZVP', 0.978, 3, 'R', 0.991, 3, 'R', 0.952, 3, 'R'),
                         ('M05/ma-TZVP', 0.979, 1, 'R', 0.993, 1, 'R', 0.954, 1, 'R'),
                         ('M05/maug-cc-pVTZ', 0.978, 1, 'R', 0.992, 1, 'R', 0.953, 1, 'R'),
                         ('M05/MG3S', 0.977, 1, 'D', 0.989, 1, 'D', 0.951, 1, 'D'),
                         ('M052X/6-31+G(d,p)', 0.961, 1, 'D', 0.974, 1, 'D', 0.936, 1, 'D'),
                         ('M052X/aug-cc-pVTZ', 0.964, 1, 'R', 0.977, 1, 'R', 0.939, 1, 'R'),
                         ('M052X/def2TZVPP', 0.962, 1, 'D', 0.976, 1, 'D', 0.938, 1, 'D'),
                         ('M052X/ma-TZVP', 0.965, 1, 'R', 0.979, 1, 'R', 0.94, 1, 'R'),
                         ('M052X/maug-cc-pVTZ', 0.964, 1, 'R', 0.977, 1, 'R', 0.939, 1, 'R'),
                         ('M052X/MG3S', 0.962, 1, 'D', 0.975, 1, 'D', 0.937, 1, 'D'),
                         ('M06/6-31+G(d,p)', 0.98, 1, 'D', 0.989, 1, 'D', 0.95, 1, 'D'),
                         ('M06/aug-cc-pVTZ', 0.984, 1, 'R', 0.998, 1, 'R', 0.958, 1, 'R'),
                         ('M06/def2TZVP', 0.982, 3, 'R', 0.996, 3, 'R', 0.956, 3, 'R'),
                         ('M06/def2TZVPP', 0.979, 1, 'D', 0.992, 1, 'D', 0.953, 1, 'D'),
                         ('M06/ma-TZVP', 0.982, 1, 'R', 0.996, 1, 'R', 0.956, 1, 'R'),
                         ('M06/maug-cc-pVTZ', 0.982, 1, 'R', 0.996, 1, 'R', 0.956, 1, 'R'),
                         ('M06/MG3S', 0.981, 1, 'D', 0.994, 1, 'D', 0.955, 1, 'D'),
                         ('M062X/6-31+G(d,p)', 0.967, 1, 'D', 0.979, 1, 'D', 0.94, 1, 'D'),
                         ('M062X/6-311+G(d,p)', 0.97, 5, 'D', 0.983, 5, 'R', 0.944, 5, 'R'),
                         ('M062X/6-311++G(d,p)', 0.97, 5, 'D', 0.983, 5, 'R', 0.944, 5, 'R'),
                         ('M062X/aug-cc-pVDZ', 0.979, 14, 'D', 0.993, 14, 'R', 0.954, 14, 'R'),
                         ('M062X/aug-cc-pVTZ', 0.971, 1, 'D', 0.985, 1, 'D', 0.946, 1, 'D'),
                         ('M062X/def2TZVP', 0.971, 7, 'D', 0.984, 7, 'D', 0.946, 7, 'D'),
                         ('M062X/def2QZVP', 0.97, 7, 'D', 0.983, 7, 'D', 0.945, 7, 'D'),
                         ('M062X/def2TZVPP', 0.97, 1, 'D', 0.983, 1, 'D', 0.945, 1, 'D'),
                         ('M062X/jul-cc-pVDZ', 0.977, 14, 'D', 0.991, 14, 'R', 0.952, 14, 'R'),
                         ('M062X/jul-cc-pVTZ', 0.971, 14, 'D', 0.985, 14, 'R', 0.946, 14, 'R'),
                         ('M062X/jun-cc-pVDZ', 0.976, 14, 'D', 0.99, 14, 'R', 0.951, 14, 'R'),
                         ('M062X/jun-cc-pVTZ', 0.971, 14, 'D', 0.985, 14, 'R', 0.946, 14, 'R'),
                         ('M062X/ma-TZVP', 0.972, 1, 'R', 0.986, 1, 'R', 0.947, 1, 'R'),
                         ('M062X/maug-cc-pV(T+d)Z', 0.971, 1, 'D', 0.984, 1, 'D', 0.945, 1, 'D'),
                         ('M062X/MG3S', 0.97, 1, 'D', 0.982, 1, 'D', 0.944, 1, 'D'),
                         ('M06HF/6-31+G(d,p)', 0.954, 1, 'D', 0.969, 1, 'D', 0.931, 1, 'D'),
                         ('M06HF/aug-cc-pVTZ', 0.961, 1, 'R', 0.974, 1, 'R', 0.936, 1, 'R'),
                         ('M06HF/def2TZVPP', 0.958, 1, 'D', 0.97, 1, 'D', 0.932, 1, 'D'),
                         ('M06HF/ma-TZVP', 0.957, 1, 'R', 0.97, 1, 'R', 0.932, 1, 'R'),
                         ('M06HF/maug-cc-pVTZ', 0.959, 1, 'R', 0.972, 1, 'R', 0.934, 1, 'R'),
                         ('M06HF/MG3S', 0.955, 1, 'D', 0.967, 1, 'D', 0.93, 1, 'D'),
                         ('M06L/6-31G(d,p)', 0.977, 15, 'D', 0.991, 15, 'R', 0.952, 15, 'R'),
                         ('M06L/6-31+G(d,p)', 0.978, 1, 'D', 0.992, 1, 'D', 0.953, 1, 'D'),
                         ('M06L/aug-cc-pVTZ', 0.98, 1, 'R', 0.994, 1, 'R', 0.955, 1, 'R'),
                         ('M06L/aug-cc-pV(T+d)Z', 0.98, 9, 'R', 0.994, 9, 'R', 0.955, 9, 'R'),
                         ('M06L/aug-cc-pVTZ-pp', 0.98, 9, 'R', 0.994, 9, 'R', 0.955, 9, 'R'),
                         ('M06L(DKH2)/aug-cc-pwcVTZ-DK', 0.985, 1, 'D', 0.999, 1, 'R', 0.959, 1, 'R'),
                         ('M06L/def2TZVP', 0.976, 3, 'R', 0.99, 3, 'R', 0.951, 3, 'R'),
                         ('M06L/def2TZVPP', 0.976, 1, 'D', 0.995, 1, 'D', 0.956, 1, 'D'),
                         ('M06L/ma-TZVP', 0.977, 1, 'R', 0.991, 1, 'R', 0.952, 1, 'R'),
                         ('M06L/maug-cc-pVTZ', 0.977, 1, 'R', 0.991, 1, 'R', 0.952, 1, 'R'),
                         ('M06L/MG3S', 0.978, 1, 'D', 0.996, 1, 'D', 0.958, 1, 'D'),
                         ('M08HX/6-31+G(d,p)', 0.972, 1, 'D', 0.983, 1, 'D', 0.944, 1, 'D'),
                         ('M08HX/aug-cc-pVTZ', 0.975, 1, 'R', 0.989, 1, 'R', 0.95, 1, 'R'),
                         ('M08HX/cc-pVTZ+', 0.974, 1, 'D', 0.985, 1, 'D', 0.946, 1, 'D'),
                         ('M08HX/def2TZVPP', 0.973, 1, 'D', 0.984, 1, 'D', 0.945, 1, 'D'),
                         ('M08HX/jun-cc-pVTZ', 0.974, 6, 'D', 0.986, 6, 'D', 0.947, 6, 'D'),
                         ('M08HX/ma-TZVP', 0.976, 1, 'R', 0.99, 1, 'R', 0.951, 1, 'R'),
                         ('M08HX/maug-cc-pVTZ', 0.976, 1, 'R', 0.99, 1, 'R', 0.951, 1, 'R'),
                         ('M08HX/MG3S', 0.973, 1, 'D', 0.984, 1, 'D', 0.946, 1, 'D'),
                         ('M08SO/6-31+G(d,p)', 0.979, 1, 'D', 0.989, 1, 'D', 0.951, 1, 'D'),
                         ('M08SO/aug-cc-pVTZ', 0.985, 1, 'R', 0.999, 1, 'R', 0.959, 1, 'R'),
                         ('M08SO/cc-pVTZ+', 0.982, 1, 'D', 0.995, 1, 'D', 0.956, 1, 'D'),
                         ('M08SO/def2TZVPP', 0.98, 1, 'D', 0.993, 1, 'D', 0.954, 1, 'D'),
                         ('M08SO/ma-TZVP', 0.984, 1, 'R', 0.998, 1, 'R', 0.958, 1, 'R'),
                         ('M08SO/maug-cc-pVTZ', 0.983, 1, 'R', 0.997, 1, 'R', 0.957, 1, 'R'),
                         ('M08SO/MG3', 0.984, 4, 'D', 0.998, 4, 'R', 0.959, 4, 'R'),
                         ('M08SO/MG3S', 0.983, 1, 'D', 0.995, 1, 'D', 0.956, 1, 'D'),
                         ('M08SO/MG3SXP', 0.984, 1, 'D', 0.996, 1, 'D', 0.957, 1, 'D'),
                         ('M11L/maug-cc-pVTZ', 0.988, 16, 'D', 1.002, 16, 'R', 0.962, 16, 'R'),
                         ('MN11-L/MG3S', 0.985, 16, 'D', 0.999, 16, 'R', 0.959, 16, 'R'),
                         ('MN12L/jul-cc-pVDZ', 0.974, 14, 'R', 0.988, 14, 'R', 0.95, 14, 'R'),
                         ('MN12L/MG3S', 0.968, 6, 'D', 0.981, 6, 'D', 0.943, 6, 'D'),
                         ('MN12SX/6-311++G(d,p)', 0.976, 6, 'D', 0.986, 6, 'D', 0.947, 6, 'D'),
                         ('MN12SX/jul-cc-pVDZ', 0.979, 14, 'R', 0.993, 14, 'R', 0.954, 14, 'R'),
                         ('MN15L/MG3S', 0.977, 1, 'D', 0.991, 1, 'R', 0.952, 1, 'R'),
                         ('MN15L/maug-cc-pVTZ', 0.979, 1, 'D', 0.993, 1, 'R', 0.954, 1, 'R'),
                         ('MC3BB', 0.965, 1, 'C', 0.979, 1, 'R', 0.94, 1, 'R'),
                         ('MC3MPW', 0.964, 1, 'C', 0.977, 1, 'R', 0.939, 1, 'R'),
                         ('MC-QCISD/3', 0.992, 1, 'C', 1.006, 1, 'R', 0.966, 1, 'R'),
                         ('MOHLYP/ma-TZVP', 1.027, 1, 'R', 1.041, 1, 'R', 1.0, 1, 'R'),
                         ('MOHLYP/MG3S', 1.022, 1, 'R', 1.036, 1, 'R', 0.995, 1, 'R'),
                         ('MP2(FC)/6-31+G(d,p)', 0.968, 1, 'C', 0.982, 1, 'R', 0.943, 1, 'R'),
                         ('MP2(FC)/6-311G(d,p)', 0.97, 1, 'R', 0.984, 1, 'R', 0.945, 1, 'R'),
                         ('MP2(FC)/6-31G(d)', 0.964, 1, 'R', 0.977, 1, 'R', 0.939, 1, 'R'),
                         ('MP2(FC)/6-31G(d,p)', 0.958, 1, 'R', 0.971, 1, 'R', 0.933, 1, 'R'),
                         ('MP2(FC)/cc-pVDZ', 0.977, 1, 'C', 0.991, 1, 'R', 0.952, 1, 'R'),
                         ('MP2(FC)/cc-pVTZ', 0.975, 1, 'D', 0.992, 1, 'D', 0.953, 1, 'D'),
                         ('MP2(FULL)/6-31G(d)', 0.963, 1, 'R', 0.976, 1, 'R', 0.938, 1, 'R'),
                         ('MP4(SDQ)/jul-cc-pVTZ', 0.973, 1, 'R', 0.987, 1, 'R', 0.948, 1, 'R'),
                         ('MPW1B95/6-31+G(d,p)', 0.97, 1, 'C', 0.984, 1, 'R', 0.945, 1, 'R'),
                         ('MPW1B95/MG3', 0.97, 1, 'C', 0.984, 1, 'R', 0.945, 1, 'R'),
                         ('MPW1B95/MG3S', 0.972, 1, 'C', 0.986, 1, 'R', 0.947, 1, 'R'),
                         ('MPW1K/6-31+G(d,p)', 0.949, 1, 'C', 0.962, 1, 'R', 0.924, 1, 'R'),
                         ('MPW1K/aug-cc-PDTZ', 0.959, 14, 'R', 0.972, 14, 'R', 0.934, 14, 'R'),
                         ('MPW1K/aug-cc-PVTZ', 0.955, 14, 'R', 0.968, 14, 'R', 0.93, 14, 'R'),
                         ('MPW1K/jul-cc-pVDZ', 0.957, 14, 'R', 0.97, 14, 'R', 0.932, 14, 'R'),
                         ('MPW1K/jul-cc-pVTZ', 0.954, 14, 'R', 0.967, 14, 'R', 0.929, 14, 'R'),
                         ('MPW1K/jun-cc-pVDZ', 0.955, 14, 'R', 0.968, 14, 'R', 0.93, 14, 'R'),
                         ('MPW1K/jun-cc-pVTZ', 0.954, 14, 'R', 0.967, 14, 'R', 0.929, 14, 'R'),
                         ('MPW1K/ma-TZVP', 0.956, 1, 'R', 0.969, 1, 'R', 0.931, 1, 'R'),
                         ('MPW1K/MG3', 0.953, 1, 'C', 0.966, 1, 'R', 0.928, 1, 'R'),
                         ('MPW1K/MG3S', 0.956, 1, 'C', 0.969, 1, 'R', 0.931, 1, 'R'),
                         ('MPW1K/MIDI!', 0.953, 1, 'R', 0.966, 1, 'R', 0.928, 1, 'R'),
                         ('MPW1K/MIDIY', 0.947, 1, 'R', 0.96, 1, 'R', 0.922, 1, 'R'),
                         ('MPW3LYP/6-31+G(d,p)', 0.98, 1, 'C', 0.994, 1, 'R', 0.955, 1, 'R'),
                         ('MPW3LYP/6-311+G(2d,p)', 0.986, 1, 'R', 1.0, 1, 'R', 0.96, 1, 'R'),
                         ('MPW3LYP/6-31G(d)', 0.976, 1, 'R', 0.99, 1, 'R', 0.951, 1, 'R'),
                         ('MPW3LYP/ma-TZVP', 0.986, 1, 'R', 1.0, 1, 'R', 0.96, 1, 'R'),
                         ('MPW3LYP/MG3S', 0.982, 1, 'C', 0.996, 1, 'R', 0.956, 1, 'R'),
                         ('MPW74/6-311+G(d,p)', 0.912, 2, 'C', 0.888, 2, 'R', 0.925, 2, 'R'),
                         ('MPW76/MG3S', 0.909, 2, 'C', 0.885, 2, 'R', 0.922, 2, 'R'),
                         ('MPWB1K/6-31+G(d,p)', 0.951, 1, 'C', 0.964, 1, 'R', 0.926, 1, 'R'),
                         ('MPWB1K/MG3S', 0.954, 1, 'C', 0.967, 1, 'R', 0.929, 1, 'R'),
                         ('MPWLYP1M/ma-TZVP', 1.009, 1, 'R', 1.023, 1, 'R', 0.983, 1, 'R'),
                         ('MW3.2//CCSD(T)-F12a/TZF', 0.984, 13, 'R', 0.998, 13, 'R', 0.958, 13, 'R'),
                         ('OreLYP/ma-TZVP', 1.01, 7, 'D', 1.024, 7, 'D', 0.984, 7, 'D'),
                         ('OreLYP/def2TZVP', 1.008, 7, 'D', 1.023, 7, 'D', 0.982, 7, 'D'),
                         ('PBE/def2TZVP', 1.011, 3, 'R', 1.026, 3, 'R', 0.985, 3, 'R'),
                         ('PBE/MG3S', 1.01, 1, 'D', 1.025, 1, 'D', 0.985, 1, 'D'),
                         ('PBE/ma-TZVP', 1.014, 1, 'D', 1.028, 1, 'D', 0.987, 1, 'D'),
                         ('PBE0/MG3S', 0.975, 1, 'D', 0.989, 1, 'D', 0.95, 1, 'D'),
                         ('PBE1KCIS/MG3', 0.981, 1, 'C', 0.995, 1, 'R', 0.955, 1, 'R'),
                         ('PBE1KCIS/MG3S', 0.981, 1, 'C', 0.995, 1, 'R', 0.955, 1, 'R'),
                         ('PM3', 0.94, 1, 'R', 0.953, 1, 'R', 0.916, 1, 'R'),
                         ('PM6', 1.078, 1, 'R', 1.093, 1, 'R', 1.05, 1, 'R'),
                         ('PW6B95/def2TZVP', 0.974, 8, 'R', 0.988, 8, 'R', 0.949, 8, 'R'),
                         ('PWB6K/cc-pVDZ', 0.953, 12, 'D', 0.966, 12, 'R', 0.928, 12, 'R'),
                         ('QCISD/cc-pVTZ', 0.975, 11, 'R', 0.989, 11, 'R', 0.95, 11, 'R'),
                         ('QCISD/MG3S', 0.978, 10, 'R', 0.992, 10, 'R', 0.953, 10, 'R'),
                         ('QCISD(FC)/6-31G(d)', 0.973, 1, 'R', 0.987, 1, 'R', 0.948, 1, 'R'),
                         ('QCISD(T)/aug-cc-pVQZ', 0.989, 10, 'R', 1.003, 10, 'R', 0.963, 10, 'R'),
                         ('revTPSS/def2TZVP', 0.998, 7, 'D', 1.012, 7, 'D', 0.972, 7, 'D'),
                         ('revTPSS/ma-TZVP', 0.999, 7, 'D', 1.013, 7, 'D', 0.973, 7, 'D'),
                         ('SOGGA/ma-TZVP', 1.017, 1, 'R', 1.031, 1, 'R', 0.991, 1, 'R'),
                         ('THCTHhyb/ma-TZVP', 0.989, 1, 'R', 1.003, 1, 'R', 0.963, 1, 'R'),
                         ('TPSS1KCIS/def2TZVP', 0.982, 1, 'R', 0.996, 1, 'R', 0.956, 1, 'R'),
                         ('TPSS1KCIS/ma-TZVP', 0.983, 1, 'R', 0.997, 1, 'R', 0.957, 1, 'R'),
                         ('TPSSh/MG3S', 0.984, 1, 'D', 1.002, 1, 'D', 0.963, 1, 'D'),
                         ('VSXC/MG3S', 0.986, 1, 'D', 1.001, 1, 'D', 0.962, 1, 'D'),
                         ('wB97/def2TZVP', 0.969, 1, 'R', 0.983, 1, 'R', 0.944, 1, 'R'),
                         ('wB97/ma-TZVP', 0.97, 1, 'R', 0.984, 1, 'R', 0.945, 1, 'R'),
                         ('wB97X/def2TZVP', 0.97, 1, 'R', 0.984, 1, 'R', 0.945, 1, 'R'),
                         ('wB97X/ma-TZVP', 0.971, 1, 'R', 0.985, 1, 'R', 0.946, 1, 'R'),
                         ('wB97XD/def2TZVP', 0.975, 1, 'R', 0.989, 1, 'R', 0.95, 1, 'R'),
                         ('wB97XD/ma-TZVP', 0.975, 1, 'R', 0.989, 1, 'R', 0.95, 1, 'R'),
                         ('wB97XD/maug-cc-pVTZ', 0.974, 1, 'R', 0.988, 1, 'R', 0.949, 1, 'R'),
                         ('W3X//CCSD(T)-F12a/TZF', 0.984, 13, 'R', 0.998, 13, 'R', 0.958, 13, 'R'),
                         ('W3XL//CCSD(T)-F12a/TZF', 0.984, 13, 'R', 0.998, 13, 'R', 0.958, 13, 'R'),
                         ('W3XL//QCISD/STZ', 0.973, 13, 'R', 0.987, 13, 'R', 0.948, 13, 'R'),
                         ('X1B95/6-31+G(d,p)', 0.968, 1, 'C', 0.982, 1, 'R', 0.943, 1, 'R'),
                         ('X1B95/MG3S', 0.971, 1, 'C', 0.985, 1, 'R', 0.946, 1, 'R'),
                         ('XB1K/6-31+G(d,p)', 0.952, 1, 'C', 0.965, 1, 'R', 0.927, 1, 'R'),
                         ('XB1K/MG3S', 0.955, 1, 'C', 0.968, 1, 'R', 0.93, 1, 'R')],
                        dtype=[('level', Str_char % 40), ('zpe_fac', 'f4'), ('zpe_ref', 'i4'),
                               ('zpe_meth', Str_char % 1), ('harm_fac', 'f4'), ('harm_ref', 'i4'),
                               ('harm_meth', Str_char % 1), ('fund_fac', 'f4'), ('fund_ref', 'i4'),
                               ('fund_meth', Str_char % 1)])

PERIODIC_TABLE = ["", "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si",
                  "P", "S", "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
                  "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd",
                  "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm",
                  "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt",
                  "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu",
                  "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds",
                  "Rg", "Uub", "Uut", "Uuq", "Uup", "Uuh", "Uus", "Uuo"]


# Radii used to determine connectivity in symmetry corrections
# Covalent radii taken from Cambridge Structural Database
RADII = {'H': 0.32, 'He': 0.93, 'Li': 1.23, 'Be': 0.90, 'B': 0.82, 'C': 0.77, 'N': 0.75, 'O': 0.73, 'F': 0.72,
         'Ne': 0.71, 'Na': 1.54, 'Mg': 1.36, 'Al': 1.18, 'Si': 1.11, 'P': 1.06, 'S': 1.02, 'Cl': 0.99, 'Ar': 0.98,
         'K': 2.03, 'Ca': 1.74, 'Sc': 1.44, 'Ti': 1.32, 'V': 1.22, 'Cr': 1.18, 'Mn': 1.17, 'Fe': 1.17, 'Co': 1.16,
         'Ni': 1.15, 'Cu': 1.17, 'Zn': 1.25, 'Ga': 1.26, 'Ge': 1.22, 'As': 1.20, 'Se': 1.16, 'Br': 1.14, 'Kr': 1.12,
         'Rb': 2.16, 'Sr': 1.91, 'Y': 1.62, 'Zr': 1.45, 'Nb': 1.34, 'Mo': 1.30, 'Tc': 1.27, 'Ru': 1.25, 'Rh': 1.25,
         'Pd': 1.28, 'Ag': 1.34, 'Cd': 1.48, 'In': 1.44, 'Sn': 1.41, 'Sb': 1.40, 'Te': 1.36, 'I': 1.33, 'Xe': 1.31,
         'Cs': 2.35, 'Ba': 1.98, 'La': 1.69, 'Lu': 1.60, 'Hf': 1.44, 'Ta': 1.34, 'W': 1.30, 'Re': 1.28, 'Os': 1.26,
         'Ir': 1.27, 'Pt': 1.30, 'Au': 1.34, 'Hg': 1.49, 'Tl': 1.48, 'Pb': 1.47, 'Bi': 1.46, 'X': 0}
# Bondi van der Waals radii for all atoms from: Bondi, A. J. Phys. Chem. 1964, 68, 441-452,
# except hydrogen, which is taken from Rowland, R. S.; Taylor, R. J. Phys. Chem. 1996, 100, 7384-7391.
# Radii unavailable in either of these publications are set to 2 Angstrom
# (Unfinished)
BONDI = {'H': 1.09, 'He': 1.40, 'Li': 1.82, 'Be': 2.00, 'B': 2.00, 'C': 1.70, 'N': 1.55, 'O': 1.52, 'F': 1.47,
         'Ne': 1.54}


def element_id(mass_num, num=False):
    try:
        if num:
            return PERIODIC_TABLE.index(mass_num)
        return PERIODIC_TABLE[mass_num]
    except IndexError:
        return "XX"


class CalcBBE:
    # The function to compute the "black box" entropy and enthalpy values
    # along with all other thermochemical quantities
    def __init__(self, file, qs, qh, s_freq_cutoff, h_freq_cutoff, temperature, conc, freq_scale_factor,
                 zpe_scale_factor, solv='none', spc=False, invert=False,
                 d3_energy=0.0, ssymm=False, cosmo=None, mm_freq_scale_factor=False):

        h_damp, u_vib_qrrho, qh_u_vib = None, None, None   # make IDE happy

        # List of frequencies and default values
        im_freq_cutoff, frequency_wn, im_frequency_wn, rot_temp, roconst, linear_mol, link, freq_loc, linkmax, sym_no, \
            self.cpu, inverted_freqs = 0.0, [], [], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], 0, 0, 0, 0, 1, [0, 0, 0, 0, 0], []
        s_vib_rrqho = []  # make IDE happy
        linear_warning = False
        if mm_freq_scale_factor:
            fract_model_sys = []
            freq_scale_factor = [freq_scale_factor, mm_freq_scale_factor]
        else:
            fract_model_sys = False
        self.xyz = GetOutData(file)
        self.job_type = job_type(file)
        # Parse some useful information from the file
        self.sp_energy, self.program, self.version_program, self.solvation_model, self.file, self.charge, \
            self.empirical_dispersion, self.multiplicity = parse_data(file)
        with open(file) as f:
            g_output = f.readlines()
        self.cosmo_qhg = 0.0
        # Read any single point energies if requested
        if spc and spc != 'link':
            name, ext = os.path.splitext(file)
            try:
                self.sp_energy, self.sp_program, self.sp_version_program, self.sp_solvation_model, self.sp_file, \
                    self.sp_charge, self.sp_empirical_dispersion, self.sp_multiplicity = parse_data(name + '_' +
                                                                                                    spc + ext)
                self.cpu = sp_cpu(name)
            except ValueError:
                self.sp_energy = '!'
                pass
        elif spc == 'link':
            self.sp_energy, self.sp_program, self.sp_version_program, self.sp_solvation_model, self.sp_file, \
                self.sp_charge, self.sp_empirical_dispersion, self.sp_multiplicity = parse_data(file)
        # Count number of links
        for g_line in g_output:
            # Only read first link + freq not other link jobs
            if "Normal termination" in g_line:
                linkmax += 1
            else:
                frequency_wn = []
            if 'Frequencies --' in g_line:
                freq_loc = linkmax

        # Iterate over output
        molecular_mass = None  # make IDE happy
        if freq_loc == 0:
            freq_loc = len(g_output)
        for i, g_line in enumerate(g_output):
            # Link counter
            g_line = g_line.strip()
            if "Normal termination" in g_line:
                link += 1
                # Reset frequencies if in final freq link
                if link == freq_loc:
                    frequency_wn = []
                    im_frequency_wn = []
                    if mm_freq_scale_factor:
                        fract_model_sys = []
            # If spc specified will take last Energy from file, otherwise will break after freq calc
            if link > freq_loc:
                break
            # Iterate over output: look out for low frequencies
            if g_line.startswith('Frequencies -- '):
                new_line = None  # make IDE happy
                if mm_freq_scale_factor:
                    new_line = g_output[i + 3]
                for j in range(2, 5):
                    try:
                        x = float(g_line.split()[j])
                        # If given MM freq scale factor fill the fract_model_sys array:
                        if mm_freq_scale_factor:
                            y = float(new_line.strip().split()[j]) / 100.0
                            y = float('{:.6f}'.format(y))
                        else:
                            y = 1.0
                        # Only deal with real frequencies
                        if x > 0.00:
                            frequency_wn.append(x)
                            if mm_freq_scale_factor:
                                fract_model_sys.append(y)
                        # Check if we want to make any low lying imaginary frequencies positive
                        elif x < -1 * im_freq_cutoff:
                            if invert:
                                if x > float(invert):
                                    frequency_wn.append(x * -1.)
                                    inverted_freqs.append(x)
                                else:
                                    im_frequency_wn.append(x)
                            else:
                                im_frequency_wn.append(x)
                    except IndexError:
                        pass
            # For QM calculations look for SCF energies, last one will be the optimized energy
            elif g_line.startswith('SCF Done:'):
                self.scf_energy = float(g_line.split()[4])
            # For Counterpoise calculations the corrected energy value will be taken
            elif g_line.startswith('Counterpoise corrected energy'):
                self.scf_energy = float(g_line.split()[4])
            # For MP2 calculations replace with EUMP2
            elif 'EUMP2 =' in g_line:
                self.scf_energy = float((g_line.split()[5]).replace('D', 'E'))
            # For ONIOM calculations use the extrapolated value rather than SCF value
            elif "ONIOM: extrapolated energy" in g_line:
                self.scf_energy = float(g_line.split()[4])
            # For Semi-empirical or Molecular Mechanics calculations
            elif "Energy= " in g_line and "Predicted" not in g_line and "Thermal" not in g_line:
                self.scf_energy = float(g_line.split()[1])
            # Look for thermal corrections, paying attention to point group symmetry
            elif g_line.startswith('Zero-point correction='):
                self.zero_point_corr = float(g_line.split()[2])
            # Grab Multiplicity
            elif 'Multiplicity' in g_line:
                try:
                    self.mult = int(g_line.split('=')[-1].strip().split()[0])
                except IndexError:  # not really sure what error might happen, or why the below would fix it
                    self.mult = int(g_line.split()[-1])
            # Grab molecular mass
            elif g_line.startswith('Molecular mass:'):
                molecular_mass = float(g_line.split()[2])
            # Grab rational symmetry number
            elif g_line.startswith('Rotational symmetry number'):
                sym_no = int(g_line.split()[3].split(".")[0])
            # Grab point group
            elif g_line.startswith('Full point group'):
                if g_line.split()[3] == 'D*H' or g_line.split()[3] == 'C*V':
                    linear_mol = 1
            # Grab rotational constants
            elif g_line.startswith('Rotational constants (GHZ):'):
                try:
                    split_line = g_line.replace(':', ' ').split()
                    self.roconst = [float(split_line[3]), float(split_line[4]), float(split_line[5])]
                except ValueError:
                    if g_line.find('********'):
                        linear_warning = True
                        split_line = g_line.replace(':', ' ').split()
                        self.roconst = [float(split_line[4]), float(split_line[5])]
            # Grab rotational temperatures
            elif g_line.startswith('Rotational temperature '):
                rot_temp = [float(g_line.split()[3])]
            elif g_line.startswith('Rotational temperatures'):
                try:
                    split_line = g_line.split()
                    rot_temp = [float(split_line[3]), float(split_line[4]), float(split_line[5])]
                except ValueError:
                    rot_temp = None
                    if g_line.find('********'):
                        linear_warning = True
                        rot_temp = [float(g_line.split()[4]), float(g_line.split()[5])]
            if "Job cpu time" in g_line:
                split_line = g_line.split()
                # noinspection PyUnresolvedReferences
                days = int(split_line[3]) + self.cpu[0]
                hours = int(split_line[5]) + self.cpu[1]
                mins = int(split_line[7]) + self.cpu[2]
                secs = 0 + self.cpu[3]
                msecs = int(float(split_line[9]) * 1000.0) + self.cpu[4]
                self.cpu = [days, hours, mins, secs, msecs]
        self.inverted_freqs = inverted_freqs
        # Skip the calculation if unable to parse the frequencies or zpe from the output file
        if hasattr(self, "zero_point_corr") and rot_temp:
            cutoffs = [s_freq_cutoff] * len(frequency_wn)

            # Translational and electronic contributions to the energy and entropy do not depend on frequencies
            u_trans = calc_translational_energy(temperature)
            s_trans = calc_translational_entropy(molecular_mass, conc, temperature, solv)
            s_elec = calc_electronic_entropy(self.mult)

            # Rotational and Vibrational contributions to the energy entropy
            if len(frequency_wn) > 0:
                zpe = calc_zeropoint_energy(frequency_wn, zpe_scale_factor, fract_model_sys)
                u_rot = calc_rotational_energy(self.zero_point_corr, temperature, linear_mol)
                u_vib = calc_vibrational_energy(frequency_wn, temperature, freq_scale_factor, fract_model_sys)
                s_rot = calc_rotational_entropy(self.zero_point_corr, linear_mol, sym_no, rot_temp, temperature)

                # Calculate harmonic entropy, free-rotor entropy and damping function for each frequency
                s_vib_rrho = calc_rrho_entropy(frequency_wn, temperature, freq_scale_factor, fract_model_sys)

                if s_freq_cutoff > 0.0:
                    s_vib_rrqho = calc_rrho_entropy(cutoffs, temperature, freq_scale_factor, fract_model_sys)
                s_vib_free_rot = calc_free_rot_entropy(frequency_wn, temperature, freq_scale_factor, fract_model_sys)
                s_damp = calc_damp(frequency_wn, s_freq_cutoff)

                # check for qh
                if qh:
                    u_vib_qrrho = calc_q_rrho_energy(frequency_wn, temperature, freq_scale_factor)
                    h_damp = calc_damp(frequency_wn, h_freq_cutoff)

                # Compute entropy (cal/mol/K) using the two values and damping function
                vib_entropy = []
                vib_energy = []
                for j in range(0, len(frequency_wn)):
                    # Entropy correction
                    if qs == "grimme":
                        vib_entropy.append(s_vib_rrho[j] * s_damp[j] + (1 - s_damp[j]) * s_vib_free_rot[j])
                    elif qs == "truhlar":
                        if s_freq_cutoff > 0.0:
                            if frequency_wn[j] > s_freq_cutoff:
                                vib_entropy.append(s_vib_rrho[j])
                            else:
                                # yes, this is different from above, by one letter
                                vib_entropy.append(s_vib_rrqho[j])
                        else:
                            vib_entropy.append(s_vib_rrho[j])
                    # Enthalpy correction
                    if qh:
                        vib_energy.append(h_damp[j] * u_vib_qrrho[j] + (1 - h_damp[j]) * 0.5 *
                                          GAS_CONSTANT * temperature)

                qh_s_vib, h_s_vib = sum(vib_entropy), sum(s_vib_rrho)
                if qh:
                    qh_u_vib = sum(vib_energy)
            else:
                zpe, u_rot, u_vib, qh_u_vib, s_rot, h_s_vib, qh_s_vib = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0

            # The D3 term is added to the energy term here. If not requested then this term is zero
            # It is added to the SPC energy if defined (instead of the SCF energy)
            if spc:
                self.sp_energy += d3_energy
            else:
                self.scf_energy += d3_energy

            # Add terms (converted to au) to get Free energy - perform separately
            # for harmonic and quasi-harmonic values out of interest
            self.enthalpy = self.scf_energy + (u_trans + u_rot + u_vib + GAS_CONSTANT * temperature) / AU_TO_J
            self.qh_enthalpy = 0.0
            if qh:
                self.qh_enthalpy = self.scf_energy + (u_trans + u_rot + qh_u_vib + GAS_CONSTANT * temperature) / AU_TO_J
            # Single point correction replaces energy from optimization with single point value
            if spc:
                try:
                    self.enthalpy = self.enthalpy - self.scf_energy + self.sp_energy
                except TypeError:
                    pass
                if qh:
                    try:
                        self.qh_enthalpy = self.qh_enthalpy - self.scf_energy + self.sp_energy
                    except TypeError:
                        pass

            self.zpe = zpe / AU_TO_J
            self.entropy = (s_trans + s_rot + h_s_vib + s_elec) / AU_TO_J
            self.qh_entropy = (s_trans + s_rot + qh_s_vib + s_elec) / AU_TO_J

            # Symmetry - entropy correction for molecular symmetry
            if ssymm:
                sym_entropy_correction, p_group = self.sym_correction(file.split('.')[0].replace('/', '_'))
                self.point_group = p_group
                self.entropy += sym_entropy_correction
                self.qh_entropy += sym_entropy_correction

            # Calculate Free Energy
            if qh:
                self.gibbs_free_energy = self.enthalpy - temperature * self.entropy
                self.qh_gibbs_free_energy = self.qh_enthalpy - temperature * self.qh_entropy
            else:
                self.gibbs_free_energy = self.enthalpy - temperature * self.entropy
                self.qh_gibbs_free_energy = self.enthalpy - temperature * self.qh_entropy

            if cosmo:
                self.cosmo_qhg = self.qh_gibbs_free_energy + cosmo
            self.im_freq = []
            for freq in im_frequency_wn:
                if freq < -1 * im_freq_cutoff:
                    self.im_freq.append(freq)
        self.frequency_wn = frequency_wn
        self.im_frequency_wn = im_frequency_wn
        self.linear_warning = linear_warning

    def int_sym(self):
        self.xyz.get_connectivity()
        cap = [1, 9, 17]
        neighbor = [5, 6, 7, 8, 14, 15, 16]
        int_sym = 1

        for i, xyz_row in enumerate(self.xyz.connectivity):
            if self.xyz.atom_nums[i] != 6:
                continue
            a_array = np.array(self.xyz.atom_nums)[xyz_row]
            if len(a_array == 4):
                neighbors = [x for x in a_array if x in neighbor]
                caps = [x for x in a_array if x in cap]
                if (len(neighbors) == 1) and (len(set(caps)) == 1):
                    int_sym *= 3
        return int_sym

    def sym_correction(self, file):
        ex_sym, p_group = self.ex_sym(file)
        int_sym = self.int_sym()
        sym_num = ex_sym * int_sym
        sym_correction = (-GAS_CONSTANT * np.log(sym_num)) / AU_TO_J
        return sym_correction, p_group


# Read molecule data from a compchem output file
# Currently supports Gaussian and ORCA output types
#
class GetOutData:
    def __init__(self, file):
        self.atom_types = None
        with open(file) as f:
            out_data = f.readlines()
        program = 'none'

        for o_line in out_data:
            if "Gaussian" in o_line:
                program = "Gaussian"
                break
            if "* O   R   C   A *" in o_line:
                program = "Orca"
                break

        # noinspection PyShadowingNames
        def get_freqs(self, outlines, n_atoms, f_format):
            self.FREQS = []
            self.REDUCED_MASS = []
            self.FORCE_CONST = []
            self.NORMAL_MODE = []
            freqs_so_far = 0
            if f_format == "Gaussian":
                for i in range(0, len(outlines)):
                    if outlines[i].find(" Frequencies -- ") > -1:
                        n_freqs = len(outlines[i].split())
                        for j in range(2, n_freqs):
                            self.FREQS.append(float(outlines[i].split()[j]))
                            self.NORMAL_MODE.append([])
                        for j in range(3, n_freqs + 1):
                            self.REDUCED_MASS.append(float(outlines[i + 1].split()[j]))
                        for j in range(3, n_freqs + 1):
                            self.FORCE_CONST.append(float(outlines[i + 2].split()[j]))

                        for j in range(0, n_atoms):
                            for k in range(0, n_freqs - 2):
                                self.NORMAL_MODE[(freqs_so_far + k)].append(
                                    [float(outlines[i + 5 + j].split()[3 * k + 2]),
                                     float(outlines[i + 5 + j].split()[3 * k + 3]),
                                     float(outlines[i + 5 + j].split()[3 * k + 4])])
                        freqs_so_far = freqs_so_far + n_freqs - 2

        # noinspection PyShadowingNames
        def getatom_types(self, outlines, program_name):
            if program_name == "Gaussian":
                for i, p_line in enumerate(outlines):
                    if "Input orientation" in p_line or "Standard orientation" in p_line:
                        self.atom_nums, self.atom_types, self.cartesians, self.atomic_types, carts = [], [], [], [], \
                                                                                                    outlines[i + 5:]
                        for j, c_line in enumerate(carts):
                            if "-------" in c_line:
                                break
                            split_line = c_line.split()
                            self.atom_nums.append(int(split_line[1]))
                            self.atom_types.append(element_id(int(split_line[1])))
                            self.atomic_types.append(int(split_line[2]))
                            if len(split_line) > 5:
                                self.cartesians.append([float(split_line[3]), float(split_line[4]),
                                                        float(split_line[5])])
                            else:
                                self.cartesians.append([float(split_line[2]), float(split_line[3]),
                                                        float(split_line[4])])
            if program_name == "Orca":
                for i, r_line in enumerate(outlines):
                    if "*" in r_line and ">" in r_line and "xyz" in r_line:
                        self.atom_nums, self.atom_types, self.cartesians, carts = [], [], [], outlines[i + 1:]
                        for j, c_line in enumerate(carts):
                            if ">" in c_line and "*" in c_line:
                                break
                            split_line = r_line.split()
                            if len(split_line) > 5:
                                self.cartesians.append([float(split_line[3]), float(split_line[4]),
                                                        float(split_line[5])])
                                self.atom_types.append(split_line[2])
                                self.atom_nums.append(element_id(split_line[2], num=True))
                            else:
                                self.cartesians.append([float(split_line[2]), float(split_line[3]),
                                                        float(split_line[4])])
                                self.atom_types.append(split_line[1])
                                self.atom_nums.append(element_id(split_line[1], num=True))

        getatom_types(self, out_data, program)
        # noinspection PyTypeChecker
        n_atoms = len(self.atom_types)

        # noinspection PyTypeChecker
        get_freqs(self, out_data, n_atoms, program)

    # Obtain molecule connectivity to be used for internal symmetry determination
    def get_connectivity(self):
        connectivity = []
        tolerance = 0.2

        for i, ai in enumerate(self.atom_types):
            atom_row = []
            for j, aj in enumerate(self.atom_types):
                if i == j:
                    continue
                cutoff = RADII[ai] + RADII[aj] + tolerance
                distance = np.linalg.norm(np.array(self.cartesians[i]) - np.array(self.cartesians[j]))
                if distance < cutoff:
                    atom_row.append(j)
            connectivity.append(atom_row)
            # noinspection PyAttributeOutsideInit
            self.connectivity = connectivity


def calc_rrho_entropy(frequency_wn, temperature, freq_scale_factor, fract_model_sys):
    """
    Rigid rotor harmonic oscillator (RRHO) entropy evaluation - this is the default treatment
    Entropic contributions (J/(mol*K)) according to a rigid-rotor
    harmonic-oscillator description for a list of vibrational modes
    Sv = RSum(hv/(kT(e^(hv/kT)-1) - ln(1-e^(-hv/kT)))
    """
    factor = get_factors(fract_model_sys, freq_scale_factor, frequency_wn, temperature)
    entropy = [entry * GAS_CONSTANT / (np.exp(entry) - 1) - GAS_CONSTANT * np.log(1 - np.exp(-entry))
               for entry in factor]
    return entropy


def get_factors(fract_model_sys, freq_scale_factor, frequency_wn, temperature=1.0):
    if fract_model_sys:
        freq_scale_factor = [freq_scale_factor[0] * fract_model_sys[i] + freq_scale_factor[1] *
                             (1.0 - fract_model_sys[i]) for i in range(len(fract_model_sys))]
        factor = [(PLANCK_CONST_JS * frequency_wn[i] * SPEED_OF_LIGHT * freq_scale_factor[i]) /
                  (KB * temperature) for i in range(len(frequency_wn))]
    else:
        factor = [(PLANCK_CONST_JS * freq * SPEED_OF_LIGHT * freq_scale_factor) / (KB * temperature)
                  for freq in frequency_wn]
    return factor


def job_type(file):
    # Read output for the level of theory and basis set used
    job = ''
    with open(file) as f:
        job_data = f.readlines()
    for line in job_data:
        if line.strip().find('\\SP\\') > -1:
            job += 'SP'
        if line.strip().find('\\FOpt\\') > -1:
            job += 'GS'
        if line.strip().find('\\FTS\\') > -1:
            job += 'TS'
        if line.strip().find('\\Freq\\') > -1:
            job += 'Freq'
    return job


def read_file_contents(f_name):
    """
    Checks that the file is either '.out' or '.log', and reads the contents if so
    :param f_name: str, the file name
    :return: f_contents, str, contents of the file
    """
    if os.path.exists(os.path.splitext(f_name)[0] + '.log') or os.path.exists(os.path.splitext(f_name)[0] + '.out'):
        with open(f_name) as f:
            f_contents = f.readlines()
    elif os.path.exists(f_name):
        raise ValueError(f"Expected file name to end in '.out' or '.out' for file: {f_name}")
    else:
        raise ValueError("File {} does not exist".format(f_name))
    return f_contents


# noinspection DuplicatedCode
def parse_data(file):
    # Read Gaussian output and obtain single point energy, program type,
    # program version, solvation_model, charge, empirical_dispersion, multiplicity
    spe, program = 'none', 'none'
    version_program, solvation_model, keyword_line = '', '', ''
    charge, multiplicity = None, None
    file_contents = read_file_contents(file)

    for line in file_contents:
        if "Gaussian" in line:
            program = "Gaussian"
            break
        if "* O   R   C   A *" in line:
            program = "Orca"
            break
    repeated_link1 = 0
    for line in file_contents:
        if program == "Gaussian":
            if line.strip().startswith('SCF Done:'):
                spe = float(line.strip().split()[4])
            if line.strip().startswith('Counterpoise corrected energy'):
                spe = float(line.strip().split()[4])
            # For MP2 calculations replace with EUMP2
            if 'EUMP2 =' in line.strip():
                spe = float((line.strip().split()[5]).replace('D', 'E'))
            # For ONIOM calculations use the extrapolated value rather than SCF value
            if "ONIOM: extrapolated energy" in line.strip():
                spe = (float(line.strip().split()[4]))
            # For Semi-empirical or Molecular Mechanics calculations
            if "Energy= " in line.strip() and "Predicted" not in line.strip() and "Thermal" not in line.strip():
                spe = (float(line.strip().split()[1]))
            if "Gaussian" in line and "Revision" in line and repeated_link1 == 0:
                for i in range(len(line.strip(",").split(",")) - 1):
                    version_program += line.strip(",").split(",")[i]
                    repeated_link1 = 1
                version_program = version_program[1:]
            if "Charge" in line.strip() and "Multiplicity" in line.strip():
                charge = line.strip("=").split()[2]
                multiplicity = line.strip('=').split()[5]
        if program == "Orca":
            if line.strip().startswith('FINAL SINGLE POINT ENERGY'):
                spe = float(line.strip().split()[4])
            if 'Program Version' in line.strip():
                version_program = "ORCA version " + line.split()[2]
            if "Total Charge" in line.strip() and "...." in line.strip():
                charge = int(line.strip("=").split()[-1])
            if "Multiplicity" in line.strip() and "...." in line.strip():
                multiplicity = int(line.strip("=").split()[-1])

    # Solvation model and empirical dispersion detection
    empirical_dispersion = None
    sorted_solvation_model = None
    display_solvation_model = None
    if 'Gaussian' in version_program.strip():
        for i, line in enumerate(file_contents):
            if '#' in line.strip():
                for j, d_line in enumerate(file_contents[i:i + 10]):
                    if '--' in line.strip():
                        break
                    else:
                        for k in range(len(d_line.strip().split("\n"))):
                            keyword_line += d_line.strip().split("\n")[k]
        keyword_line = keyword_line.lower()
        if 'scrf' not in keyword_line.strip():
            solvation_model = "gas phase"
        else:
            start_scrf = keyword_line.strip().find('scrf') + 4
            if '(' in keyword_line[start_scrf:start_scrf + 4]:
                start_scrf += keyword_line[start_scrf:start_scrf + 4].find('(') + 1
                end_scrf = keyword_line.find(")", start_scrf)
                display_solvation_model = "scrf=(" + ','.join(
                    keyword_line[start_scrf:end_scrf].lower().split(',')) + ')'
                sorted_solvation_model = "scrf=(" + ','.join(
                    sorted(keyword_line[start_scrf:end_scrf].lower().split(','))) + ')'
            else:
                if ' = ' in keyword_line[start_scrf:start_scrf + 4]:
                    start_scrf += keyword_line[start_scrf:start_scrf + 4].find(' = ') + 3
                elif ' =' in keyword_line[start_scrf:start_scrf + 4]:
                    start_scrf += keyword_line[start_scrf:start_scrf + 4].find(' =') + 2
                elif '=' in keyword_line[start_scrf:start_scrf + 4]:
                    start_scrf += keyword_line[start_scrf:start_scrf + 4].find('=') + 1
                end_scrf = keyword_line.find(" ", start_scrf)
                if end_scrf == -1:
                    display_solvation_model = "scrf=(" + ','.join(keyword_line[start_scrf:].lower().split(',')) + ')'
                    sorted_solvation_model = "scrf=(" + ','.join(
                        sorted(keyword_line[start_scrf:].lower().split(','))) + ')'
                else:
                    display_solvation_model = "scrf=(" + ','.join(
                        keyword_line[start_scrf:end_scrf].lower().split(',')) + ')'
                    sorted_solvation_model = "scrf=(" + ','.join(
                        sorted(keyword_line[start_scrf:end_scrf].lower().split(','))) + ')'
        if solvation_model != "gas phase":
            solvation_model = [sorted_solvation_model, display_solvation_model]
        empirical_dispersion = ''
        if keyword_line.strip().find('empiricaldispersion') == -1 and keyword_line.strip().find(
                'emp=') == -1 and keyword_line.strip().find('emp =') == -1 and keyword_line.strip().find('emp(') == -1:
            empirical_dispersion = "No empirical dispersion detected"
        elif keyword_line.strip().find('empiricaldispersion') > -1:
            start_emp_disp = keyword_line.strip().find('empiricaldispersion') + 19
            empirical_dispersion = get_emp_dispersion(keyword_line, start_emp_disp)
        elif keyword_line.strip().find('emp=') > -1 or keyword_line.strip().find(
                'emp =') > -1 or keyword_line.strip().find('emp(') > -1:
            # Check for temp keyword
            temp, emp_e, emp_p = False, False, False
            check_temp = keyword_line.strip().find('emp=')
            start_emp_disp = keyword_line.strip().find('emp=')
            if check_temp == -1:
                check_temp = keyword_line.strip().find('emp =')
                start_emp_disp = keyword_line.strip().find('emp =')
            if check_temp == -1:
                check_temp = keyword_line.strip().find('emp=(')
                start_emp_disp = keyword_line.strip().find('emp(')
            check_temp += -1
            if keyword_line[check_temp].lower() == 't':
                temp = True  # Look for a new one
                if keyword_line.strip().find('emp=', check_temp + 5) > -1:
                    emp_e = True
                    start_emp_disp = keyword_line.strip().find('emp=', check_temp + 5) + 3
                elif keyword_line.strip().find('emp =', check_temp + 5) > -1:
                    emp_e = True
                    start_emp_disp = keyword_line.strip().find('emp =', check_temp + 5) + 3
                elif keyword_line.strip().find('emp(', check_temp + 5) > -1:
                    emp_p = True
                    start_emp_disp = keyword_line.strip().find('emp(', check_temp + 5) + 3
                else:
                    empirical_dispersion = "No empirical dispersion detected"
            else:
                start_emp_disp += 3
            if (temp and emp_e) or (not temp and keyword_line.strip().find('emp=') > -1) or (
                    not temp and keyword_line.strip().find('emp =')):
                empirical_dispersion = get_emp_dispersion(keyword_line, start_emp_disp)
            elif (temp and emp_p) or (not temp and keyword_line.strip().find('emp(') > -1):
                start_emp_disp += keyword_line[start_emp_disp:start_emp_disp + 4].find('(') + 1
                end_emp_disp = keyword_line.find(")", start_emp_disp)
                empirical_dispersion = 'empiricaldispersion=(' + ','.join(
                    sorted(keyword_line[start_emp_disp:end_emp_disp].lower().split(','))) + ')'
    if 'ORCA' in version_program.strip():
        keyword_line_1 = "gas phase"
        keyword_line_2 = ''
        keyword_line_3 = ''
        for i, line in enumerate(file_contents):
            if 'CPCM SOLVATION MODEL' in line.strip():
                keyword_line_1 = "CPCM,"
            if 'SMD CDS free energy correction energy' in line.strip():
                keyword_line_2 = "SMD,"
            if "Solvent:              " in line.strip():
                keyword_line_3 = line.strip().split()[-1]
        solvation_model = keyword_line_1 + keyword_line_2 + keyword_line_3
        empirical_dispersion1 = 'No empirical dispersion detected'
        empirical_dispersion2 = ''
        empirical_dispersion3 = ''
        if keyword_line.strip().find('DFT DISPERSION CORRECTION') > -1:
            empirical_dispersion1 = ''
        if keyword_line.strip().find('DFTD3') > -1:
            empirical_dispersion2 = "D3"
        if keyword_line.strip().find('USING zero damping') > -1:
            empirical_dispersion3 = ' with zero damping'
        empirical_dispersion = empirical_dispersion1 + empirical_dispersion2 + empirical_dispersion3
    return spe, program, version_program, solvation_model, file, charge, empirical_dispersion, multiplicity


# noinspection DuplicatedCode
def get_emp_dispersion(keyword_line, start_emp_disp):
    if '(' in keyword_line[start_emp_disp:start_emp_disp + 4]:
        start_emp_disp += keyword_line[start_emp_disp:start_emp_disp + 4].find('(') + 1
        end_emp_disp = keyword_line.find(")", start_emp_disp)
        empirical_dispersion = 'empiricaldispersion=(' + ','.join(
            sorted(keyword_line[start_emp_disp:end_emp_disp].lower().split(','))) + ')'
    else:
        if ' = ' in keyword_line[start_emp_disp:start_emp_disp + 4]:
            start_emp_disp += keyword_line[start_emp_disp:start_emp_disp + 4].find(' = ') + 3
        elif ' =' in keyword_line[start_emp_disp:start_emp_disp + 4]:
            start_emp_disp += keyword_line[start_emp_disp:start_emp_disp + 4].find(' =') + 2
        elif '=' in keyword_line[start_emp_disp:start_emp_disp + 4]:
            start_emp_disp += keyword_line[start_emp_disp:start_emp_disp + 4].find('=') + 1
        end_emp_disp = keyword_line.find(" ", start_emp_disp)
        if end_emp_disp == -1:
            empirical_dispersion = "empiricaldispersion=(" + ','.join(
                sorted(keyword_line[start_emp_disp:].lower().split(','))) + ')'
        else:
            empirical_dispersion = "empiricaldispersion=(" + ','.join(
                sorted(keyword_line[start_emp_disp:end_emp_disp].lower().split(','))) + ')'
    return empirical_dispersion


def sp_cpu(f_name):
    # Read single-point output for cpu time
    spe, program, cpu = None, None, None
    file_contents = read_file_contents(f_name)

    for line in file_contents:
        if line.find("Gaussian") > -1:
            program = "Gaussian"
            break
        if line.find("* O   R   C   A *") > -1:
            program = "Orca"
            break

    for line in file_contents:
        if program == "Gaussian":
            if line.strip().find("Job cpu time") > -1:
                days = int(line.split()[3])
                hours = int(line.split()[5])
                mins = int(line.split()[7])
                secs = 0
                msecs = int(float(line.split()[9]) * 1000.0)
                cpu = [days, hours, mins, secs, msecs]
        if program == "Orca":
            if line.strip().find("TOTAL RUN TIME") > -1:
                days = int(line.split()[3])
                hours = int(line.split()[5])
                mins = int(line.split()[7])
                secs = int(line.split()[9])
                msecs = float(line.split()[11])
                cpu = [days, hours, mins, secs, msecs]
    return cpu


def calc_translational_energy(temperature):
    """
    # Translational energy evaluation
    # Depends on temperature

    Calculates the translational energy (J/mol) of an ideal gas
    i.e. non-interacting molecules so molar energy = Na * atomic energy.
    This approximation applies to all energies and entropies computed within
    e_trans = 3/2 RT!
    """
    energy = 1.5 * GAS_CONSTANT * temperature
    return energy


def calc_rotational_energy(zpe, temperature, linear):
    """
    # Rotational energy evaluation
    # Depends on molecular shape and temperature
    Calculates the rotational energy (J/mol)
    E_trans = 0 (atomic) ; RT (linear); 3/2 RT (non-linear)
    """
    if zpe == 0.0:
        energy = 0.0
    elif linear == 1:
        energy = GAS_CONSTANT * temperature
    else:
        energy = 1.5 * GAS_CONSTANT * temperature
    return energy


def calc_zeropoint_energy(frequency_wn, scale_factor, fract_model_sys):
    """
    # Vibrational Zero point energy evaluation
    # Depends on frequencies and scaling factor: default = 1.0
    Calculates the vibrational ZPE (J/mol)
    E_ZPE = Sum(0.5 hv/k)
    """
    factor = get_factors(fract_model_sys, scale_factor, frequency_wn)
    energy = [0.5 * entry * GAS_CONSTANT for entry in factor]
    return sum(energy)


def calc_translational_entropy(molecular_mass, conc, temperature, solv):
    """
    # Translational entropy evaluation
    # Depends on mass, concentration, temperature, solvent free space: default = 1000.0
    Calculates the translational entropic contribution (J/(mol*K)) of an ideal gas.
    Needs the molecular mass. Convert mass in amu to kg; conc in mol/l to number per m^3
    s_trans = R(Ln(2 pi mkT/h^2)^3/2(1/C)) + 1 + 3/2)
    """
    e_lambda = ((2.0 * np.pi * molecular_mass * AMU_TO_KG * KB * temperature) ** 0.5) / PLANCK_CONST_JS
    freespace = get_free_space(solv)
    n_dens = conc * 1000 * AVOGADRO_CONST / (freespace / 1000.0)
    entropy = GAS_CONSTANT * (2.5 + np.log(e_lambda ** 3 / n_dens))
    return entropy


def calc_electronic_entropy(multiplicity):
    """
    # Electronic entropy evaluation
    # Depends on multiplicity
    Calculates the electronic entropic contribution (J/(mol*K)) of the molecule
    S_elec = R(Ln(multiplicity)
    """
    entropy = GAS_CONSTANT * (np.log(multiplicity))
    return entropy


def calc_vibrational_energy(frequency_wn, temperature, freq_scale_factor, fract_model_sys):
    """
    # Vibrational energy evaluation
    # Depends on frequencies, temperature and scaling factor: default = 1.0
    Calculates the vibrational energy contribution (J/mol).
    Includes ZPE (0K) and thermal contributions
    E_vib = R * Sum(0.5 hv/k + (hv/k)/(e^(hv/KT)-1))
    """
    factor = get_factors(fract_model_sys, freq_scale_factor, frequency_wn, temperature)
    # Error occurs if T is too low when performing np.exp
    for entry in factor:
        if entry > np.log(sys.float_info.max):
            raise InvalidDataError("Temperature may be too low to calculate vibrational energy. "
                                   "Please adjust using the `-t` option and try again.\n")

    energy = [entry * GAS_CONSTANT * temperature * (0.5 + (1.0 / (np.exp(entry) - 1.0))) for entry in factor]

    return sum(energy)


def calc_rotational_entropy(zpe, linear, sym_no, rot_temp, temperature):
    """
    # Rotational entropy evaluation
    # Depends on molecular shape and temp.
    Calculates the rotational entropy (J/(mol*K))
    S_trans = 0 (atomic) ; R(Ln(q)+1) (linear); R(Ln(q)+3/2) (non-linear)
    """
    q_rot = None  # make IDE happy
    if rot_temp == [0.0, 0.0, 0.0] or zpe == 0.0:  # Monatomic
        entropy = 0.0
    else:
        if len(rot_temp) == 1:  # Diatomic or linear molecules
            linear = 1
            q_rot = temperature / rot_temp[0]
        elif len(rot_temp) == 2:  # Possible gaussian problem with linear triatomic
            linear = 2
        else:
            q_rot = np.pi * temperature ** 3 / (rot_temp[0] * rot_temp[1] * rot_temp[2])
            q_rot = q_rot ** 0.5
        if linear == 1:
            entropy = GAS_CONSTANT * (np.log(q_rot / sym_no) + 1)
        elif linear == 2:
            entropy = 0.0
        else:
            entropy = GAS_CONSTANT * (np.log(q_rot / sym_no) + 1.5)
    return entropy


def calc_q_rrho_energy(frequency_wn, temperature, freq_scale_factor):
    """
    # Quasi-rigid rotor harmonic oscillator energy evaluation
    # used for calculating quasi-harmonic enthalpy
    Head-Gordon RRHO-vibrational energy contribution (J/mol*K) of
    vibrational modes described by a rigid-rotor harmonic approximation
    V_RRHO = 1/2(Nhv) + RT(hv/kT)e^(-hv/kT)/(1-e^(-hv/kT))
    """
    factor = [PLANCK_CONST_JS * freq * SPEED_OF_LIGHT * freq_scale_factor for freq in frequency_wn]
    energy = [0.5 * AVOGADRO_CONST * entry + GAS_CONSTANT * temperature * entry / KB
              / temperature * np.exp(-entry / KB / temperature) /
              (1 - np.exp(-entry / KB / temperature)) for entry in factor]
    return energy


def calc_free_rot_entropy(frequency_wn, temperature, freq_scale_factor, fract_model_sys):
    """
    # Free rotor entropy evaluation
    # used for low frequencies below the cut-off if qs=grimme is specified
    Entropic contributions (J/(mol*K)) according to a free-rotor
    description for a list of vibrational modes
    Sr = R(1/2 + 1/2ln((8pi^3u'kT/h^2))
    """
    # This is the average moment of inertia used by Grimme
    bav = 1.00e-44
    if fract_model_sys:
        freq_scale_factor = [freq_scale_factor[0] * fract_model_sys[i] + freq_scale_factor[1] *
                             (1.0 - fract_model_sys[i]) for i in range(len(fract_model_sys))]
        mu = [PLANCK_CONST_JS / (8 * np.pi ** 2 * frequency_wn[i] * SPEED_OF_LIGHT * freq_scale_factor[i]) for i in
              range(len(frequency_wn))]
    else:
        mu = [PLANCK_CONST_JS / (8 * np.pi ** 2 * freq * SPEED_OF_LIGHT * freq_scale_factor) for freq in frequency_wn]
    mu_primed = [entry * bav / (entry + bav) for entry in mu]
    factor = [8 * np.pi ** 3 * entry * KB * temperature / PLANCK_CONST_JS ** 2 for entry in mu_primed]
    entropy = [(0.5 + np.log(entry ** 0.5)) * GAS_CONSTANT for entry in factor]
    return entropy


# A damping function to interpolate between RRHO and free rotor vibrational entropy values
def calc_damp(frequency_wn, freq_cutoff):
    alpha = 4
    damp = [1 / (1 + (freq_cutoff / entry) ** alpha) for entry in frequency_wn]
    return damp


def get_free_space(solv):
    """
    # Computed the amount of accessible free space (ml per L) in solution
    # accessible to a solute immersed in bulk solvent, i.e. this is the volume
    # not occupied by solvent molecules, calculated using literature values for
    # molarity and B3LYP/6-31G* computed molecular volumes.
    Calculates the free space in a litre of bulk solvent, based on
    Shakhnovich and Whitesides (J. Org. Chem. 1998, 63, 3821-3830)
    """
    # todo: add more solvents???
    solvent_list = ["none", "H2O", "toluene", "DMF", "AcOH", "chloroform"]
    molarity = [1.0, 55.6, 9.4, 12.9, 17.4, 12.5]  # mol/l
    molecular_vol = [1.0, 27.944, 149.070, 77.442, 86.10, 97.0]  # Angstrom^3

    n_solv = 0
    for i in range(0, len(solvent_list)):
        if solv == solvent_list[i]:
            n_solv = i
    solv_molarity = molarity[n_solv]
    solv_volume = molecular_vol[n_solv]
    if n_solv > 0:
        v_free = 8 * ((1e27 / (solv_molarity * AVOGADRO_CONST)) ** 0.333333 - solv_volume ** 0.333333) ** 3
        freespace = v_free * solv_molarity * AVOGADRO_CONST * 1e-24
    else:
        freespace = 1000.0
    return freespace
