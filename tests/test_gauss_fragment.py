import unittest
import os
from gaussian_wrangler.gauss_fragment import main
from common_wrangler.common import diff_lines, silent_remove, capture_stdout, capture_stderr
import logging

# logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger(__name__)
DISABLE_REMOVE = logger.isEnabledFor(logging.DEBUG)

__author__ = 'hmayes'

TEST_DIR = os.path.dirname(__file__)
MAIN_DIR = os.path.dirname(TEST_DIR)
DATA_DIR = os.path.join(os.path.dirname(__file__), 'test_data')
SUB_DATA_DIR = os.path.join(DATA_DIR, 'gauss_fragment')
TEMP_DIR = os.path.join(SUB_DATA_DIR, 'temp_dir')

MISSING_FILE_INI = os.path.join(SUB_DATA_DIR, 'ghost_frag.ini')

DEF_INI = os.path.join(SUB_DATA_DIR, 'gausscom_fragment.ini')
F1_15_14_OUT = os.path.join(SUB_DATA_DIR, 'pet_mono_1_tzvp_15_14_f1.com')
F2_15_14_OUT = os.path.join(SUB_DATA_DIR, 'pet_mono_1_tzvp_15_14_f2.com')
CP_15_14_OUT = os.path.join(SUB_DATA_DIR, 'pet_mono_1_tzvp_15_14_cp.com')
GOOD_CP_15_14_OUT = os.path.join(SUB_DATA_DIR, 'pet_mono_1_tzvp_cp_good.com')
CP_14_15_OUT = os.path.join(SUB_DATA_DIR, 'pet_mono_1_tzvp_14_15_cp.com')
GOOD_CP_14_15_OUT = os.path.join(SUB_DATA_DIR, 'pet_mono_1_tzvp_14_15_cp_good.com')
F1_14_15_OUT = os.path.join(SUB_DATA_DIR, 'pet_mono_1_tzvp_14_15_f1.com')
GOOD_F1_14_15_OUT = os.path.join(SUB_DATA_DIR, 'pet_mono_1_tzvp_14_15_f1_good.com')
F2_14_15_OUT = os.path.join(SUB_DATA_DIR, 'pet_mono_1_tzvp_14_15_f2.com')
GOOD_F2_14_15_OUT = os.path.join(SUB_DATA_DIR, 'pet_mono_1_tzvp_14_15_f2_good.com')

LONELY_INI = os.path.join(SUB_DATA_DIR, 'gausscom_lonely_fragments.ini')
CP_20_21_OUT = os.path.join(SUB_DATA_DIR, 'pet_mono_1_tzvp_20_21_cp.com')
GOOD_CP_20_21_OUT = os.path.join(SUB_DATA_DIR, 'pet_mono_1_tzvp_20_21_cp_good.com')
F2_20_21_OUT = os.path.join(SUB_DATA_DIR, 'pet_mono_1_tzvp_20_21_f2.com')
GOOD_F2_20_21_OUT = os.path.join(SUB_DATA_DIR, 'pet_mono_1_tzvp_20_21_f2_good.com')

CP_22_23_OUT = os.path.join(SUB_DATA_DIR, 'pet_mono_1_tzvp_22_23_cp.com')
GOOD_CP_22_23_OUT = os.path.join(SUB_DATA_DIR, 'pet_mono_1_tzvp_22_23_cp_good.com')
F2_22_23_OUT = os.path.join(SUB_DATA_DIR, 'pet_mono_1_tzvp_22_23_f2.com')
GOOD_F2_22_23_OUT = os.path.join(SUB_DATA_DIR, 'pet_mono_1_tzvp_22_23_f2_good.com')

CP_24_18_OUT = os.path.join(SUB_DATA_DIR, 'pet_mono_1_tzvp_24_18_cp.com')
GOOD_CP_24_18_OUT = os.path.join(SUB_DATA_DIR, 'pet_mono_1_tzvp_24_18_cp_good.com')
F2_24_18_OUT = os.path.join(SUB_DATA_DIR, 'pet_mono_1_tzvp_24_18_f2.com')
GOOD_F2_24_18_OUT = os.path.join(SUB_DATA_DIR, 'pet_mono_1_tzvp_24_18_f2_good.com')

DIMER_INI = os.path.join(SUB_DATA_DIR, 'gauss_frag_dimer.ini')
CP_DI_18_24_OUT = os.path.join(SUB_DATA_DIR, 'pet_dimer_tzvp_18_24_cp.com')
GOOD_CP_DI_18_24_OUT = os.path.join(SUB_DATA_DIR, 'pet_dimer_tzvp_18_24_cp_good.com')
F2_DI_18_24_OUT = os.path.join(SUB_DATA_DIR, 'pet_dimer_tzvp_18_24_f2.com')
GOOD_F2_DI_18_24_OUT = os.path.join(SUB_DATA_DIR, 'pet_dimer_tzvp_18_24_f2_good.com')
CP_DI_18_16_OUT = os.path.join(SUB_DATA_DIR, 'pet_dimer_tzvp_18_16_cp.com')
GOOD_CP_DI_18_16_OUT = os.path.join(SUB_DATA_DIR, 'pet_dimer_tzvp_18_16_cp_good.com')
F1_DI_18_16_OUT = os.path.join(SUB_DATA_DIR, 'pet_dimer_tzvp_18_16_f1.com')
GOOD_F1_DI_18_16_OUT = os.path.join(SUB_DATA_DIR, 'pet_dimer_tzvp_18_16_f1_good.com')
F2_DI_18_16_OUT = os.path.join(SUB_DATA_DIR, 'pet_dimer_tzvp_18_16_f2.com')
GOOD_F2_DI_18_16_OUT = os.path.join(SUB_DATA_DIR, 'pet_dimer_tzvp_18_16_f2_good.com')

CC_INI = os.path.join(SUB_DATA_DIR, 'tbut_frag.ini')
CP_CC_12_11_OUT = os.path.join(SUB_DATA_DIR, 'tbut_12_11_cp.com')
GOOD_CP_CC_12_11_OUT = os.path.join(SUB_DATA_DIR, 'tbut_12_11_cp_good.com')
F1_CC_12_11_OUT = os.path.join(SUB_DATA_DIR, 'tbut_12_11_f1.com')
GOOD_F1_CC_12_11_OUT = os.path.join(SUB_DATA_DIR, 'tbut_12_11_f1_good.com')
F2_CC_12_11_OUT = os.path.join(SUB_DATA_DIR, 'tbut_12_11_f2.com')
GOOD_F2_CC_12_11_OUT = os.path.join(SUB_DATA_DIR, 'tbut_12_11_f2_good.com')
CP_CC_12_38_OUT = os.path.join(SUB_DATA_DIR, 'tbut_12_38_cp.com')
GOOD_CP_CC_12_38_OUT = os.path.join(SUB_DATA_DIR, 'tbut_12_38_cp_good.com')
F2_CC_12_38_OUT = os.path.join(SUB_DATA_DIR, 'tbut_12_38_f2.com')
GOOD_F2_CC_12_38_OUT = os.path.join(SUB_DATA_DIR, 'tbut_12_38_f2_good.com')
CP_CC_13_12_OUT = os.path.join(SUB_DATA_DIR, 'tbut_13_12_cp.com')
GOOD_CP_CC_13_12_OUT = os.path.join(SUB_DATA_DIR, 'tbut_13_12_cp_good.com')
F1_CC_13_12_OUT = os.path.join(SUB_DATA_DIR, 'tbut_13_12_f1.com')
GOOD_F1_CC_13_12_OUT = os.path.join(SUB_DATA_DIR, 'tbut_13_12_f1_good.com')
F2_CC_13_12_OUT = os.path.join(SUB_DATA_DIR, 'tbut_13_12_f2.com')
GOOD_F2_CC_13_12_OUT = os.path.join(SUB_DATA_DIR, 'tbut_13_12_f2_good.com')

IGNORE_MAX_INI = os.path.join(SUB_DATA_DIR, 'gauss_frag_ignore_max_dist.ini')
CP_TPA_OUT = os.path.join(SUB_DATA_DIR, 'tpaegh1ats_ts_ircf_opt_1_2_cp.com')
GOOD_CP_TPA_OUT = os.path.join(SUB_DATA_DIR, 'tpaegh1ats_ts_ircf_opt_1_2_cp_good.com')
F2_TPA_OUT = os.path.join(SUB_DATA_DIR, 'tpaegh1ats_ts_ircf_opt_1_2_f2.com')
GOOD_F2_TPA_OUT = os.path.join(SUB_DATA_DIR, 'tpaegh1ats_ts_ircf_opt_1_2_f2_good.com')

N_TRIPLE_BOND_INI = os.path.join(SUB_DATA_DIR, 'iso_frag.ini')
N_TRIPLE_BOND_SUB_DIR_INI = os.path.join(SUB_DATA_DIR, 'iso_frag_sub_dir.ini')
CP_N_OUT = os.path.join(SUB_DATA_DIR, 'initrile_16_8_cp.com')
GOOD_CP_N_OUT = os.path.join(SUB_DATA_DIR, 'initrile_16_8_cp_good.com')
F2_N_OUT = os.path.join(SUB_DATA_DIR, 'initrile_16_8_f2.com')
GOOD_F2_N_OUT = os.path.join(SUB_DATA_DIR, 'initrile_16_8_f2_good.com')
CP_N_OUT_SUB_DIR = os.path.join(TEMP_DIR, 'initrile_16_8_cp.com')
F2_N_OUT_SUB_DIR = os.path.join(TEMP_DIR, 'initrile_16_8_f2.com')

METAL_BOND_INI = os.path.join(SUB_DATA_DIR, 'gauss_frag_w_metal.ini')
CP_METAL_OUT = os.path.join(SUB_DATA_DIR, 'tieg5ipatse_ts_ircr_optts_37_38_cp.com')
GOOD_CP_METAL_OUT = os.path.join(SUB_DATA_DIR, 'tieg5ipatse_ts_ircr_optts_37_38_cp_good.com')
F1_METAL_OUT = os.path.join(SUB_DATA_DIR, 'tieg5ipatse_ts_ircr_optts_37_38_f1.com')
GOOD_F1_METAL_OUT = os.path.join(SUB_DATA_DIR, 'tieg5ipatse_ts_ircr_optts_37_38_f1_good.com')
F2_METAL_OUT = os.path.join(SUB_DATA_DIR, 'tieg5ipatse_ts_ircr_optts_37_38_f2.com')
GOOD_F2_METAL_OUT = os.path.join(SUB_DATA_DIR, 'tieg5ipatse_ts_ircr_optts_37_38_f2_good.com')

ADD_CP_FOOTER_INI = os.path.join(SUB_DATA_DIR, 'gauss_frag_add_end_basis.ini')
CP_FOOTER_CP_OUT = os.path.join(MAIN_DIR, '2011shi_fig5cts_origts_ircf_opt_19_40_cp.com')
CP_FOOTER_F1_OUT = os.path.join(MAIN_DIR, '2011shi_fig5cts_origts_ircf_opt_19_40_f1.com')
CP_FOOTER_F2_OUT = os.path.join(MAIN_DIR, '2011shi_fig5cts_origts_ircf_opt_19_40_f2.com')
GOOD_CP_FOOTER_CP_OUT = os.path.join(SUB_DATA_DIR, '2011shi_fig5cts_origts_ircf_opt_19_40_cp_good.com')
GOOD_CP_FOOTER_F1_OUT = os.path.join(SUB_DATA_DIR, '2011shi_fig5cts_origts_ircf_opt_19_40_f1_good.com')
GOOD_CP_FOOTER_F2_OUT = os.path.join(SUB_DATA_DIR, '2011shi_fig5cts_origts_ircf_opt_19_40_f2_good.com')

LOG_AS_COM_INI = os.path.join(SUB_DATA_DIR, 'gauss_log_as_com.ini')
LOG_AND_COM_INI = os.path.join(SUB_DATA_DIR, 'gauss_log_and_com.ini')
NO_LOG_OR_COM_INI = os.path.join(SUB_DATA_DIR, 'gauss_no_log_or_com.ini')
MORE_THAN_TW0_ATOM_PAIR_INI = os.path.join(SUB_DATA_DIR, 'gauss_frag_too_many_atoms.ini')
INVALID_ATOM_ID_INI = os.path.join(SUB_DATA_DIR, 'gauss_frag_invalid_atom_id.ini')
FAIL_TO_IGNORE_MAX_INI = os.path.join(SUB_DATA_DIR, 'gauss_frag_no_ignore_max_dist.ini')


class TestGausscomFragNoOut(unittest.TestCase):
    # These all test failure cases
    def testNoArgs(self):
        with capture_stderr(main, []) as output:
            self.assertTrue("WARNING:  Problems reading file: Could not read file" in output)
        with capture_stdout(main, []) as output:
            self.assertTrue("optional arguments" in output)

    def testHelp(self):
        test_input = ['-h']
        if logger.isEnabledFor(logging.DEBUG):
            main(test_input)
        with capture_stderr(main, test_input) as output:
            self.assertFalse(output)
        with capture_stdout(main, test_input) as output:
            self.assertTrue("optional arguments" in output)

    def testMissingFile(self):
        test_input = ['-c', MISSING_FILE_INI]
        if logger.isEnabledFor(logging.DEBUG):
            main(test_input)
        with capture_stderr(main, test_input) as output:
            self.assertTrue("No such file or directory" in output)

    def testLogAsCom(self):
        test_input = ['-c', LOG_AS_COM_INI]
        # main(test_input)
        with capture_stderr(main, test_input) as output:
            self.assertTrue("where charge and multiplicity are expected" in output)

    def testLogAndCom(self):
        test_input = ['-c', LOG_AND_COM_INI]
        # main(test_input)
        with capture_stderr(main, test_input) as output:
            self.assertTrue("Both an" in output)

    def testNoLogOrCom(self):
        test_input = ['-c', NO_LOG_OR_COM_INI]
        # main(test_input)
        with capture_stderr(main, test_input) as output:
            self.assertTrue("requires either" in output)

    def testMoreThan2AtomsInPair(self):
        test_input = ['-c', MORE_THAN_TW0_ATOM_PAIR_INI]
        # main(test_input)
        with capture_stderr(main, test_input) as output:
            self.assertTrue("should be sets of two atoms" in output)

    def testInvalidAtomID(self):
        test_input = ['-c', INVALID_ATOM_ID_INI]
        # main(test_input)
        with capture_stderr(main, test_input) as output:
            self.assertTrue("in 'cut_atoms'" in output)

    def testTooDistant(self):
        test_input = ['-c', FAIL_TO_IGNORE_MAX_INI]
        # main(test_input)
        with capture_stderr(main, test_input) as output:
            self.assertTrue("Angstroms apart" in output)


class TestGausscomFrag(unittest.TestCase):
    # These test/demonstrate different options
    def testDefIni(self):
        test_input = ["-c", DEF_INI]
        try:
            main(test_input)
            self.assertFalse(diff_lines(CP_15_14_OUT, GOOD_CP_15_14_OUT))
            self.assertFalse(diff_lines(CP_14_15_OUT, GOOD_CP_14_15_OUT))
            self.assertFalse(diff_lines(F1_14_15_OUT, GOOD_F1_14_15_OUT))
            self.assertFalse(diff_lines(F2_14_15_OUT, GOOD_F2_14_15_OUT))
            self.assertEqual(len(diff_lines(F1_15_14_OUT, GOOD_F1_14_15_OUT)), 2)
            self.assertEqual(len(diff_lines(F2_15_14_OUT, GOOD_F2_14_15_OUT)), 2)
        finally:
            silent_remove(CP_15_14_OUT, disable=DISABLE_REMOVE)
            silent_remove(CP_14_15_OUT, disable=DISABLE_REMOVE)
            silent_remove(F1_14_15_OUT, disable=DISABLE_REMOVE)
            silent_remove(F2_14_15_OUT, disable=DISABLE_REMOVE)
            silent_remove(F1_15_14_OUT, disable=DISABLE_REMOVE)
            silent_remove(F2_15_14_OUT, disable=DISABLE_REMOVE)
            pass

    def testLonelyFragments(self):
        test_input = ["-c", LONELY_INI]
        try:
            main(test_input)
            self.assertFalse(diff_lines(CP_20_21_OUT, GOOD_CP_20_21_OUT))
            self.assertFalse(diff_lines(CP_22_23_OUT, GOOD_CP_22_23_OUT))
            self.assertFalse(diff_lines(CP_24_18_OUT, GOOD_CP_24_18_OUT))
            self.assertFalse(diff_lines(F2_20_21_OUT, GOOD_F2_20_21_OUT))
            self.assertFalse(diff_lines(F2_22_23_OUT, GOOD_F2_22_23_OUT))
            self.assertFalse(diff_lines(F2_24_18_OUT, GOOD_F2_24_18_OUT))
        finally:
            silent_remove(CP_20_21_OUT, disable=DISABLE_REMOVE)
            silent_remove(CP_22_23_OUT, disable=DISABLE_REMOVE)
            silent_remove(CP_24_18_OUT, disable=DISABLE_REMOVE)
            silent_remove(F2_20_21_OUT, disable=DISABLE_REMOVE)
            silent_remove(F2_22_23_OUT, disable=DISABLE_REMOVE)
            silent_remove(F2_24_18_OUT, disable=DISABLE_REMOVE)
            pass

    def testDubO(self):
        test_input = ["-c", DIMER_INI]
        try:
            main(test_input)
            self.assertFalse(diff_lines(CP_DI_18_24_OUT, GOOD_CP_DI_18_24_OUT))
            self.assertFalse(diff_lines(F2_DI_18_24_OUT, GOOD_F2_DI_18_24_OUT))
            self.assertFalse(diff_lines(CP_DI_18_16_OUT, GOOD_CP_DI_18_16_OUT))
            self.assertFalse(diff_lines(F1_DI_18_16_OUT, GOOD_F1_DI_18_16_OUT))
            self.assertFalse(diff_lines(F2_DI_18_16_OUT, GOOD_F2_DI_18_16_OUT))
        finally:
            silent_remove(CP_DI_18_24_OUT, disable=DISABLE_REMOVE)
            silent_remove(F2_DI_18_24_OUT, disable=DISABLE_REMOVE)
            silent_remove(CP_DI_18_16_OUT, disable=DISABLE_REMOVE)
            silent_remove(F1_DI_18_16_OUT, disable=DISABLE_REMOVE)
            silent_remove(F2_DI_18_16_OUT, disable=DISABLE_REMOVE)
            pass

    def testDubCC(self):
        test_input = ["-c", CC_INI]
        try:
            main(test_input)
            self.assertFalse(diff_lines(CP_CC_12_11_OUT, GOOD_CP_CC_12_11_OUT))
            self.assertFalse(diff_lines(F1_CC_12_11_OUT, GOOD_F1_CC_12_11_OUT))
            self.assertFalse(diff_lines(F2_CC_12_11_OUT, GOOD_F2_CC_12_11_OUT))
            self.assertFalse(diff_lines(CP_CC_12_38_OUT, GOOD_CP_CC_12_38_OUT))
            self.assertFalse(diff_lines(F2_CC_12_38_OUT, GOOD_F2_CC_12_38_OUT))
            self.assertFalse(diff_lines(CP_CC_13_12_OUT, GOOD_CP_CC_13_12_OUT))
            self.assertFalse(diff_lines(F1_CC_13_12_OUT, GOOD_F1_CC_13_12_OUT))
            self.assertFalse(diff_lines(F2_CC_13_12_OUT, GOOD_F2_CC_13_12_OUT))
        finally:
            silent_remove(CP_CC_12_11_OUT, disable=DISABLE_REMOVE)
            silent_remove(F1_CC_12_11_OUT, disable=DISABLE_REMOVE)
            silent_remove(F2_CC_12_11_OUT, disable=DISABLE_REMOVE)
            silent_remove(CP_CC_12_38_OUT, disable=DISABLE_REMOVE)
            silent_remove(F2_CC_12_38_OUT, disable=DISABLE_REMOVE)
            silent_remove(CP_CC_13_12_OUT, disable=DISABLE_REMOVE)
            silent_remove(F1_CC_13_12_OUT, disable=DISABLE_REMOVE)
            silent_remove(F2_CC_13_12_OUT, disable=DISABLE_REMOVE)
            pass

    def testNTripleBond(self):
        test_input = ["-c", N_TRIPLE_BOND_INI]
        try:
            main(test_input)
            self.assertFalse(diff_lines(CP_N_OUT, GOOD_CP_N_OUT))
            self.assertFalse(diff_lines(F2_N_OUT, GOOD_F2_N_OUT))
        finally:
            silent_remove(CP_N_OUT, disable=DISABLE_REMOVE)
            silent_remove(F2_N_OUT, disable=DISABLE_REMOVE)
            pass

    def testNTripleBondSubDir(self):
        # Test that output files were created in subdir
        test_input = ["-c", N_TRIPLE_BOND_SUB_DIR_INI]
        try:
            silent_remove(TEMP_DIR, dir_with_files=True)
            main(test_input)
            self.assertFalse(diff_lines(CP_N_OUT_SUB_DIR, GOOD_CP_N_OUT))
            self.assertFalse(diff_lines(F2_N_OUT_SUB_DIR, GOOD_F2_N_OUT))
        finally:
            silent_remove(TEMP_DIR, disable=DISABLE_REMOVE, dir_with_files=True)
            pass

    def testSepMolecules(self):
        # to help when there is a reactant or product complex that is two different molecules; ignores the
        #   larger distance between the atoms belonging to two different molecules
        test_input = ["-c", IGNORE_MAX_INI]
        try:
            main(test_input)
            self.assertFalse(diff_lines(CP_TPA_OUT, GOOD_CP_TPA_OUT))
            self.assertFalse(diff_lines(F2_TPA_OUT, GOOD_F2_TPA_OUT))
        finally:
            silent_remove(CP_TPA_OUT, disable=DISABLE_REMOVE)
            silent_remove(F2_TPA_OUT, disable=DISABLE_REMOVE)
            pass

    def testMetalBond(self):
        # Uses a longer bond cut-off when one of the atoms is a metal
        test_input = ["-c", METAL_BOND_INI]
        try:
            main(test_input)
            self.assertFalse(diff_lines(CP_METAL_OUT, GOOD_CP_METAL_OUT))
            self.assertFalse(diff_lines(F1_METAL_OUT, GOOD_F1_METAL_OUT))
            self.assertFalse(diff_lines(F2_METAL_OUT, GOOD_F2_METAL_OUT))
        finally:
            silent_remove(CP_METAL_OUT, disable=DISABLE_REMOVE)
            silent_remove(F1_METAL_OUT, disable=DISABLE_REMOVE)
            silent_remove(F2_METAL_OUT, disable=DISABLE_REMOVE)
            pass

    def testCPFooter(self):
        # Checks adding a multiline footer
        # Also checks making the output directory the current directory using 'output_directory =' in the ini
        test_input = ["-c", ADD_CP_FOOTER_INI]
        try:
            main(test_input)
            self.assertFalse(diff_lines(CP_FOOTER_CP_OUT, GOOD_CP_FOOTER_CP_OUT))
            self.assertFalse(diff_lines(CP_FOOTER_F1_OUT, GOOD_CP_FOOTER_F1_OUT))
            self.assertFalse(diff_lines(CP_FOOTER_F2_OUT, GOOD_CP_FOOTER_F2_OUT))
        finally:
            silent_remove(CP_FOOTER_CP_OUT, disable=DISABLE_REMOVE)
            silent_remove(CP_FOOTER_F1_OUT, disable=DISABLE_REMOVE)
            silent_remove(CP_FOOTER_F2_OUT, disable=DISABLE_REMOVE)
            pass
