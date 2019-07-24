import unittest
import os
from nrel_tools.gauss_fragment import main
from nrel_tools.common import diff_lines, silent_remove, capture_stdout, capture_stderr
import logging

# logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger(__name__)
DISABLE_REMOVE = logger.isEnabledFor(logging.DEBUG)

__author__ = 'hmayes'

TEST_DIR = os.path.dirname(__file__)
MAIN_DIR = os.path.dirname(TEST_DIR)
DATA_DIR = os.path.join(os.path.dirname(__file__), 'test_data')
SUB_DATA_DIR = os.path.join(DATA_DIR, 'gauss_fragment')

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
