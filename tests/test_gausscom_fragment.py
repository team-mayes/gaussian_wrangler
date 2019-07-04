import unittest
import os
from nrel_tools.gausscom_fragment import main
from nrel_tools.common import diff_lines, silent_remove, capture_stdout, capture_stderr
import logging

# logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger(__name__)
DISABLE_REMOVE = logger.isEnabledFor(logging.DEBUG)

__author__ = 'hmayes'

TEST_DIR = os.path.dirname(__file__)
MAIN_DIR = os.path.dirname(TEST_DIR)
DATA_DIR = os.path.join(os.path.dirname(__file__), 'test_data')
SUB_DATA_DIR = os.path.join(DATA_DIR, 'gausscom_fragment')

DEF_INI = os.path.join(SUB_DATA_DIR, 'gausscom_fragment.ini')

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
    def testCPTpl(self):
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
