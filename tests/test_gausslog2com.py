import unittest
import os
from nrel_tools.gausslog2com import main
from nrel_tools.common import diff_lines, silent_remove, capture_stdout, capture_stderr
import logging

# logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger(__name__)
DISABLE_REMOVE = logger.isEnabledFor(logging.DEBUG)

__author__ = 'hmayes'

TEST_DIR = os.path.dirname(__file__)
MAIN_DIR = os.path.dirname(TEST_DIR)
DATA_DIR = os.path.join(os.path.dirname(__file__), 'test_data')
SUB_DATA_DIR = os.path.join(DATA_DIR, 'gausslog2com')

COM_TPL = os.path.join(SUB_DATA_DIR, 'cp.tpl')
LOG_LIST = os.path.join(SUB_DATA_DIR, 'log_list.txt')
COM1_OUT = os.path.join(SUB_DATA_DIR, 'pet_mono_843_tzvp_cp.com')
GOOD_COM1_OUT = os.path.join(SUB_DATA_DIR, 'pet_mono_843_tzvp_cp_good.com')
COM2_OUT = os.path.join(SUB_DATA_DIR, 'pet_mono_901_tzvp_cp.com')
GOOD_COM2_OUT = os.path.join(SUB_DATA_DIR, 'pet_mono_901_tzvp_cp_good.com')

LOG_FILE = os.path.join(SUB_DATA_DIR, 'pet_mono_1_tzvp.log')
COM3_OUT = os.path.join(SUB_DATA_DIR, 'pet_mono_1_tzvp_cp.com')
GOOD_COM3_OUT = os.path.join(SUB_DATA_DIR, 'pet_mono_1_tzvp_cp_good.com')


class Testgausslog2comNoOut(unittest.TestCase):
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


class Testgausslog2com(unittest.TestCase):
    # These test/demonstrate different options
    def testCPTpl(self):
        test_input = ["-t", COM_TPL, "-l", LOG_LIST]
        try:
            main(test_input)
            self.assertFalse(diff_lines(COM1_OUT, GOOD_COM1_OUT))
            self.assertFalse(diff_lines(COM2_OUT, GOOD_COM2_OUT))
        finally:
            silent_remove(COM1_OUT, disable=DISABLE_REMOVE)
            silent_remove(COM2_OUT, disable=DISABLE_REMOVE)
            pass

    def testFileCPTpl(self):
        test_input = ["-t", COM_TPL, "-f", LOG_FILE]
        try:
            main(test_input)
            self.assertFalse(diff_lines(COM3_OUT, GOOD_COM3_OUT))
        finally:
            silent_remove(COM3_OUT, disable=DISABLE_REMOVE)
            pass