import unittest
import os
from shutil import copyfile
from nrel_tools.check_gauss import main
from nrel_tools.common import diff_lines, silent_remove, capture_stdout, capture_stderr
import logging

# logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger(__name__)
DISABLE_REMOVE = logger.isEnabledFor(logging.DEBUG)

__author__ = 'hmayes'

TEST_DIR = os.path.dirname(__file__)
MAIN_DIR = os.path.dirname(TEST_DIR)
DATA_DIR = os.path.join(os.path.dirname(__file__), 'test_data')
SUB_DATA_DIR = os.path.join(DATA_DIR, 'check_gauss')

NORM_TERM_LOG = os.path.join(SUB_DATA_DIR,'pet_mono_637_tzvp.tpl')
TEMP_NORM_TERM_LOG = os.path.join(SUB_DATA_DIR,'pet_mono_637_tzvp.log')


class TestCheckGaussNoOut(unittest.TestCase):
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


class TestCheckGauss(unittest.TestCase):
    # These test/demonstrate different options
    def testBasicUse(self):
        test_input = ["-d", SUB_DATA_DIR]
        copyfile(NORM_TERM_LOG, TEMP_NORM_TERM_LOG)
        try:
            main(test_input)
        finally:
            pass

