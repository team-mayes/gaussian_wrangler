import unittest
import os
from shutil import copyfile
from nrel_tools.check_gauss import main
from nrel_tools.common import capture_stdout, capture_stderr, diff_lines, silent_remove
import logging

# logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger(__name__)
DISABLE_REMOVE = logger.isEnabledFor(logging.DEBUG)

__author__ = 'hmayes'

TEST_DIR = os.path.dirname(__file__)
MAIN_DIR = os.path.dirname(TEST_DIR)
DATA_DIR = os.path.join(TEST_DIR, 'test_data')
SUB_DATA_DIR = os.path.join(DATA_DIR, 'check_gauss')

NORM_TERM_LOG = os.path.join(SUB_DATA_DIR, 'pet_mono_637_tzvp.tpl')
TEMP_NORM_TERM_LOG = os.path.join(SUB_DATA_DIR, 'pet_mono_637_tzvp.log')
FOR_HARTREE_DIR = os.path.join(MAIN_DIR, 'for_hartree')
MOVED_FILE = os.path.join(FOR_HARTREE_DIR, 'pet_mono_637_tzvp.log')

GOOD_OUT = "The following files may have failed:\n" \
           "    /Users/hmayes/bee/code/python/nrel_tools/tests/test_data/check_gauss/pet_mono_819_tzvp.log\n" \
           "    /Users/hmayes/bee/code/python/nrel_tools/tests/test_data/check_gauss/pet_mono_872_tzvp.log\n" \
           "The following files may still be running:\n" \
           "    /Users/hmayes/bee/code/python/nrel_tools/tests/test_data/check_gauss/pet_mono_674_tzvp.log\n" \
           "    /Users/hmayes/bee/code/python/nrel_tools/tests/test_data/check_gauss/pet_mono_671_tzvp.log\n"


class TestCheckGaussNoOut(unittest.TestCase):
    # These test failure cases that do not produce output files
    def testNoArgs(self):
        test_input = []
        with capture_stderr(main, test_input) as output:
            self.assertTrue("WARNING:  Problems reading data: Could not find files" in output)

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
        silent_remove(MOVED_FILE, disable=DISABLE_REMOVE)
        try:
            # main(test_input)
            with capture_stdout(main, test_input) as output:
                self.assertTrue(output == GOOD_OUT)
            self.assertFalse(diff_lines(MOVED_FILE, NORM_TERM_LOG))
        finally:
            silent_remove(MOVED_FILE, disable=DISABLE_REMOVE)
            pass
