import unittest
import os
from shutil import copyfile
from nrel_tools.check_gauss import main
from common_wrangler.common import capture_stdout, capture_stderr, diff_lines, silent_remove
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
SINGLE_FILE = os.path.join(SUB_DATA_DIR, 'me2propprpnt_7.log')
LIST_FILE = os.path.join(SUB_DATA_DIR, 'list.txt')

GOOD_OUT = "The following files completed normally:\n" \
           "    tests/test_data/check_gauss/pet_mono_637_tzvp.log\n" \
           "The following files may have failed:\n" \
           "    tests/test_data/check_gauss/me2propprpnt_7.log\n" \
           "    tests/test_data/check_gauss/pet_mono_674_tzvp.log\n" \
           "    tests/test_data/check_gauss/pet_mono_819_tzvp.log\n" \
           "    tests/test_data/check_gauss/pet_mono_872_tzvp.log\n" \
           "The following files may still be running:\n" \
           "    tests/test_data/check_gauss/pet_mono_671_tzvp.log\n"


class TestCheckGaussNoOut(unittest.TestCase):
    # These test failure cases that do not produce output files
    def testNoArgs(self):
        test_input = []
        # main(test_input)
        with capture_stderr(main, test_input) as output:
            self.assertTrue("Could not find files" in output)
        silent_remove(FOR_HARTREE_DIR, disable=DISABLE_REMOVE)

    def testHelp(self):
        test_input = ['-h']
        if logger.isEnabledFor(logging.DEBUG):
            main(test_input)
        with capture_stderr(main, test_input) as output:
            self.assertFalse(output)
        with capture_stdout(main, test_input) as output:
            self.assertTrue("optional arguments" in output)
        silent_remove(FOR_HARTREE_DIR, disable=DISABLE_REMOVE)

    def testWrongKey(self):
        test_input = ['-ghost']
        if logger.isEnabledFor(logging.DEBUG):
            main(test_input)
        with capture_stderr(main, test_input) as output:
            self.assertTrue("unrecognized arguments" in output)
        with capture_stdout(main, test_input) as output:
            self.assertTrue("optional arguments" in output)
        silent_remove(FOR_HARTREE_DIR, disable=DISABLE_REMOVE)

    def testWrongDir(self):
        test_input = ["-d", "ghost"]
        if logger.isEnabledFor(logging.DEBUG):
            main(test_input)
        with capture_stderr(main, test_input) as output:
            self.assertTrue("Could not find" in output)
        silent_remove(FOR_HARTREE_DIR, disable=DISABLE_REMOVE)

    def testConflictingOptions(self):
        test_input = ["-s", "-z"]
        if logger.isEnabledFor(logging.DEBUG):
            main(test_input)
        with capture_stderr(main, test_input) as output:
            self.assertTrue("Choose either" in output)
        silent_remove(FOR_HARTREE_DIR, disable=DISABLE_REMOVE)

    def testDirInsteadOfFile(self):
        test_input = ["-f", SUB_DATA_DIR, "-z"]
        try:
            with capture_stderr(main, test_input) as output:
                self.assertTrue("Is a directory" in output)
        finally:
            silent_remove(FOR_HARTREE_DIR, disable=DISABLE_REMOVE)
            pass


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
            with capture_stderr(main, test_input) as output:
                self.assertTrue("not read" in output)
        finally:
            silent_remove(MOVED_FILE, disable=DISABLE_REMOVE)
            silent_remove(FOR_HARTREE_DIR, disable=DISABLE_REMOVE)
            pass

    def testSingleFinalConvergence(self):
        test_input = ["-f", SINGLE_FILE, "-z"]
        good_out = 'File                                 Convergence Convergence_Error\n'\
                   'me2propprpnt_7.log                      111.4981 True\n'
        with capture_stdout(main, test_input) as output:
            self.assertTrue(output == good_out)

    def testListFinalConvergence(self):
        test_input = ["-l", LIST_FILE, "-z"]
        main(test_input)
        good_out = 'File                                 Convergence Convergence_Error\n' \
                   'hexyl_acrylate_419.log                    0.0706 False\n'\
                   'hexyl_acrylate_239.log                    1.1100 False\n'
        with capture_stdout(main, test_input) as output:
            self.assertTrue(output == good_out)
