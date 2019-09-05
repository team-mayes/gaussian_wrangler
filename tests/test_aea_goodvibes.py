import unittest
import os
from nrel_tools.aea_goodvibes import main
from nrel_tools.common import diff_lines, silent_remove, capture_stdout, capture_stderr
import logging

# logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger(__name__)
DISABLE_REMOVE = logger.isEnabledFor(logging.DEBUG)

__author__ = 'hmayes'

TEST_DIR = os.path.dirname(__file__)
MAIN_DIR = os.path.dirname(TEST_DIR)
DATA_DIR = os.path.join(os.path.dirname(__file__), 'test_data')
SUB_DATA_DIR = os.path.join(DATA_DIR, 'aea_goodvibes')

GOODVIBES_DAT = os.path.abspath(os.path.join(TEST_DIR, '..', 'Goodvibes_output.dat'))
GOODVIBES_CSV = os.path.abspath(os.path.join(TEST_DIR, '..', 'Goodvibes_output.csv'))

UNI_REACT = os.path.join(SUB_DATA_DIR, 'ipaegh1dts_t_ircr_opt.log')
UNI_TS = os.path.join(SUB_DATA_DIR, 'ipaegh1dts.log')
FILE_LIST = os.path.join(SUB_DATA_DIR, 'list.txt')
FILE_LIST_MISMATCH_SOLV = os.path.join(SUB_DATA_DIR, 'list_mismatch_solv.txt')

AE_OUT = os.path.join(SUB_DATA_DIR, 'aea_out.csv')
GOOD_AE_OUT = os.path.join(SUB_DATA_DIR, 'aea_out_good.csv')

BI_LIST = os.path.join(SUB_DATA_DIR, 'list_bimolec.txt')
GOOD_AE_BI_OUT = os.path.join(SUB_DATA_DIR, 'aea_out_bi_good.csv')


class TestAEaGoodVibesNoOut(unittest.TestCase):
    # These all test failure cases
    def testNoArgs(self):
        test_input = []
        main(test_input)
        # with capture_stderr(main, test_input) as output:
        #     self.assertTrue("WARNING:  Problems reading file: Could not read file" in output)
        # with capture_stdout(main, test_input) as output:
        #     self.assertTrue("optional arguments" in output)

    def testHelp(self):
        test_input = ['-h']
        if logger.isEnabledFor(logging.DEBUG):
            main(test_input)
        with capture_stderr(main, test_input) as output:
            self.assertFalse(output)
        with capture_stdout(main, test_input) as output:
            self.assertTrue("optional arguments" in output)


class TestAEaGoodVibesInputError(unittest.TestCase):
    def testMixMatchSolvent(self):
        test_input = ["-l", FILE_LIST_MISMATCH_SOLV]
        try:
            # main(test_input)
            with capture_stderr(main, test_input) as output:
                self.assertTrue("GoodVibes error" in output)
        finally:
            silent_remove(GOODVIBES_DAT, disable=DISABLE_REMOVE)
            silent_remove(GOODVIBES_CSV, disable=DISABLE_REMOVE)
            pass

    def testMixMatchTheory(self):
        test_input = [UNI_REACT, UNI_TS]
        try:
            # main(test_input)
            with capture_stderr(main, test_input) as output:
                self.assertTrue("GoodVibes error" in output)
        finally:
            silent_remove(GOODVIBES_DAT, disable=DISABLE_REMOVE)
            silent_remove(GOODVIBES_CSV, disable=DISABLE_REMOVE)
            pass


class TestAEaGoodVibes(unittest.TestCase):
    # These test/demonstrate different options
    def testTwoUni(self):
        test_input = ["-l", FILE_LIST, "-d", SUB_DATA_DIR]
        try:
            main(test_input)
            self.assertFalse(diff_lines(AE_OUT, GOOD_AE_OUT))
        finally:
            silent_remove(GOODVIBES_DAT, disable=DISABLE_REMOVE)
            silent_remove(AE_OUT, disable=DISABLE_REMOVE)
            pass

    def testBiomol(self):
        test_input = ["-l", BI_LIST, "-d", SUB_DATA_DIR]
        try:
            main(test_input)
            self.assertFalse(diff_lines(AE_OUT, GOOD_AE_BI_OUT))
        finally:
            silent_remove(GOODVIBES_DAT, disable=DISABLE_REMOVE)
            silent_remove(AE_OUT, disable=DISABLE_REMOVE)
            pass
