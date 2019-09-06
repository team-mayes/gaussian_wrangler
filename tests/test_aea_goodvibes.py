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
FILE_LIST_MISMATCH_STOICH = os.path.join(SUB_DATA_DIR, 'list_mismatch_stoich.txt')

AE_OUT = os.path.join(SUB_DATA_DIR, 'aea_out.csv')
GOOD_AE_OUT = os.path.join(SUB_DATA_DIR, 'aea_out_good.csv')

BI_LIST = os.path.join(SUB_DATA_DIR, 'list_bimolec.txt')
GOOD_AE_BI_OUT = os.path.join(SUB_DATA_DIR, 'aea_out_bi_good.csv')
BI_VIBES_OUT1 = os.path.join(SUB_DATA_DIR, 'ethygly2_tzvp_vibes.dat')
BI_VIBES_OUT2 = os.path.join(SUB_DATA_DIR, 'pdc2_h_vibes.dat')
BI_VIBES_OUT3 = os.path.join(SUB_DATA_DIR, 'pdc2_eghtsct_vibes.dat')

TI_LIST = os.path.join(SUB_DATA_DIR, 'list_ti.txt')
GOOD_AE_TI_OUT = os.path.join(SUB_DATA_DIR, 'aea_out_ti_good.csv')
TI_VIBES_OUT = os.path.join(SUB_DATA_DIR, 'aea_out_vibes.dat')


class TestAEaGoodVibesNoOut(unittest.TestCase):
    # These all test failure cases
    def testNoArgs(self):
        test_input = []
        # main(test_input)
        with capture_stderr(main, test_input) as output:
            self.assertTrue("No files" in output)
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
    def testNoSuchFile(self):
        test_input = ["ghost.log"]
        try:
            # main(test_input)
            with capture_stderr(main, test_input) as output:
                self.assertTrue("Problems reading file" in output)
        finally:
            silent_remove(GOODVIBES_DAT, disable=DISABLE_REMOVE)
            silent_remove(GOODVIBES_CSV, disable=DISABLE_REMOVE)
            pass

    def testMisMatchStoich(self):
        test_input = ["-l", FILE_LIST_MISMATCH_STOICH]
        try:
            # main(test_input)
            with capture_stderr(main, test_input) as output:
                self.assertTrue("Check stoichiometries" in output)
        finally:
            silent_remove(GOODVIBES_DAT, disable=DISABLE_REMOVE)
            silent_remove(GOODVIBES_CSV, disable=DISABLE_REMOVE)
            pass

    def testMixMatchSolvent(self):
        test_input = ["-l", FILE_LIST_MISMATCH_SOLV]
        try:
            # main(test_input)
            with capture_stderr(main, test_input) as output:
                self.assertTrue("Different solvents" in output)
        finally:
            silent_remove(GOODVIBES_DAT, disable=DISABLE_REMOVE)
            silent_remove(GOODVIBES_CSV, disable=DISABLE_REMOVE)
            pass

    def testMixMatchTheory(self):
        test_input = [UNI_REACT, UNI_TS]
        try:
            # main(test_input)
            with capture_stderr(main, test_input) as output:
                self.assertTrue("Different basis sets" in output)
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
        # checks a bimolecular reaction and also saving GoodVibes output for each file
        test_input = ["-l", BI_LIST, "-d", SUB_DATA_DIR, "-s"]
        # make sure files not left from a previous run
        for fname in [BI_VIBES_OUT1, BI_VIBES_OUT2, BI_VIBES_OUT3]:
            silent_remove(fname)
        try:
            main(test_input)
            self.assertFalse(diff_lines(AE_OUT, GOOD_AE_BI_OUT))
            for fname in [BI_VIBES_OUT1, BI_VIBES_OUT2, BI_VIBES_OUT3]:
                self.assertTrue(os.path.exists(fname))
        finally:
            for fname in [BI_VIBES_OUT1, BI_VIBES_OUT2, BI_VIBES_OUT3]:
                silent_remove(fname, disable=DISABLE_REMOVE)
            silent_remove(GOODVIBES_DAT, disable=DISABLE_REMOVE)
            silent_remove(AE_OUT, disable=DISABLE_REMOVE)
            pass

    def testTi(self):
        # check handles it when not all atoms in are in all molecules
        # also checks saving GoodVibes output together
        test_input = ["-l", TI_LIST, "-d", SUB_DATA_DIR, "-t"]
        silent_remove(TI_VIBES_OUT)
        try:
            main(test_input)
            self.assertTrue(os.path.exists(TI_VIBES_OUT))
        finally:
            silent_remove(GOODVIBES_DAT, disable=DISABLE_REMOVE)
            silent_remove(AE_OUT, disable=DISABLE_REMOVE)
            silent_remove(TI_VIBES_OUT, disable=DISABLE_REMOVE)
            pass
