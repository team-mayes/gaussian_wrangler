import unittest
import os
from gaussian_wrangler.goodvibes_hm import main
from common_wrangler.common import (capture_stdout, capture_stderr)
import logging

# logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger(__name__)
DISABLE_REMOVE = logger.isEnabledFor(logging.DEBUG)

__author__ = 'hmayes'

TEST_DIR = os.path.dirname(__file__)
MAIN_DIR = os.path.dirname(TEST_DIR)
DATA_DIR = os.path.join(os.path.dirname(__file__), 'test_data')
SUB_DATA_DIR = os.path.join(DATA_DIR, 'goodvibes_helper')

GOODVIBES_DAT = os.path.abspath(os.path.join(TEST_DIR, '..', 'Goodvibes_output.dat'))

TEST_LOG1 = os.path.join(SUB_DATA_DIR, 'ethygly2_tzvp.log')
TEST_LOG2 = os.path.join(SUB_DATA_DIR, 'tpaegh1ats_ts.log')
TEST_LOG3 = os.path.join(SUB_DATA_DIR, 'ts3b_ircr_opt_gas.log')
TEST_LOG4 = os.path.join(SUB_DATA_DIR, 'co_gas.log')
TEST_LOG5 = os.path.join(SUB_DATA_DIR, 'lmethyllactate_1_8_cp.log')

INCOMPLETE_LOG = os.path.join(SUB_DATA_DIR, 'ipah_d_incomplete.log')
FAILED_LOG = os.path.join(SUB_DATA_DIR, 'co_fail_gas.log')


class TestGoodVibesNoOut(unittest.TestCase):
    # These all test failure cases
    def testNoArgs(self):
        test_input = []
        # main(test_input)
        with capture_stderr(main, test_input) as output:
            self.assertTrue("No file" in output)

    def testHelp(self):
        test_input = ['-h']
        # main(test_input)
        if logger.isEnabledFor(logging.DEBUG):
            main(test_input)
        with capture_stderr(main, test_input) as output:
            self.assertFalse(output)
        with capture_stdout(main, test_input) as output:
            self.assertTrue("optional arguments" in output)

    def testMissingValue(self):
        test_input = [TEST_LOG1, "-c"]
        # main(test_input)
        if logger.isEnabledFor(logging.DEBUG):
            main(test_input)
        with capture_stderr(main, test_input) as output:
            self.assertTrue("expected one argument" in output)

    def testNoSuchFile(self):
        test_input = ["ghost.log"]
        # main(test_input)
        if logger.isEnabledFor(logging.DEBUG):
            main(test_input)
        with capture_stderr(main, test_input) as output:
            self.assertTrue("No file names" in output)

    def testIncompleteFile(self):
        test_input = [INCOMPLETE_LOG]
        # main(test_input)
        with capture_stderr(main, test_input) as output:
            self.assertTrue('may not have terminated' in output)

    def testFailedLog(self):
        test_input = [FAILED_LOG]
        # main(test_input)
        with capture_stderr(main, test_input) as output:
            self.assertTrue('Error termination' in output)

    def testCP(self):
        test_input = [TEST_LOG5, "--cpu"]
        # main(test_input)
        with capture_stderr(main, test_input) as output:
            self.assertTrue("Could not find frequency" in output)
        with capture_stdout(main, test_input) as output:
            self.assertTrue("0 days  0 hrs 27 mins 11 secs" in output)

    def testCPWithSPC(self):
        test_input = [TEST_LOG5, "--spc"]
        main(test_input)
        with capture_stderr(main, test_input) as output:
            self.assertTrue("SPC calculation" in output)


class TestGoodVibesHM(unittest.TestCase):
    # These test/demonstrate different options
    def testFileName(self):
        test_input = [TEST_LOG1]
        # main(test_input)
        with capture_stdout(main, test_input) as output:
            self.assertTrue("-230.257454   0.083435   -230.167674   0.033971   0.033946   -230.201645   "
                            "-230.201620" in output)

    def testTempVib(self):
        test_input = [TEST_LOG1, "-t", "788.15", "-v", "0.984"]
        # main(test_input)
        with capture_stdout(main, test_input) as output:
            self.assertTrue("-230.257454   0.084552   -230.144250   0.122864   0.122499   -230.267114   "
                            "-230.266749" in output)

    def testTempVibFreq(self):
        test_input = [TEST_LOG3, "-t", "788.15", "-v", "0.984", "-f", '30']
        main(test_input)
        with capture_stdout(main, test_input) as output:
            self.assertTrue("-460.884947   0.180092   -460.642788   0.204516   0.204169   -460.847303   "
                            "-460.846956" in output)

    def testLinearMolecule(self):
        test_input = [TEST_LOG4, "--qs", "Truhlar"]
        # main(test_input)
        with capture_stdout(main, test_input) as output:
            self.assertTrue("-113.322295   0.005041   -113.313948   0.022414   0.022414   -113.336362   "
                            "-113.336362" in output)

    def testTempVibConc(self):
        test_input = [TEST_LOG1, "-t", "788.15", "-v", "0.984", "-c", "1"]
        # main(test_input)
        with capture_stdout(main, test_input) as output:
            self.assertTrue("-230.257454   0.084552   -230.144250   0.112457   0.112093   -230.256708   "
                            "-230.256343" in output)

    def testTempRangeVib(self):
        test_input = [TEST_LOG1, "-t", "788.15", "-v", "0.984", "--ti", "688.15,888.15,25"]
        # main(test_input)
        with capture_stdout(main, test_input) as output:
            self.assertTrue("688.15              -230.149803   0.102094   0.101814   -230.251897   "
                            "-230.251617" in output)
            self.assertTrue("788.15              -230.144250   0.122864   0.122499   -230.267114   "
                            "-230.266749" in output)
            self.assertTrue("888.15              -230.138332   0.144728   0.144273   -230.283060   "
                            "-230.282605" in output)

    def testTempRangeVibConc(self):
        test_input = [TEST_LOG1, "-t", "788.15", "-v", "0.984", "-c", "1", "--ti", "688.15,888.15,25"]
        # main(test_input)
        with capture_stdout(main, test_input) as output:
            self.assertTrue("688.15              -230.149803   0.093304   0.093023   -230.243107   "
                            "-230.242827" in output)
            self.assertTrue("788.15              -230.144250   0.112457   0.112093   -230.256708   "
                            "-230.256343" in output)
            self.assertTrue("888.15              -230.138332   0.132665   0.132210   -230.270997   "
                            "-230.270543" in output)

    def testTempRangeVibQ(self):
        test_input = [TEST_LOG1, "-t", "788.15", "-v", "0.984", "--ti", "688.15,888.15,25", "-q"]
        # main(test_input)
        with capture_stdout(main, test_input) as output:
            self.assertTrue("688.15              -230.149803   -230.150097   0.102094   0.101814   -230.251897   "
                            "-230.251911" in output)
            self.assertTrue("788.15              -230.144250   -230.144585   0.122864   0.122499   -230.267114   "
                            "-230.267084" in output)
            self.assertTrue("888.15              -230.138332   -230.138708   0.144728   0.144273   -230.283060   "
                            "-230.282981" in output)

    def testTempRangeVibConcQ(self):
        test_input = [TEST_LOG1, "-t", "788.15", "-v", "0.984", "-c", "1", "--ti", "688.15,888.15,25", "-q"]
        # main(test_input)
        with capture_stdout(main, test_input) as output:
            self.assertTrue("688.15              -230.149803   -230.150097   0.093304   0.093023   -230.243107   "
                            "-230.243121" in output)
            self.assertTrue("788.15              -230.144250   -230.144585   0.112457   0.112093   -230.256708   "
                            "-230.256677" in output)
            self.assertTrue("888.15              -230.138332   -230.138708   0.132665   0.132210   -230.270997   "
                            "-230.270918" in output)

    def testTempRangeVibConcQAltInput(self):
        test_input = [TEST_LOG2, "-t", "788.15", "-v", "0.984", "-c", "1", "--ti", "688.15,888.15,25", "-q"]
        # main(test_input)
        with capture_stdout(main, test_input) as output:
            self.assertTrue("688.15              -839.735116   -839.741128   0.212003   0.200263   -839.947119   "
                            "-839.941391" in output)
            self.assertTrue("788.15              -839.717088   -839.723965   0.261737   0.247365   -839.978825   "
                            "-839.971331" in output)
            self.assertTrue("888.15              -839.698000   -839.705743   0.314850   0.297737   -840.012851   "
                            "-840.003480" in output)

    def testTempRangeVibConcQCheck(self):
        test_input = [TEST_LOG1, TEST_LOG2, "-t", "788.15", "-v", "0.984", "-c", "1", "--ti", "688.15,888.15,25", "-q",
                      "--check"]
        good_output = "o  Using Gaussian 16 Revision B.01 in all calculations.\n" \
                      "o  Using M062X/def2TZVP in all calculations.\n" \
                      "o  Using scrf=(solvent=1,2-ethanediol,read) in all calculations.\n" \
                      "o  Using a standard concentration of 1 M for solvent phase.\n" \
                      "o  No duplicates or enantiomers found\n" \
                      "-  No linear molecules found.\n" \
                      "-  No empirical dispersion detected in any of the calculations."

        good_error = "WARNING:  Different charge and multiplicity found: \n" \
                     "WARNING:          0 1: tests/test_data/goodvibes_helper/ethygly2_tzvp.log\n" \
                     "WARNING:          1 1: tests/test_data/goodvibes_helper/tpaegh1ats_ts.log"
        # main(test_input)
        with capture_stdout(main, test_input) as output:
            self.assertTrue(good_output in output)
        with capture_stderr(main, test_input) as output:
            self.assertTrue(good_error in output)
