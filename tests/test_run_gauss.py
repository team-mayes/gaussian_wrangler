import unittest
import os
from nrel_tools.run_gauss import main
from nrel_tools.common import diff_lines, silent_remove, capture_stdout, capture_stderr
import logging

# logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger(__name__)
DISABLE_REMOVE = logger.isEnabledFor(logging.DEBUG)

__author__ = 'hmayes'

TEST_DIR = os.path.dirname(__file__)
PARENT_DIR = os.path.join(TEST_DIR, os.pardir)
MAIN_DIR = os.path.dirname(TEST_DIR)
DATA_DIR = os.path.join(os.path.dirname(__file__), 'test_data')
SUB_DATA_DIR = os.path.join(DATA_DIR, 'run_gauss')

DEF_INI = os.path.join(SUB_DATA_DIR, 'run_gauss_bde.ini')
DEF_SH_OUT = os.path.join(PARENT_DIR, 'ethylrad.sh')
GOOD_SH_OUT = os.path.join(SUB_DATA_DIR, 'good_ethylrad.sh')
DEF_LOG_OUT = os.path.join(PARENT_DIR, 'ethylrad.log')
ONE_JOB_INI = os.path.join(SUB_DATA_DIR, 'run_gauss_bde_one_job.ini')
GOOD_ONE_SH_OUT = os.path.join(SUB_DATA_DIR, 'good_ethylrad_one.sh')

MISSING_TPL_INI = os.path.join(SUB_DATA_DIR, 'run_gauss_missing_tpl.ini')
ONE_NEW_JOB_INI = os.path.join(SUB_DATA_DIR, 'run_gauss_one.ini')


class TestRunGaussBDENoOut(unittest.TestCase):
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

    def testMissingListIni(self):
        test_input = ["tests/test_data/run_gauss/ethylrad", "-c", MISSING_TPL_INI]
        if logger.isEnabledFor(logging.DEBUG):
            main(test_input)
        with capture_stderr(main, test_input) as output:
            self.assertTrue("could not find a template file" in output)


class TestRunGaussBDE(unittest.TestCase):
    # These test/demonstrate different options
    def testDefIni(self):
        test_input = ["tests/test_data/run_gauss/ethylrad", "-c", DEF_INI]
        try:
            main(test_input)
            self.assertFalse(diff_lines(DEF_SH_OUT, GOOD_SH_OUT))
        except IOError:
            pass
        finally:
            silent_remove(DEF_SH_OUT, disable=DISABLE_REMOVE)
            silent_remove(DEF_LOG_OUT, disable=DISABLE_REMOVE)
            pass

    def testOneJobIni(self):
        test_input = ["tests/test_data/run_gauss/ethylrad", "-c", ONE_JOB_INI]
        try:
            main(test_input)
            self.assertFalse(diff_lines(DEF_SH_OUT, GOOD_ONE_SH_OUT))
        finally:
            silent_remove(DEF_SH_OUT, disable=DISABLE_REMOVE)
            silent_remove(DEF_LOG_OUT, disable=DISABLE_REMOVE)
            pass
