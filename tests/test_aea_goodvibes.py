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

UNI_REACT = os.path.join(SUB_DATA_DIR, 'ipaegh1dts_t_ircr_opt.log')
UNI_TS = os.path.join(SUB_DATA_DIR, 'ipaegh1dts.log')
FILE_LIST = os.path.join(SUB_DATA_DIR, 'list.txt')


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


class TestAEaGoodVibes(unittest.TestCase):
    # These test/demonstrate different options
    def testCPTpl(self):
        test_input = [UNI_REACT, UNI_TS, "-l", FILE_LIST]
        try:
            main(test_input)
            # self.assertFalse(diff_lines(COM1_OUT, GOOD_COM1_OUT))
        finally:
            # silent_remove(COM1_OUT, disable=DISABLE_REMOVE)
            # silent_remove(COM2_OUT, disable=DISABLE_REMOVE)
            pass

