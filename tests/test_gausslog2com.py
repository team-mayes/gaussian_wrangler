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

CP_COM_TPL = os.path.join(SUB_DATA_DIR, 'cp.tpl')
LOG_LIST = os.path.join(SUB_DATA_DIR, 'log_list.txt')
COM1_OUT = os.path.join(SUB_DATA_DIR, 'pet_mono_843_tzvp_cp.com')
GOOD_COM1_OUT = os.path.join(SUB_DATA_DIR, 'pet_mono_843_tzvp_cp_good.com')
COM2_OUT = os.path.join(SUB_DATA_DIR, 'pet_mono_901_tzvp_cp.com')
GOOD_COM2_OUT = os.path.join(SUB_DATA_DIR, 'pet_mono_901_tzvp_cp_good.com')

LOG_FILE = os.path.join(SUB_DATA_DIR, 'pet_mono_1_tzvp.log')
COM3_OUT = os.path.join(SUB_DATA_DIR, 'pet_mono_1_tzvp_cp.com')
GOOD_COM3_OUT = os.path.join(SUB_DATA_DIR, 'pet_mono_1_tzvp_cp_good.com')

LOG_LOW_E_FILE = os.path.join(SUB_DATA_DIR, 'nylon66_dimer_49_69_f2_fb.log')
LOW_E_TPL = os.path.join(SUB_DATA_DIR, 'fresh_fb.tpl')
COM_LOW_E_OUT = os.path.join(SUB_DATA_DIR, 'nylon66_dimer_49_69_f2_fb_fresh_fb.com')
GOOD_COM_LOW_E_FILE = os.path.join(SUB_DATA_DIR, 'nylon66_dimer_49_69_f2_fb_fresh_fb_good.com')

TYPE_MATCH_TPL = os.path.join(SUB_DATA_DIR, 's12but_17_84_f2.com')
TYPE_MATCH_LOG = os.path.join(SUB_DATA_DIR, 's12but_17_84_f2.log')
TYPE_MATCH_OUT = os.path.join(SUB_DATA_DIR, 's12but_17_84_f2_s12but_17_84_f2.com')
GOOD_TYPE_MATCH_OUT = os.path.join(SUB_DATA_DIR, 's12but_17_84_f2_good.com')

TYPE_NUM_MATCH_TPL = os.path.join(SUB_DATA_DIR, 's12but_17_84_f2_tpl.com')
TYPE_NUM_MATCH_OUT = os.path.join(MAIN_DIR, 's12but_17_84_f2_t.com')
GOOD_TYPE_NUM_MATCH_OUT = os.path.join(SUB_DATA_DIR, 's12but_17_84_f2_tpl_good.com')


class TestGausslog2comNoOut(unittest.TestCase):
    # These all test failure cases
    def testNoArgs(self):
        test_input = []
        with capture_stderr(main, test_input) as output:
            self.assertTrue("No files to process" in output)

    def testHelp(self):
        test_input = ['-h']
        if logger.isEnabledFor(logging.DEBUG):
            main(test_input)
        with capture_stderr(main, test_input) as output:
            self.assertFalse(output)
        with capture_stdout(main, test_input) as output:
            self.assertTrue("optional arguments" in output)

    def testMissingTpl(self):
        test_input = ["-f", LOG_FILE, "-c"]
        # main(test_input)
        if logger.isEnabledFor(logging.DEBUG):
            main(test_input)
        with capture_stderr(main, test_input) as output:
            self.assertTrue("No template file" in output)

    def testMissingFile(self):
        test_input = ["-t", "ghost.tpl", "-f", LOG_FILE, "-c"]
        if logger.isEnabledFor(logging.DEBUG):
            main(test_input)
        with capture_stderr(main, test_input) as output:
            self.assertTrue("Problems reading file" in output)

    def testMisMatchLogCom(self):
        test_input = ["-t", TYPE_MATCH_TPL, "-f", LOG_LOW_E_FILE]
        # main(test_input)
        with capture_stderr(main, test_input) as output:
            self.assertTrue("has atom type" in output)

    def testNotTpl(self):
        test_input = ["-t", LOG_LOW_E_FILE, "-f", LOG_LOW_E_FILE]
        # main(test_input)
        with capture_stderr(main, test_input) as output:
            self.assertTrue("Problems reading data" in output)

    def testNotLog(self):
        test_input = ["-t", TYPE_MATCH_TPL, "-f", TYPE_MATCH_TPL, "-e"]
        # main(test_input)
        with capture_stderr(main, test_input) as output:
            self.assertTrue("log file has coordinates to use" in output)


class TestGausslog2com(unittest.TestCase):
    # These test/demonstrate different options
    def testCPTpl(self):
        test_input = ["-t", CP_COM_TPL, "-l", LOG_LIST, "-c"]
        try:
            main(test_input)
            self.assertFalse(diff_lines(COM1_OUT, GOOD_COM1_OUT))
            self.assertFalse(diff_lines(COM2_OUT, GOOD_COM2_OUT))
        finally:
            silent_remove(COM1_OUT, disable=DISABLE_REMOVE)
            silent_remove(COM2_OUT, disable=DISABLE_REMOVE)
            pass

    def testFileCPTpl(self):
        test_input = ["-t", CP_COM_TPL, "-f", LOG_FILE, "-c"]
        try:
            main(test_input)
            self.assertFalse(diff_lines(COM3_OUT, GOOD_COM3_OUT))
        finally:
            silent_remove(COM3_OUT, disable=DISABLE_REMOVE)
            pass

    def testLowEn(self):
        test_input = ["-t", LOW_E_TPL, "-f", LOG_LOW_E_FILE, "-e"]
        try:
            main(test_input)
            self.assertFalse(diff_lines(COM_LOW_E_OUT, GOOD_COM_LOW_E_FILE))
        finally:
            silent_remove(COM_LOW_E_OUT, disable=DISABLE_REMOVE)
            pass

    def testTemplateTypeMatching(self):
        test_input = ["-t", TYPE_MATCH_TPL, "-f", TYPE_MATCH_LOG, "-e"]
        try:
            main(test_input)
            self.assertFalse(diff_lines(TYPE_MATCH_OUT, GOOD_TYPE_MATCH_OUT))
        finally:
            silent_remove(TYPE_MATCH_OUT, disable=DISABLE_REMOVE)
            pass

    def testTemplateNumberMatching(self):
        test_input = ["-t", TYPE_NUM_MATCH_TPL, "-f", TYPE_MATCH_LOG, "-e", "-o", "s12but_17_84_f2_t.com"]
        try:
            main(test_input)
            self.assertFalse(diff_lines(TYPE_NUM_MATCH_OUT, GOOD_TYPE_NUM_MATCH_OUT))
        finally:
            silent_remove(TYPE_NUM_MATCH_OUT, disable=DISABLE_REMOVE)
            pass
