import unittest
import os
from nrel_tools.pdbs2gausscoms import main
from nrel_tools.common import diff_lines, silent_remove, capture_stdout, capture_stderr
import logging

# logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger(__name__)
DISABLE_REMOVE = logger.isEnabledFor(logging.DEBUG)

__author__ = 'hmayes'

TEST_DIR = os.path.dirname(__file__)
MAIN_DIR = os.path.dirname(TEST_DIR)
DATA_DIR = os.path.join(os.path.dirname(__file__), 'test_data')
SUB_DATA_DIR = os.path.join(DATA_DIR, 'pdbs2gausscoms')

DEF_INI = os.path.join(SUB_DATA_DIR, 'pdb2gau.ini')
GAU_OUT1 = os.path.join(SUB_DATA_DIR, 'pet_mono_1.com')
GOOD_GAU_OUT1 = os.path.join(SUB_DATA_DIR, 'pet_mono_1_good.com')
GAU_OUT2 = os.path.join(SUB_DATA_DIR, 'pet_mono_2.com')
GOOD_GAU_OUT2 = os.path.join(SUB_DATA_DIR, 'pet_mono_2_good.com')
GAU_OUT3 = os.path.join(SUB_DATA_DIR, 'pet_mono_3.com')
GOOD_GAU_OUT3 = os.path.join(SUB_DATA_DIR, 'pet_mono_3_good.com')

REMOVE_H_INI = os.path.join(SUB_DATA_DIR, 'pdb2gau_h.ini')
REMOVE_H_OUT1 = os.path.join(SUB_DATA_DIR, 'pet_mono_f1hs_1.com')
GOOD_REMOVE_H_OUT1 = os.path.join(SUB_DATA_DIR, 'pet_mono_f1hs_1_good.com')
REMOVE_H_OUT2 = os.path.join(SUB_DATA_DIR, 'pet_mono_f1hs_2.com')
GOOD_REMOVE_H_OUT2 = os.path.join(SUB_DATA_DIR, 'pet_mono_f1hs_2_good.com')
REMOVE_H_OUT3 = os.path.join(SUB_DATA_DIR, 'pet_mono_f1hs_3.com')
GOOD_REMOVE_H_OUT3 = os.path.join(SUB_DATA_DIR, 'pet_mono_f1hs_3_good.com')


class Testpdbs2gausscomsNoOut(unittest.TestCase):
    # These all test failure cases
    def testNoArgs(self):
        with capture_stderr(main, []) as output:
            self.assertTrue("WARNING:  Problems reading file: Could not read file" in output)
        with capture_stdout(main, []) as output:
            self.assertTrue("optional arguments" in output)

    # def testHelp(self):
    #     test_input = ['-h']
    #     if logger.isEnabledFor(logging.DEBUG):
    #         main(test_input)
    #     with capture_stderr(main, test_input) as output:
    #         self.assertFalse(output)
    #     with capture_stdout(main, test_input) as output:
    #         self.assertTrue("optional arguments" in output)


class Testpdbs2gausscoms(unittest.TestCase):
    # These test/demonstrate different options
    def testDefIni(self):
        test_input = ["-c", DEF_INI]
        if logger.isEnabledFor(logging.DEBUG):
            main(test_input)
        try:
            main(test_input)
            self.assertFalse(diff_lines(GAU_OUT1, GOOD_GAU_OUT1))
            self.assertFalse(diff_lines(GAU_OUT2, GOOD_GAU_OUT2))
            self.assertFalse(diff_lines(GAU_OUT3, GOOD_GAU_OUT3))
        finally:
            silent_remove(GAU_OUT1, disable=DISABLE_REMOVE)
            silent_remove(GAU_OUT2, disable=DISABLE_REMOVE)
            silent_remove(GAU_OUT3, disable=DISABLE_REMOVE)

    def testRemoveHIni(self):
        test_input = ["-c", REMOVE_H_INI]
        if logger.isEnabledFor(logging.DEBUG):
            main(test_input)
        try:
            main(test_input)
            self.assertFalse(diff_lines(REMOVE_H_OUT1, GOOD_REMOVE_H_OUT1))
            self.assertFalse(diff_lines(REMOVE_H_OUT2, GOOD_REMOVE_H_OUT2))
            self.assertFalse(diff_lines(REMOVE_H_OUT3, GOOD_REMOVE_H_OUT3))
        finally:
            silent_remove(REMOVE_H_OUT1, disable=DISABLE_REMOVE)
            silent_remove(REMOVE_H_OUT2, disable=DISABLE_REMOVE)
            silent_remove(REMOVE_H_OUT3, disable=DISABLE_REMOVE)

