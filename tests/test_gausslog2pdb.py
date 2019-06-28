import unittest
import os
from nrel_tools.gausslog2pdb import main
from nrel_tools.common import diff_lines, silent_remove, capture_stdout, capture_stderr
import logging

# logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger(__name__)
DISABLE_REMOVE = logger.isEnabledFor(logging.DEBUG)

__author__ = 'hmayes'

TEST_DIR = os.path.dirname(__file__)
MAIN_DIR = os.path.dirname(TEST_DIR)
DATA_DIR = os.path.join(os.path.dirname(__file__), 'test_data')
SUB_DATA_DIR = os.path.join(DATA_DIR, 'gausslog2pdb')

DEF_INI = os.path.join(SUB_DATA_DIR, 'gausslog2pdb.ini')
# noinspection PyUnresolvedReferences
PDB_OUT = os.path.join(SUB_DATA_DIR, 'pet_mono_f1hs_1.pdb')
GOOD_PDB_ALL_OUT = os.path.join(SUB_DATA_DIR, 'pet_mono_f1hs_1_all_good.pdb')

LAST_INI = os.path.join(SUB_DATA_DIR, 'gausslog2pdb_last.ini')
GOOD_PDB_LAST_OUT = os.path.join(SUB_DATA_DIR, 'pet_mono_f1hs_1_last_good.pdb')

MULT_INI = os.path.join(SUB_DATA_DIR, 'gausslog2pdb_mult.ini')
PDB2_OUT = os.path.join(SUB_DATA_DIR, 'pet_mono_f1hs_2.pdb')
GOOD_PDB2_LAST_OUT = os.path.join(SUB_DATA_DIR, 'pet_mono_f1hs_2_last_good.pdb')

COMB_INI = os.path.join(SUB_DATA_DIR, 'gausslog2pdb_comb.ini')
PDB12_OUT = os.path.join(SUB_DATA_DIR, 'pet_mono_f1hs_comb.pdb')
GOOD_PDB12_LAST_OUT = os.path.join(SUB_DATA_DIR, 'pet_mono_f1hs_comb_good.pdb')


class Testgausslog2pdbNoOut(unittest.TestCase):
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


class Testgausslog2pdb(unittest.TestCase):
    # These test/demonstrate different options
    def testDefIni(self):
        test_input = ["-c", DEF_INI]
        if logger.isEnabledFor(logging.DEBUG):
            main(test_input)
        try:
            main(test_input)
            self.assertFalse(diff_lines(PDB_OUT, GOOD_PDB_ALL_OUT))
        finally:
            silent_remove(PDB_OUT, disable=DISABLE_REMOVE)

    def testLastIni(self):
        test_input = ["-c", LAST_INI]
        if logger.isEnabledFor(logging.DEBUG):
            main(test_input)
        try:
            main(test_input)
            self.assertFalse(diff_lines(PDB_OUT, GOOD_PDB_LAST_OUT))
        finally:
            silent_remove(PDB_OUT, disable=DISABLE_REMOVE)

    def testMultIni(self):
        test_input = ["-c", MULT_INI]
        if logger.isEnabledFor(logging.DEBUG):
            main(test_input)
        try:
            main(test_input)
            self.assertFalse(diff_lines(PDB_OUT, GOOD_PDB_LAST_OUT))
            self.assertFalse(diff_lines(PDB2_OUT, GOOD_PDB2_LAST_OUT))
        finally:
            silent_remove(PDB_OUT, disable=DISABLE_REMOVE)
            silent_remove(PDB2_OUT, disable=DISABLE_REMOVE)
            pass

    def testCombIni(self):
        test_input = ["-c", COMB_INI]
        if logger.isEnabledFor(logging.DEBUG):
            main(test_input)
        try:
            main(test_input)
            self.assertFalse(diff_lines(PDB12_OUT, GOOD_PDB12_LAST_OUT))
        finally:
            silent_remove(PDB12_OUT, disable=DISABLE_REMOVE)
            pass

