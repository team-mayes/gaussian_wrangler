import unittest
import os
from gaussian_wrangler.gausscom2pdb import main
from common_wrangler.common import diff_lines, silent_remove, capture_stdout, capture_stderr
import logging

# logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger(__name__)
DISABLE_REMOVE = logger.isEnabledFor(logging.DEBUG)

__author__ = 'hmayes'

TEST_DIR = os.path.dirname(__file__)
MAIN_DIR = os.path.dirname(TEST_DIR)
DATA_DIR = os.path.join(os.path.dirname(__file__), 'test_data')
SUB_DATA_DIR = os.path.join(DATA_DIR, 'gausscom2pdb')
DEF_INI = os.path.join(SUB_DATA_DIR, 'gausscom2pdb.ini')

NO_TPL_INI = os.path.join(SUB_DATA_DIR, 'gausscom2pdb_no_tpl.ini')

MISSING_PDB_INI = os.path.join(SUB_DATA_DIR, 'gausscom2pdb_missing_tpl.ini')
BAD_ATOM_INI = os.path.join(SUB_DATA_DIR, 'gausscom2pdb_bad_atom.ini')
EMPTY_LIST_INI = os.path.join(SUB_DATA_DIR, 'gausscom2pdb_empty_list.ini')

# noinspection PyUnresolvedReferences
PDB_OUT = os.path.join(SUB_DATA_DIR, 'pet_cp1_def2_end.pdb')
GOOD_PDB_OUT = os.path.join(SUB_DATA_DIR, 'pet_cp1_def2_end_good.pdb')

GOOD_NO_TPL_OUT = os.path.join(SUB_DATA_DIR, 'pet_cp1_def2_end_no_tpl_good.pdb')


class TestGausscom2pdbNoOut(unittest.TestCase):
    # These all test failure cases
    def testNoArgs(self):
        with capture_stderr(main, []) as output:
            self.assertTrue("WARNING:  Problems reading file: Could not read file" in output)
        with capture_stdout(main, []) as output:
            self.assertTrue("optional arguments" in output)

    # def testTypoIni(self):
    #     with capture_stderr(main, ["-c", TYPO_INI]) as output:
    #         self.assertTrue("Unexpected key" in output)
    #
    # def testMissReqKeyIni(self):
    #     with capture_stderr(main, ["-c", MISS_INI]) as output:
    #         self.assertTrue("Missing config val for key 'pdb_tpl_file'" in output)
    #
    def testNoFiles(self):
        with capture_stderr(main, ["-c", EMPTY_LIST_INI]) as output:
            self.assertTrue("No files to process" in output)

    def testNoFilesInList(self):
        with capture_stderr(main, ["-c", MISSING_PDB_INI]) as output:
            self.assertTrue("Problems reading file" in output)

    def testNotIni(self):
        # gracefully fail if give the wrong file to the -c option
        test_input = ["-c", GOOD_PDB_OUT]
        if logger.isEnabledFor(logging.DEBUG):
            main(test_input)
        with capture_stderr(main, test_input) as output:
            self.assertTrue("WARNING:  File contains no section headers" in output)

    def testHelp(self):
        test_input = ['-h']
        if logger.isEnabledFor(logging.DEBUG):
            main(test_input)
        with capture_stderr(main, test_input) as output:
            self.assertFalse(output)
        with capture_stdout(main, test_input) as output:
            self.assertTrue("optional arguments" in output)


class TestGaussCom2pdb(unittest.TestCase):
    # These test/demonstrate different options
    def testDefIni(self):
        test_input = ["-c", DEF_INI]
        if logger.isEnabledFor(logging.DEBUG):
            main(test_input)
        try:
            main(test_input)
            self.assertFalse(diff_lines(PDB_OUT, GOOD_PDB_OUT))
        finally:
            # silent_remove(PDB_OUT, disable=DISABLE_REMOVE)
            pass

    def testNoTplIni(self):
        test_input = ["-c", NO_TPL_INI]
        if logger.isEnabledFor(logging.DEBUG):
            main(test_input)
        try:
            main(test_input)
            self.assertFalse(diff_lines(PDB_OUT, GOOD_NO_TPL_OUT))
        finally:
            # silent_remove(PDB_TPL_OUT)
            silent_remove(PDB_OUT, disable=DISABLE_REMOVE)

    def testBadAtomIni(self):
        test_input = ["-c", BAD_ATOM_INI]
        if logger.isEnabledFor(logging.DEBUG):
            main(test_input)
        with capture_stderr(main, test_input) as output:
            self.assertTrue("Atom types do not match" in output)
        silent_remove(PDB_OUT, disable=DISABLE_REMOVE)
