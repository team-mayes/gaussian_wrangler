import unittest
import os
from gaussian_wrangler.smi2gausscom import main
from common_wrangler.common import diff_lines, silent_remove, capture_stdout, capture_stderr
import logging

# logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger(__name__)
DISABLE_REMOVE = logger.isEnabledFor(logging.DEBUG)

__author__ = 'hmayes'

TEST_DIR = os.path.dirname(__file__)
MAIN_DIR = os.path.dirname(TEST_DIR)
DATA_DIR = os.path.join(os.path.dirname(__file__), 'test_data')
SUB_DATA_DIR = os.path.join(DATA_DIR, 'smi2gausscom')
TEMP_DIR = os.path.join(SUB_DATA_DIR, 'temp')
GAU_TPL = os.path.join(SUB_DATA_DIR, "gau.tpl")


class TestMainNoOutput(unittest.TestCase):
    # These all test failure cases
    def testNoArgs(self):
        test_input = []
        # main(test_input)
        with capture_stderr(main, test_input) as output:
            self.assertTrue("required" in output)
        with capture_stdout(main, test_input) as output:
            self.assertTrue("optional arguments" in output)

    def testHelp(self):
        test_input = ['-h']
        # main(test_input)
        with capture_stderr(main, test_input) as output:
            self.assertFalse(output)
        with capture_stdout(main, test_input) as output:
            self.assertTrue("optional arguments" in output)

    def testNoSMI(self):
        test_input = ['-t', GAU_TPL]
        # main(test_input)
        with capture_stderr(main, test_input) as output:
            self.assertTrue("SMILES " in output)
        with capture_stdout(main, test_input) as output:
            self.assertTrue("optional arguments" in output)

    def testTplMissingField(self):
        missing_field_tpl = os.path.join(SUB_DATA_DIR, "gau_missing_field.tpl")
        test_input = ["-s", "CCCO", "-t", missing_field_tpl, "-o", TEMP_DIR]
        # main(test_input)
        with capture_stderr(main, test_input) as output:
            self.assertTrue("required string " in output)

    def testNoSuchTpl(self):
        test_input = ["-s", "CCCO", "-t", "ghost", "-o", TEMP_DIR]
        # main(test_input)
        with capture_stderr(main, test_input) as output:
            self.assertTrue("Could not locate" in output)

    def testNoSuchListFile(self):
        test_input = ["-l", "ghost", "-t", GAU_TPL, "-o", TEMP_DIR]
        # main(test_input)
        with capture_stderr(main, test_input) as output:
            self.assertTrue("No such" in output)

    def testBadSmi(self):
        test_input = ["-s", "CCC0,CCCO", "-t", GAU_TPL, "-o", TEMP_DIR]
        an_output = os.path.join(TEMP_DIR, "cid_1031_0.com")
        try:
            silent_remove(TEMP_DIR, dir_with_files=True)
            # main(test_input)
            with capture_stderr(main, test_input) as output:
                self.assertTrue("Skipping" in output)
            self.assertTrue(os.path.isfile(an_output))
        finally:
            silent_remove(TEMP_DIR, dir_with_files=True)
            pass


class TestMain(unittest.TestCase):
    # These test/demonstrate different options
    def testOneSmi(self):
        test_input = ["-s", "CCCO", "-t", GAU_TPL, "-o", TEMP_DIR]
        try:
            silent_remove(TEMP_DIR, dir_with_files=True)
            main(test_input)
            num_created_files = len(os.listdir(TEMP_DIR))
            self.assertEqual(7, num_created_files)
            for conf_id in range(3):
                good_out_fname = os.path.join(SUB_DATA_DIR, f"cid_1031_{conf_id}_good.com")
                out_fname = os.path.join(TEMP_DIR, f"cid_1031_{conf_id}.com")
                self.assertFalse(diff_lines(good_out_fname, out_fname))
        finally:
            # silent_remove(TEMP_DIR, dir_with_files=True)
            pass

    def testSmiWithSpecialCharacters(self):
        smi = "CC(C)O,C=C1C=CC(=O)C=C1,C([C@@H]1[C@H]([C@@H]([C@H]([C@@H](O1)O)O)O)O)O," \
              "C1=CC=C(C2COC3COCC32)C=C1"
        test_input = ["-s", smi, "-t", GAU_TPL, "-o", TEMP_DIR, "-m", "10"]
        try:
            silent_remove(TEMP_DIR, dir_with_files=True)
            main(test_input)
            num_created_files = len(os.listdir(TEMP_DIR))
            self.assertEqual(24, num_created_files)
            for base_name in ["cid_136328", "cid_64689", "cid_3776", "C1eCCeCpC2COC3COCC32qCeC1"]:
                self.assertTrue(os.path.isfile(os.path.join(TEMP_DIR, f"{base_name}_0.com")))
        finally:
            silent_remove(TEMP_DIR, dir_with_files=True)
            pass

    def testSmiList(self):
        # smoke test
        list_fname = os.path.join(SUB_DATA_DIR, "smi_list.txt")
        test_input = ["-l", list_fname, "-t", GAU_TPL, "-o", TEMP_DIR, "-m", "10"]
        try:
            silent_remove(TEMP_DIR, dir_with_files=True)
            main(test_input)
            num_created_files = len(os.listdir(TEMP_DIR))
            self.assertEqual(60, num_created_files)
        finally:
            silent_remove(TEMP_DIR, dir_with_files=True)
            pass

    def testSmiListDup(self):
        # check that making a set removes duplicate input
        list_fname = os.path.join(SUB_DATA_DIR, "smi_list_with_dup.txt")
        test_input = ["-l", list_fname, "-t", GAU_TPL, "-o", TEMP_DIR, "-m", "10"]
        try:
            silent_remove(TEMP_DIR, dir_with_files=True)
            # main(test_input)
            with capture_stdout(main, test_input) as output:
                self.assertEqual(output.count("cid_811"), 1)
        finally:
            silent_remove(TEMP_DIR, dir_with_files=True)
            pass
