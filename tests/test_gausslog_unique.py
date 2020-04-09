import unittest
import os
from common_wrangler.common import capture_stdout, capture_stderr, DIHES
from gaussian_wrangler.gausslog_unique import main, compare_gausslog_info, print_results
from gaussian_wrangler.gw_common import process_gausslog_file, CONVERG_ERR, TS, CONVERG
import logging

# logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger(__name__)
DISABLE_REMOVE = logger.isEnabledFor(logging.DEBUG)

__author__ = 'hmayes'

TEST_DIR = os.path.dirname(__file__)
MAIN_DIR = os.path.dirname(TEST_DIR)
DATA_DIR = os.path.join(os.path.dirname(__file__), 'test_data')
SUB_DATA_DIR = os.path.join(DATA_DIR, 'gausslog_unique')

LOG_LIST = os.path.join(SUB_DATA_DIR, 'list.txt')

HEADER = 'File,Convergence,Energy,Enthalpy\n'
ALPHA_FIRST = 'lme2acetoxprpnt_ts3_ircf_opt.log,7.5369,-535.578434,-535.403311\n'
MISSING_FREQ = 'lme2acetoxprpnt_ts3_ircf_opt_no_freq.log,0.9367,-535.578434,None\n'
ENERGY_FIRST = 'lme2acetoxypropionate_25_t.log,0.1303,-535.578443,-535.403293\n'

EMPTY_LIST = os.path.join(SUB_DATA_DIR, 'empty_list.txt')
LIST_NO_FREQ = os.path.join(SUB_DATA_DIR, 'list_no_freq.txt')
MISSING_FILE_LIST = os.path.join(SUB_DATA_DIR, 'list_with_missing_files.txt')
TWO_MOL_LIST = os.path.join(SUB_DATA_DIR, 'list_two_molecules.txt')
TWO_MORE_MOL_LIST = os.path.join(SUB_DATA_DIR, 'list_two_more_molecules.txt')
SIMILAR_LIST = os.path.join(SUB_DATA_DIR, 'list_similar_molecules.txt')
CALCALL_LIST = os.path.join(SUB_DATA_DIR, 'list_calcall.txt')


class TestGausslogUniqueNoOut(unittest.TestCase):
    def testHelp(self):
        test_input = ['-h']
        if logger.isEnabledFor(logging.DEBUG):
            main(test_input)
        with capture_stderr(main, test_input) as output:
            self.assertFalse(output)
        with capture_stdout(main, test_input) as output:
            self.assertTrue("optional arguments" in output)

    # These all test failure cases
    def testNoArgs(self):
        test_input = []
        # main(test_input)
        with capture_stderr(main, test_input) as output:
            self.assertTrue("Problems reading file" in output)

    def testNoFiles(self):
        test_input = ["-l", EMPTY_LIST]
        # main(test_input)
        with capture_stderr(main, test_input) as output:
            self.assertTrue("This program expects at least two files" in output)

    def testMissingFiles(self):
        test_input = ["-l", MISSING_FILE_LIST]
        # main(test_input)
        with capture_stderr(main, test_input) as output:
            self.assertTrue("Could not find the following file" in output)

    def testNoSuchOption(self):
        test_input = ["-@", LOG_LIST]
        # main(test_input)
        with capture_stderr(main, test_input) as output:
            self.assertTrue("unrecognized arguments" in output)


class TestGausslogUnique(unittest.TestCase):
    # These test/demonstrate different options
    def testStandard(self):
        test_input = ["-l", LOG_LIST]
        good_output = ''.join([HEADER, ALPHA_FIRST, ENERGY_FIRST]) + '\n'
        with capture_stdout(main, test_input) as output:
            self.assertTrue(output == good_output)
        with capture_stderr(main, test_input) as output:
            self.assertTrue('Check convergence' in output)

    # These test/demonstrate different options
    def testSortByEnergy(self):
        test_input = ["-l", LOG_LIST, "-e"]
        good_output = ''.join([HEADER, ENERGY_FIRST, ALPHA_FIRST]) + '\n'
        with capture_stdout(main, test_input) as output:
            self.assertTrue(output == good_output)
        with capture_stderr(main, test_input) as output:
            self.assertTrue('Check convergence' in output)

    def testSortByEnthalpy(self):
        test_input = ["-l", LOG_LIST, '-n']
        good_output = ''.join([HEADER, ALPHA_FIRST, ENERGY_FIRST]) + '\n'
        with capture_stdout(main, test_input) as output:
            self.assertTrue(output == good_output)
        with capture_stderr(main, test_input) as output:
            self.assertTrue('Check convergence' in output)

    def testNoFreq(self):
        # also tests that it can skip a blank line
        test_input = ["-l", LIST_NO_FREQ, "-n"]
        good_output = ''.join([HEADER, ENERGY_FIRST, MISSING_FREQ]) + '\n'
        main(test_input)
        with capture_stdout(main, test_input) as output:
            self.assertTrue(output == good_output)
        with capture_stderr(main, test_input) as output:
            self.assertFalse('Check convergence' in output)

    def testTwoMolecules(self):
        pet_843 = 'pet_mono_843_tzvp.log,1.4694,-917.071861,-916.796704\n'
        pet_1 = 'pet_mono_1_tzvp.log,0.8478,-917.069491,-916.794649\n'
        test_input = ["-l", TWO_MOL_LIST, "-n"]
        good_output = ''.join([HEADER, pet_843, pet_1, ALPHA_FIRST, ENERGY_FIRST]) + '\n'
        # main(test_input)
        with capture_stdout(main, test_input) as output:
            print(output)
            self.assertTrue(output == good_output)
        with capture_stderr(main, test_input) as output:
            self.assertTrue('lme2acetoxprpnt_ts3_ircf_opt.log' in output)

    def testTwoMoreMolecules(self):
        # This test checks that the program can handle the case when Gaussian prints '********' for convergence,
        #    when the convergence is so bad that it can't fit in the space allotted
        good_result = 'lme2acetoxprpnt_ts4_ircf_opt.log,2.0767,-535.576027,-535.401005\n' \
                      'lme2acetoxprpnt_ts4_b_ts_ircf_opt.log,8.2747,-535.575022,-535.399906\n'
        test_input = ["-l", TWO_MORE_MOL_LIST, "-n"]
        good_output = ''.join([HEADER, good_result]) + '\n'
        # main(test_input)
        with capture_stdout(main, test_input) as output:
            self.assertTrue(output == good_output)
        with capture_stderr(main, test_input) as output:
            self.assertTrue('lme2acetoxprpnt_ts4_ircf_opt.log' in output)
            self.assertTrue('lme2acetoxprpnt_ts4_b_ts_ircf_opt.log' in output)

    def testSimilarMolecules(self):
        # I was surprised that these weren't listed as the same; turns out, the difference in dihedral angle
        #  was almost, but not quite, 360; now, if within tolerance of 360 degrees, it will subtract 360, catching
        #  these similar conformations
        good_result = 'me2pheoxprpnt_30.log,0.0139,-613.945900,-613.726343\n\n'
        good_output = ''.join([HEADER, good_result])
        test_input = ["-l", SIMILAR_LIST, "-n"]
        with capture_stdout(main, test_input) as output:
            self.assertTrue(output == good_output)
        with capture_stderr(main, test_input) as output:
            self.assertTrue('' == output)

    def testCalcAllOutput(self):
        # make sure program sees the enthalpy when "CalcAll" is run instead of frequency
        good_result = 'hexyl_acrylate_239.log,1.1100,-503.005111,-502.751977\n' \
                      'hexyl_acrylate_419.log,0.0706,-503.004423,-502.751021\n\n'
        good_output = ''.join([HEADER, good_result])
        test_input = ["-l", CALCALL_LIST, "-n"]
        # main(test_input)
        with capture_stdout(main, test_input) as output:
            self.assertTrue(output == good_output)
        with capture_stderr(main, test_input) as output:
            self.assertTrue('' == output)


class TestGausslogUniqueFunctions(unittest.TestCase):
    def testOneFileMissingDihedralInfo(self):
        # testing for some specific problem encountered: when there is a freq only job (no opt) there is no
        #     dihedral info (which caused trouble when tried to compare it) and also prevented finding the
        #     convergence. Thus, checking that this specific error is resolved with the expected warning.
        current_file_list = ["tieg5ipatse_ts_ircf_opt_noeg.log", "tieg5ipatse_ts_noeg_reacta.log",
                             "tieg5ipatse_ts_noeg_reactc.log"]
        full_info_dict = {}
        check_list = []
        gausslog_loc = os.path.join(SUB_DATA_DIR, current_file_list[0])
        with capture_stderr(process_gausslog_file, gausslog_loc, find_dih=True, find_converg=True) as output:
            self.assertTrue("Requested dihedral data not found for file: tieg5ipatse_ts_ircf_opt_noeg.log" in
                            output)
        for gausslog_fname in current_file_list:
            if gausslog_fname == "":
                continue
            gausslog_loc = os.path.join(SUB_DATA_DIR, gausslog_fname)
            gausslog_content = process_gausslog_file(gausslog_loc, find_dih=True, find_converg=True)
            full_info_dict[gausslog_fname] = gausslog_content
            if gausslog_content[CONVERG_ERR]:
                check_list.append(gausslog_fname)
        self.assertTrue(full_info_dict[current_file_list[0]][DIHES] is None)
        self.assertTrue(full_info_dict[current_file_list[1]][DIHES] is not None)
        self.assertFalse(full_info_dict[current_file_list[0]][TS])
        self.assertTrue(full_info_dict[current_file_list[0]][CONVERG_ERR])
        self.assertTrue(full_info_dict[current_file_list[1]][CONVERG_ERR])
        self.assertEqual(len(check_list), 3)
        list_of_conf_lists = compare_gausslog_info(full_info_dict, 5)
        self.assertEqual(len(list_of_conf_lists), 3)
        winner_str, warn_files_str = print_results(full_info_dict, list_of_conf_lists, False, True, print_winners=False)
        warn_files_list = warn_files_str.split('\n')
        self.assertEqual(len(warn_files_list), 4)

    def testRemoveConvErr(self):
        # sometimes Gaussian ignore convergence error stating:
        #     "Optimization completed on the basis of negligible forces."
        # in those cases, let's also ignore the error and change the error flag to false
        gausslog_loc = os.path.join(SUB_DATA_DIR, "tieg5pdc2tsc_ts_ircr_opt.log")
        gausslog_content = process_gausslog_file(gausslog_loc, find_dih=True, find_converg=True)
        self.assertFalse(gausslog_content[CONVERG_ERR])
        self.assertAlmostEqual(gausslog_content[CONVERG], 1.5711111111111113)

    def testTSInfo(self):
        goodvibes_helper_folder = os.path.join(DATA_DIR, 'goodvibes_helper')
        file_list = ["co_gas.log", "hcoch3_gas.log", "tieg5ipatse_ts.log", "water.log", "h_gas_stable_t.log"]
        full_info_dict = {}
        check_list = []
        for fname in file_list:
            gausslog_loc = os.path.join(goodvibes_helper_folder, fname)
            gausslog_content = process_gausslog_file(gausslog_loc, find_dih=True, find_converg=True)
            full_info_dict[fname] = gausslog_content
            if gausslog_content[CONVERG_ERR]:
                check_list.append(fname)
        for fname in file_list:
            if fname == "tieg5ipatse_ts.log":
                self.assertTrue(full_info_dict[fname][TS])
            else:
                self.assertFalse(full_info_dict[fname][TS])
