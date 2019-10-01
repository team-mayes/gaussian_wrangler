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

ETHYLRAD = 'tests/test_data/run_gauss/ethylrad'
DEF_INI = os.path.join(SUB_DATA_DIR, 'run_gauss_bde.ini')
DEF_SH_OUT = os.path.join(MAIN_DIR, 'ethylrad.sh')
GOOD_SH_OUT = os.path.join(SUB_DATA_DIR, 'good_ethylrad.sh')
DEF_LOG_OUT = os.path.join(PARENT_DIR, 'ethylrad.log')
ONE_JOB_INI = os.path.join(SUB_DATA_DIR, 'run_gauss_bde_one_job.ini')
GOOD_ONE_SH_OUT = os.path.join(SUB_DATA_DIR, 'good_ethylrad_one.sh')

MISSING_TPL_INI = os.path.join(SUB_DATA_DIR, 'run_gauss_missing_tpl.ini')
ONE_NEW_JOB_INI = os.path.join(SUB_DATA_DIR, 'run_gauss_one.ini')
OPT_LOG_OUT = os.path.join(PARENT_DIR, 'ethylrad_opt.log')
OPT_SH_OUT = os.path.join(MAIN_DIR, 'ethylrad_opt.sh')
OPT_STABLE_LOG_OUT = os.path.join(PARENT_DIR, 'ethylrad_opt_stable.log')
OPT_STABLE_SH_OUT = os.path.join(MAIN_DIR, 'ethylrad_opt_stable.sh')
GOOD_OPT_SH_OUT = os.path.join(SUB_DATA_DIR, 'good_ethylrad_opt.sh')
GOOD_OPT_STABLE_SH_OUT = os.path.join(SUB_DATA_DIR, 'good_ethylrad_opt_stable.sh')

SPAWN_INI = os.path.join(SUB_DATA_DIR, 'run_spawn.ini')
SPAWN1_INI = os.path.join(SUB_DATA_DIR, 'run_spawn1.ini')
SPAWN2_INI = os.path.join(SUB_DATA_DIR, 'run_spawn2.ini')
SPAWN1_SLURM = os.path.join(SUB_DATA_DIR, 'run_spawn1.slurm')
SPAWN2_SLURM = os.path.join(SUB_DATA_DIR, 'run_spawn2.slurm')
GOOD_SPAWN1_INI = os.path.join(SUB_DATA_DIR, 'run_spawn1_good.ini')
GOOD_SPAWN2_INI = os.path.join(SUB_DATA_DIR, 'run_spawn2_good.ini')
GOOD_SPAWN1_SLURM = os.path.join(SUB_DATA_DIR, 'run_spawn1_good.slurm')
GOOD_SPAWN2_SLURM = os.path.join(SUB_DATA_DIR, 'run_spawn2_good.slurm')

SUBMIT_NO_CHK_CHECK_INI = os.path.join(SUB_DATA_DIR, 'run_spawn_no_chk_chk.ini')
NO_CHK_CHECK_INI = os.path.join(MAIN_DIR, 'ethylrad_restart.ini')
NO_CHK_CHECK_SLM = os.path.join(MAIN_DIR, 'ethylrad_restart.slurm')

SPAWN_GIVE_OLD_CHK_STR_INI = os.path.join(SUB_DATA_DIR, 'run_spawn_give_old_chk_str.ini')
SPAWN1_GIVE_CHK_INI = os.path.join(SUB_DATA_DIR, 'run_spawn_give_old_chk_str1.ini')
SPAWN2_GIVE_CHK_INI = os.path.join(SUB_DATA_DIR, 'run_spawn_give_old_chk_str2.ini')
SPAWN1_GIVE_CHK_SLM = os.path.join(SUB_DATA_DIR, 'run_spawn_give_old_chk_str1.slurm')
SPAWN2_GIVE_CHK_SLM = os.path.join(SUB_DATA_DIR, 'run_spawn_give_old_chk_str2.slurm')

SPAWN_ALT_OLD_CHK_STR_INI = os.path.join(SUB_DATA_DIR, 'run_spawn_diff_old_chk_str.ini')
GOOD_ALT_OPT_STABLE_SH_OUT = os.path.join(SUB_DATA_DIR, 'good_ethylrad_opt_stable_alt.sh')
SPAWN1_DIFF_CHK_INI = os.path.join(SUB_DATA_DIR, 'run_spawn_diff_old_chk_str1.ini')
SPAWN2_DIFF_CHK_INI = os.path.join(SUB_DATA_DIR, 'run_spawn_diff_old_chk_str2.ini')
SPAWN1_DIFF_CHK_SLM = os.path.join(SUB_DATA_DIR, 'run_spawn_diff_old_chk_str1.slurm')
SPAWN2_DIFF_CHK_SLM = os.path.join(SUB_DATA_DIR, 'run_spawn_diff_old_chk_str2.slurm')

SPAWN_ALL_NEW_INI = os.path.join(SUB_DATA_DIR, 'run_spawn_all_new.ini')
SPAWN0_NEW_INI = os.path.join(SUB_DATA_DIR, 'run_spawn_all_new0.ini')
SPAWN1_NEW_INI = os.path.join(SUB_DATA_DIR, 'run_spawn_all_new1.ini')
SPAWN2_NEW_INI = os.path.join(SUB_DATA_DIR, 'run_spawn_all_new2.ini')
SPAWN0_NEW_SLURM = os.path.join(SUB_DATA_DIR, 'run_spawn_all_new0.slurm')
SPAWN1_NEW_SLURM = os.path.join(SUB_DATA_DIR, 'run_spawn_all_new1.slurm')
SPAWN2_NEW_SLURM = os.path.join(SUB_DATA_DIR, 'run_spawn_all_new2.slurm')
GOOD_SPAWN0_NEW_INI = os.path.join(SUB_DATA_DIR, 'run_spawn_all_new0_good.ini')
GOOD_SPAWN0_NEW_SLURM = os.path.join(SUB_DATA_DIR, 'run_spawn_all_new0_good.slurm')
GOOD_SPAWN1_NEW_SLURM = os.path.join(SUB_DATA_DIR, 'run_spawn_all_new1_good.slurm')
GOOD_SPAWN2_NEW_SLURM = os.path.join(SUB_DATA_DIR, 'run_spawn_all_new2_good.slurm')

SETUP_SUBMIT_INI = os.path.join(SUB_DATA_DIR, 'set_up_submit.ini')
SETUP_INI_OUT = os.path.join(MAIN_DIR, 'ethylrad.ini')
SETUP_SLURM_OUT = os.path.join(MAIN_DIR, 'ethylrad.slurm')
GOOD_SETUP_INI_OUT = os.path.join(SUB_DATA_DIR, 'good_ethylrad.ini')
GOOD_SETUP_SLURM_OUT = os.path.join(SUB_DATA_DIR, 'good_ethylrad.slurm')
WATER_INI_OUT = os.path.join(MAIN_DIR, 'water.ini')
WATER_SLURM_OUT = os.path.join(MAIN_DIR, 'water.slurm')
GOOD_WATER_INI_OUT = os.path.join(SUB_DATA_DIR, 'water_good.ini')
GOOD_WATER_SLURM_OUT = os.path.join(SUB_DATA_DIR, 'water_good.slurm')

LIST = os.path.join(SUB_DATA_DIR, "list.txt")
SETUP_SUBMIT_SPAWN_INI = os.path.join(SUB_DATA_DIR, 'set_up_submit_w_spawn.ini')
GOOD_ETHYL_SPAWN_INI_OUT = os.path.join(SUB_DATA_DIR, 'good_ethylrad_spawn.ini')
GOOD_ETHYL_SPAWN_SLURM_OUT = os.path.join(SUB_DATA_DIR, 'good_ethylrad_spawn.slurm')
GOOD_WATER_SPAWN_SLURM_OUT = os.path.join(SUB_DATA_DIR, 'water_spawn_good.slurm')

SETUP_DEF_TPL_INI = os.path.join(SUB_DATA_DIR, 'submit_current_f_ts.ini')
SETUP_INI_DEF_DIR_OUT = os.path.join(MAIN_DIR, 'ethylrad.ini')
SETUP_SLURM_DEF_DIR_OUT = os.path.join(MAIN_DIR, 'ethylrad.slurm')
GOOD_DEF_DIR_INI_OUT = os.path.join(SUB_DATA_DIR, 'ethylrad_default_dir_good.ini')
GOOD_DEF_DIR_SLURM_OUT = os.path.join(SUB_DATA_DIR, 'ethylrad_default_dir_good.slurm')
GOOD_SHORT_DEF_DIR_SLURM_OUT = os.path.join(SUB_DATA_DIR, 'ethylrad_short_def_dir_good.slurm')

SETUP_IRCS_INI = os.path.join(SUB_DATA_DIR, 'submit_ircs_opt.ini')
SETUP_IRCR_INI_OUT = os.path.join(MAIN_DIR, 'ethylrad0.ini')
SETUP_IRCR_SLM_OUT = os.path.join(MAIN_DIR, 'ethylrad0.slurm')
SETUP_IRCF_INI_OUT = os.path.join(MAIN_DIR, 'ethylrad1.ini')
SETUP_IRCF_SLM_OUT = os.path.join(MAIN_DIR, 'ethylrad1.slurm')
GOOD_SETUP_IRCR_INI_OUT = os.path.join(SUB_DATA_DIR, 'ethylrad0_good.ini')
GOOD_SETUP_IRCR_SLM_OUT = os.path.join(SUB_DATA_DIR, 'ethylrad0_good.slurm')
GOOD_SETUP_IRCF_INI_OUT = os.path.join(SUB_DATA_DIR, 'ethylrad1_good.ini')
GOOD_SETUP_IRCF_SLM_OUT = os.path.join(SUB_DATA_DIR, 'ethylrad1_good.slurm')


class TestRunGaussNoOut(unittest.TestCase):
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
        test_input = [ETHYLRAD, "-c", MISSING_TPL_INI]
        # main(test_input)
        if logger.isEnabledFor(logging.DEBUG):
            main(test_input)
        with capture_stderr(main, test_input) as output:
            self.assertTrue("not find the submit template" in output)

    def testMissingComFileIni(self):
        test_input = ["ghost", "-c", DEF_INI]
        if logger.isEnabledFor(logging.DEBUG):
            main(test_input)
        with capture_stderr(main, test_input) as output:
            self.assertTrue("Problems reading file" in output)

    def testAttemptRunMultIni(self):
        # Create and submit more than one ini
        temp_file_list = ['ethylrad.com', 'ircr.tpl', 'ircf.tpl', 'opt.tpl']
        for fname in temp_file_list:
            with open(fname, 'w') as f:
                f.write("# for test only")
        test_input = ['ethylrad', "-c", SETUP_IRCS_INI, "-o", 'ethyl']
        if logger.isEnabledFor(logging.DEBUG):
            main(test_input)
        try:
            # main(test_input)
            with capture_stderr(main, test_input) as output:
                self.assertTrue("only supported" in output)
        finally:
            for fname in temp_file_list + [SETUP_INI_DEF_DIR_OUT, SETUP_SLURM_DEF_DIR_OUT]:
                silent_remove(fname, disable=DISABLE_REMOVE)
            pass

    def testNoChk(self):
        # Checking alternate input from above
        temp_file_list = ['ethylrad.com', 'f.tpl', 'ts.tpl']
        for fname in temp_file_list:
            with open(fname, 'w') as f:
                f.write("# for test only")

        test_input = ['tests/test_data/run_gauss/ethylrad_restart', "-c", SPAWN_INI, "-s"]
        try:
            # main(test_input)
            with capture_stderr(main, test_input) as output:
                self.assertTrue("old checkpoint" in output)
        finally:
            for fname in temp_file_list:
                silent_remove(fname, disable=DISABLE_REMOVE)
            pass


class TestRunGauss(unittest.TestCase):
    # These test/demonstrate different options
    def testDefIni(self):
        test_input = [ETHYLRAD, "-c", DEF_INI]
        try:
            main(test_input)
            self.assertFalse(diff_lines(DEF_SH_OUT, GOOD_SH_OUT))
        finally:
            silent_remove(DEF_SH_OUT, disable=DISABLE_REMOVE)
            silent_remove(DEF_LOG_OUT, disable=DISABLE_REMOVE)
            pass

    def testOneJobIni(self):
        test_input = [ETHYLRAD, "-c", ONE_JOB_INI, "-t"]
        try:
            main(test_input)
            self.assertFalse(diff_lines(DEF_SH_OUT, GOOD_ONE_SH_OUT))
        finally:
            silent_remove(DEF_SH_OUT, disable=DISABLE_REMOVE)
            silent_remove(DEF_LOG_OUT, disable=DISABLE_REMOVE)
            pass

    def testOneNewJobIni(self):
        test_input = [ETHYLRAD, "-c", ONE_NEW_JOB_INI, "-t"]
        try:
            main(test_input)
            self.assertFalse(diff_lines(OPT_SH_OUT, GOOD_OPT_SH_OUT))
        finally:
            silent_remove(OPT_SH_OUT, disable=DISABLE_REMOVE)
            silent_remove(OPT_LOG_OUT, disable=DISABLE_REMOVE)
            pass

    def testSpawn(self):
        test_input = [ETHYLRAD, "-c", SPAWN_INI, "-t"]
        try:
            main(test_input)
            self.assertFalse(diff_lines(OPT_STABLE_SH_OUT, GOOD_OPT_STABLE_SH_OUT))
            self.assertFalse(diff_lines(SPAWN1_INI, GOOD_SPAWN1_INI))
            self.assertFalse(diff_lines(SPAWN2_INI, GOOD_SPAWN2_INI))
            self.assertFalse(diff_lines(SPAWN1_SLURM, GOOD_SPAWN1_SLURM))
            self.assertFalse(diff_lines(SPAWN2_SLURM, GOOD_SPAWN2_SLURM))
        finally:
            for fname in [DEF_SH_OUT, OPT_SH_OUT, OPT_STABLE_SH_OUT, DEF_LOG_OUT, OPT_LOG_OUT, OPT_STABLE_LOG_OUT,
                          SPAWN1_INI, SPAWN2_INI, SPAWN1_SLURM, SPAWN2_SLURM]:
                silent_remove(fname, disable=DISABLE_REMOVE)
            pass

    def testSpawnGiveChkStr(self):
        test_input = [ETHYLRAD, "-c", SPAWN_GIVE_OLD_CHK_STR_INI, "-t"]
        try:
            main(test_input)
            self.assertFalse(diff_lines(OPT_STABLE_SH_OUT, GOOD_OPT_STABLE_SH_OUT))
        finally:
            for fname in [DEF_SH_OUT, OPT_SH_OUT, OPT_STABLE_SH_OUT, DEF_LOG_OUT, OPT_LOG_OUT, OPT_STABLE_LOG_OUT,
                          SPAWN1_GIVE_CHK_INI, SPAWN2_GIVE_CHK_INI, SPAWN1_GIVE_CHK_SLM, SPAWN2_GIVE_CHK_SLM]:
                silent_remove(fname, disable=DISABLE_REMOVE)
            pass

    def testSpawnAltChkStr(self):
        test_input = [ETHYLRAD, "-c", SPAWN_ALT_OLD_CHK_STR_INI, "-t"]
        try:
            main(test_input)
            self.assertFalse(diff_lines(OPT_STABLE_SH_OUT, GOOD_ALT_OPT_STABLE_SH_OUT))
        finally:
            for fname in [DEF_SH_OUT, OPT_SH_OUT, OPT_STABLE_SH_OUT, DEF_LOG_OUT, OPT_LOG_OUT, OPT_STABLE_LOG_OUT,
                          SPAWN1_DIFF_CHK_INI, SPAWN2_DIFF_CHK_INI, SPAWN1_DIFF_CHK_SLM, SPAWN2_DIFF_CHK_SLM]:
                silent_remove(fname, disable=DISABLE_REMOVE)
            pass

    def testNoChkChk(self):
        test_input = ['tests/test_data/run_gauss/ethylrad_restart', "-s", "-c", SUBMIT_NO_CHK_CHECK_INI, "-t"]
        try:
            main(test_input)
            # TODO: Check that you can submit a spawn job!
            # self.assertFalse(diff_lines(OPT_STABLE_SH_OUT, GOOD_ALT_OPT_STABLE_SH_OUT))
        finally:
            for fname in [NO_CHK_CHECK_INI, NO_CHK_CHECK_SLM]:
                silent_remove(fname, disable=DISABLE_REMOVE)
            pass

    def testSpawnAllNew(self):
        test_input = [ETHYLRAD, "-c", SPAWN_ALL_NEW_INI, "-t"]
        try:
            main(test_input)
            self.assertFalse(diff_lines(SPAWN0_NEW_INI, GOOD_SPAWN0_NEW_INI))
            self.assertFalse(diff_lines(SPAWN1_NEW_INI, GOOD_SPAWN1_INI))
            self.assertFalse(diff_lines(SPAWN2_NEW_INI, GOOD_SPAWN2_INI))
            self.assertFalse(diff_lines(SPAWN0_NEW_SLURM, GOOD_SPAWN0_NEW_SLURM))
            self.assertFalse(diff_lines(SPAWN1_NEW_SLURM, GOOD_SPAWN1_NEW_SLURM))
            self.assertFalse(diff_lines(SPAWN2_NEW_SLURM, GOOD_SPAWN2_NEW_SLURM))
        finally:
            for fname in [DEF_SH_OUT, OPT_SH_OUT, DEF_LOG_OUT, OPT_LOG_OUT,
                          SPAWN0_NEW_INI, SPAWN1_NEW_INI, SPAWN2_NEW_INI,
                          SPAWN0_NEW_SLURM, SPAWN1_NEW_SLURM, SPAWN2_NEW_SLURM]:
                silent_remove(fname, disable=DISABLE_REMOVE)
            pass

    def testSetupSubmit(self):
        test_input = [ETHYLRAD, "-s", "-c", SETUP_SUBMIT_INI]
        try:
            main(test_input)
            self.assertFalse(diff_lines(SETUP_INI_OUT, GOOD_SETUP_INI_OUT))
            self.assertFalse(diff_lines(SETUP_SLURM_OUT, GOOD_SETUP_SLURM_OUT))
        finally:
            for fname in [SETUP_INI_OUT, SETUP_SLURM_OUT]:
                silent_remove(fname, disable=DISABLE_REMOVE)
            pass

    def testSetupSubmitList(self):
        temp_file_list = ['ethylrad.com', 'water.com']
        for fname in temp_file_list:
            with open(fname, 'w') as f:
                f.write("# for test only\n\n")
        test_input = [LIST, "-l", "-c", SETUP_SUBMIT_INI]
        try:
            main(test_input)
            self.assertFalse(diff_lines(SETUP_INI_OUT, GOOD_SETUP_INI_OUT))
            self.assertFalse(diff_lines(SETUP_SLURM_OUT, GOOD_SETUP_SLURM_OUT))
            self.assertFalse(diff_lines(WATER_INI_OUT, GOOD_WATER_INI_OUT))
            self.assertFalse(diff_lines(WATER_SLURM_OUT, GOOD_WATER_SLURM_OUT))
        finally:
            for fname in temp_file_list + [SETUP_INI_OUT, SETUP_SLURM_OUT, WATER_INI_OUT, WATER_SLURM_OUT]:
                silent_remove(fname, disable=DISABLE_REMOVE)
            pass

    def testSetupSubmitListSpawn(self):
        test_input = [LIST, "-l", "-c", SETUP_SUBMIT_SPAWN_INI]
        try:
            main(test_input)
            self.assertFalse(diff_lines(SETUP_INI_OUT, GOOD_ETHYL_SPAWN_INI_OUT))
            self.assertFalse(diff_lines(SETUP_SLURM_OUT, GOOD_ETHYL_SPAWN_SLURM_OUT))
            self.assertFalse(diff_lines(WATER_INI_OUT, GOOD_ETHYL_SPAWN_INI_OUT))
            self.assertFalse(diff_lines(WATER_SLURM_OUT, GOOD_WATER_SPAWN_SLURM_OUT))
        finally:
            for fname in [SETUP_INI_OUT, SETUP_SLURM_OUT, WATER_INI_OUT, WATER_SLURM_OUT]:
                silent_remove(fname, disable=DISABLE_REMOVE)
            pass

    def testSubmitWithChkDefDirIni(self):
        # since this assumes a f.tpl and ts.tpl file in the main directory, create them, and delete at the end
        temp_file_list = ['f.tpl', 'ts.tpl']
        for fname in temp_file_list:
            with open(fname, 'w') as f:
                f.write("# for test only")

        test_input = ['tests/test_data/run_gauss/ethylrad.com', "-c", SETUP_DEF_TPL_INI, "-s",
                      "-o", 'tests/test_data/run_gauss/ethyl.chk']
        try:
            main(test_input)
            self.assertFalse(diff_lines(SETUP_INI_DEF_DIR_OUT, GOOD_DEF_DIR_INI_OUT))
            self.assertFalse(diff_lines(SETUP_SLURM_DEF_DIR_OUT, GOOD_DEF_DIR_SLURM_OUT))
        finally:
            for fname in temp_file_list + [SETUP_INI_DEF_DIR_OUT, SETUP_SLURM_DEF_DIR_OUT]:
                silent_remove(fname, disable=DISABLE_REMOVE)
            pass

    def testSubmitDefDirIni(self):
        # Checking alternate input from above
        temp_file_list = ['ethylrad.com', 'f.tpl', 'ts.tpl']
        for fname in temp_file_list:
            with open(fname, 'w') as f:
                f.write("# for test only")

        test_input = ['ethylrad', "-c", SETUP_DEF_TPL_INI, "-s", "-o", 'ethyl']
        try:
            with capture_stdout(main, test_input) as output:
                self.assertTrue("sbatch" in output)
            self.assertFalse(diff_lines(SETUP_INI_DEF_DIR_OUT, GOOD_DEF_DIR_INI_OUT))
            self.assertFalse(diff_lines(SETUP_SLURM_DEF_DIR_OUT, GOOD_SHORT_DEF_DIR_SLURM_OUT))
        finally:
            for fname in temp_file_list + [SETUP_INI_DEF_DIR_OUT, SETUP_SLURM_DEF_DIR_OUT]:
                silent_remove(fname, disable=DISABLE_REMOVE)
            pass

    def testSubmitMultIni(self):
        # Create and submit more than one ini
        temp_file_list = ['ethylrad.com', 'ircr.tpl', 'ircf.tpl', 'opt.tpl']
        for fname in temp_file_list:
            with open(fname, 'w') as f:
                f.write("# for test only")
        test_input = ['ethylrad', "-c", SETUP_IRCS_INI, "-s", ]
        try:
            main(test_input)
            self.assertFalse(diff_lines(SETUP_IRCR_INI_OUT, GOOD_SETUP_IRCR_INI_OUT))
            self.assertFalse(diff_lines(SETUP_IRCR_SLM_OUT, GOOD_SETUP_IRCR_SLM_OUT))
            self.assertFalse(diff_lines(SETUP_IRCF_INI_OUT, GOOD_SETUP_IRCF_INI_OUT))
            self.assertFalse(diff_lines(SETUP_IRCF_SLM_OUT, GOOD_SETUP_IRCF_SLM_OUT))
        finally:
            for fname in temp_file_list + [SETUP_IRCR_INI_OUT, SETUP_IRCR_SLM_OUT,
                                           SETUP_IRCF_INI_OUT, SETUP_IRCF_SLM_OUT]:
                silent_remove(fname, disable=DISABLE_REMOVE)
            pass
