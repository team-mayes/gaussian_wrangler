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
DEF_SH_OUT = os.path.join(SUB_DATA_DIR, 'ethylrad.sh')
GOOD_SH_OUT = os.path.join(SUB_DATA_DIR, 'good_ethylrad.sh')
DEF_LOG_OUT = os.path.join(PARENT_DIR, 'ethylrad.log')
ONE_JOB_INI = os.path.join(SUB_DATA_DIR, 'run_gauss_bde_one_job.ini')
GOOD_ONE_SH_OUT = os.path.join(SUB_DATA_DIR, 'good_ethylrad_one.sh')

MISSING_TPL_INI = os.path.join(SUB_DATA_DIR, 'run_gauss_missing_tpl.ini')
ONE_NEW_JOB_INI = os.path.join(SUB_DATA_DIR, 'run_gauss_one.ini')
OPT_LOG_OUT = os.path.join(PARENT_DIR, 'ethylrad_opt.log')
OPT_SH_OUT = os.path.join(SUB_DATA_DIR, 'ethylrad_opt.sh')
OPT_STABLE_LOG_OUT = os.path.join(PARENT_DIR, 'ethylrad_opt_stable.log')
OPT_STABLE_SH_OUT = os.path.join(SUB_DATA_DIR, 'ethylrad_opt_stable.sh')
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
        test_input = ["tests/test_data/run_gauss/ethylrad", "-c", MISSING_TPL_INI]
        if logger.isEnabledFor(logging.DEBUG):
            main(test_input)
        with capture_stderr(main, test_input) as output:
            self.assertTrue("could not find a template file" in output)

    def testMissingComFileIni(self):
        test_input = ["ghost", "-c", DEF_INI]
        if logger.isEnabledFor(logging.DEBUG):
            main(test_input)
        with capture_stderr(main, test_input) as output:
            self.assertTrue("Problems reading file" in output)


class TestRunGaussBDE(unittest.TestCase):
    # These test/demonstrate different options
    def testDefIni(self):
        test_input = ["tests/test_data/run_gauss/ethylrad", "-c", DEF_INI]
        try:
            main(test_input)
            self.assertFalse(diff_lines(DEF_SH_OUT, GOOD_SH_OUT))
        finally:
            silent_remove(DEF_SH_OUT, disable=DISABLE_REMOVE)
            silent_remove(DEF_LOG_OUT, disable=DISABLE_REMOVE)
            pass

    def testOneJobIni(self):
        test_input = ["tests/test_data/run_gauss/ethylrad", "-c", ONE_JOB_INI, "-t"]
        try:
            main(test_input)
            self.assertFalse(diff_lines(DEF_SH_OUT, GOOD_ONE_SH_OUT))
        finally:
            silent_remove(DEF_SH_OUT, disable=DISABLE_REMOVE)
            silent_remove(DEF_LOG_OUT, disable=DISABLE_REMOVE)
            pass

    def testOneNewJobIni(self):
        test_input = ["tests/test_data/run_gauss/ethylrad", "-c", ONE_NEW_JOB_INI, "-t"]
        try:
            main(test_input)
            self.assertFalse(diff_lines(OPT_SH_OUT, GOOD_OPT_SH_OUT))
        finally:
            silent_remove(OPT_SH_OUT, disable=DISABLE_REMOVE)
            silent_remove(OPT_LOG_OUT, disable=DISABLE_REMOVE)
            pass

    def testSpawn(self):
        test_input = ["tests/test_data/run_gauss/ethylrad", "-c", SPAWN_INI, "-t"]
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

    def testSpawnAllNew(self):
        test_input = ["tests/test_data/run_gauss/ethylrad", "-c", SPAWN_ALL_NEW_INI, "-t"]
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
