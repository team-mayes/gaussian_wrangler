import unittest
import os
from gaussian_wrangler.gausscom2com import main
from common_wrangler.common import diff_lines, silent_remove, capture_stdout, capture_stderr
import logging

# logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger(__name__)
DISABLE_REMOVE = logger.isEnabledFor(logging.DEBUG)

__author__ = 'hmayes'

TEST_DIR = os.path.dirname(__file__)
MAIN_DIR = os.path.dirname(TEST_DIR)
DATA_DIR = os.path.join(os.path.dirname(__file__), 'test_data')
SUB_DATA_DIR = os.path.join(DATA_DIR, 'gausscom2com')

TPL_FILE = os.path.join(SUB_DATA_DIR, 'acyl-min_ts.tpl')
LIST_FILE = os.path.join(SUB_DATA_DIR, 'list.txt')
COM_IN_FILE = os.path.join(SUB_DATA_DIR, 'acyl-min_ts_110.gjf')
COM_OUT_FILE = os.path.join(SUB_DATA_DIR, 'acyl-min_ts_110.com')
GOOD_COM_OUT_FILE = os.path.join(SUB_DATA_DIR, 'acyl-min_ts_110_good.com')
COM_OUT_FILE2 = os.path.join(SUB_DATA_DIR, 'acyl-min_ts_138.com')
GOOD_COM_OUT_FILE2 = os.path.join(SUB_DATA_DIR, 'acyl-min_ts_138_good.com')

MISSING_MULT_COM_OUT_FILE = os.path.join(SUB_DATA_DIR, 'acyl-min_ts_110_missing_mult.com')
PINNED_TPL_FILE = os.path.join(SUB_DATA_DIR, 'acyl-min_ts_pinned.tpl')
GOOD_PINNED_OUT_FILE = os.path.join(SUB_DATA_DIR, 'acyl-min_ts_pinned_good.com')

LOG_DATA_DIR = os.path.join(DATA_DIR, 'gausslog2com')
CP_COM_TPL = os.path.join(LOG_DATA_DIR, 'cp.tpl')
CP_COM_IN = os.path.join(SUB_DATA_DIR, 'pet_mono_901.gjf')
CP_COM_OUT = os.path.join(SUB_DATA_DIR, 'pet_mono_901.com')
GOOD_CP_COM_OUT = os.path.join(SUB_DATA_DIR, 'pet_mono_901_good.com')

CP_COM_IN_MISSING_ATOM = os.path.join(SUB_DATA_DIR, 'pet_mono_901_missing_atom.gjf')
COM_IN_NO_CHARGE_FILE = os.path.join(SUB_DATA_DIR, 'acyl-min_ts_110_missing_mult.gjf')

ROUTE_ONLY_TPL = os.path.join(SUB_DATA_DIR, 'route_only.tpl')
ROUTE_NO_CHARGE_TPL = os.path.join(SUB_DATA_DIR, 'route_no_charge.tpl')
GOOD_ROUTE_ONLY_COM_OUT = os.path.join(SUB_DATA_DIR, 'pet_mono_901_min_tpl_good.com')


PINNED_ATOM_IN = os.path.join(SUB_DATA_DIR, 'acyl-min_ts_pinned.gjf')
PINNED_ATOM_OUT = os.path.join(SUB_DATA_DIR, 'acyl-min_ts_pinned.com')


class TestGausscom2comNoOut(unittest.TestCase):
    # These all test failure cases
    def testNoArgs(self):
        test_input = []
        # main(test_input)
        with capture_stderr(main, test_input) as output:
            self.assertTrue("A template file" in output)

    def testHelp(self):
        test_input = ['-h']
        if logger.isEnabledFor(logging.DEBUG):
            main(test_input)
        with capture_stderr(main, test_input) as output:
            self.assertFalse(output)
        with capture_stdout(main, test_input) as output:
            self.assertTrue("optional arguments" in output)

    def testUnknownArg(self):
        test_input = ['--ghost']
        if logger.isEnabledFor(logging.DEBUG):
            main(test_input)
        with capture_stderr(main, test_input) as output:
            self.assertTrue("unrecognized arguments" in output)
        with capture_stdout(main, test_input) as output:
            self.assertTrue("optional arguments" in output)

    def testNoFilesToProcess(self):
        test_input = ["-t", TPL_FILE]
        main(test_input)
        if logger.isEnabledFor(logging.DEBUG):
            main(test_input)
        with capture_stderr(main, test_input) as output:
            self.assertTrue("No files have been specified to be read" in output)

    def testMissingFile(self):
        test_input = ["-t", TPL_FILE, "-f", "ghost.txt"]
        # main(test_input)
        if logger.isEnabledFor(logging.DEBUG):
            main(test_input)
        with capture_stderr(main, test_input) as output:
            self.assertTrue("Problems reading file" in output)

    def testMisMatchTplCom(self):
        test_input = ["-t", TPL_FILE, "-f", CP_COM_IN]
        # main(test_input)
        with capture_stderr(main, test_input) as output:
            self.assertTrue("Atom types do not match" in output)

    def testMismatchAtomNum(self):
        test_input = ["-t", CP_COM_TPL, "-f", CP_COM_IN_MISSING_ATOM]
        # main(test_input)
        with capture_stderr(main, test_input) as output:
            self.assertTrue("while the tpl file" in output)

    def testOddNumChargeMultLine(self):
        test_input = ["-t", TPL_FILE, "-f", COM_IN_NO_CHARGE_FILE, "-c"]
        # main(test_input)
        with capture_stderr(main, test_input) as output:
            self.assertTrue("read charge and multiplicity" in output)


class TestGausscom2com(unittest.TestCase):
    # These test/demonstrate different options
    def testFromListCase(self):
        test_input = ["-t", TPL_FILE, "-l", LIST_FILE, "-c"]
        try:
            main(test_input)
            self.assertFalse(diff_lines(COM_OUT_FILE, GOOD_COM_OUT_FILE))
            self.assertFalse(diff_lines(COM_OUT_FILE2, GOOD_COM_OUT_FILE2))
        finally:
            silent_remove(COM_OUT_FILE, disable=DISABLE_REMOVE)
            silent_remove(COM_OUT_FILE2, disable=DISABLE_REMOVE)
            pass

    def testFromFileCase(self):
        test_input = ["-t", TPL_FILE, "-f", COM_IN_FILE]
        try:
            main(test_input)
            self.assertFalse(diff_lines(COM_OUT_FILE, GOOD_COM_OUT_FILE))
        finally:
            silent_remove(COM_OUT_FILE, disable=DISABLE_REMOVE)
            pass

    def testPinnedTpl(self):
        test_input = ["-t", PINNED_TPL_FILE, "-f", COM_IN_FILE]
        try:
            main(test_input)
            self.assertFalse(diff_lines(COM_OUT_FILE, GOOD_PINNED_OUT_FILE))
        finally:
            silent_remove(COM_OUT_FILE, disable=DISABLE_REMOVE)
            pass

    def testFromPinnedAtomFileCase(self):
        test_input = ["-t", TPL_FILE, "-f", PINNED_ATOM_IN]
        try:
            main(test_input)
            self.assertFalse(diff_lines(PINNED_ATOM_OUT, GOOD_COM_OUT_FILE))
        finally:
            silent_remove(PINNED_ATOM_OUT, disable=DISABLE_REMOVE)
            pass

    def testFileCPTpl(self):
        test_input = ["-t", CP_COM_TPL, "-f", CP_COM_IN]
        try:
            main(test_input)
            self.assertFalse(diff_lines(CP_COM_OUT, GOOD_CP_COM_OUT))
        finally:
            silent_remove(CP_COM_OUT, disable=DISABLE_REMOVE)
            pass

    def testFromFileCase(self):
        test_input = ["-t", TPL_FILE, "-f", COM_IN_NO_CHARGE_FILE]
        try:
            main(test_input)
            self.assertFalse(diff_lines(MISSING_MULT_COM_OUT_FILE, GOOD_COM_OUT_FILE))
        finally:
            silent_remove(MISSING_MULT_COM_OUT_FILE, disable=DISABLE_REMOVE)
            pass

    def testRouteOnlyTpl(self):
        test_input = ["-t", ROUTE_ONLY_TPL, "-f", CP_COM_IN]
        try:
            main(test_input)
            self.assertFalse(diff_lines(CP_COM_OUT, GOOD_ROUTE_ONLY_COM_OUT))
        finally:
            silent_remove(CP_COM_OUT, disable=DISABLE_REMOVE)
            pass

    def testRouteNoChargeTpl(self):
        test_input = ["-t", ROUTE_NO_CHARGE_TPL, "-f", CP_COM_IN, "-c"]
        try:
            main(test_input)
            self.assertFalse(diff_lines(CP_COM_OUT, GOOD_ROUTE_ONLY_COM_OUT))
        finally:
            silent_remove(CP_COM_OUT, disable=DISABLE_REMOVE)
            pass
