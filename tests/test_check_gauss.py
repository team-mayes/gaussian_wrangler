import unittest
import os
import numpy as np
import matplotlib.pyplot as plt
from shutil import copyfile
from scipy.optimize import curve_fit
from gaussian_wrangler.check_gauss import main, plot_scan, process_scan_array, charmm_dihedral, find_good_fit, \
    find_stable_points, collect_output_scan_steps
from common_wrangler.common import (capture_stdout, capture_stderr, diff_lines, silent_remove, list_to_file, make_dir)
import logging

# logging.basicConfig(level=logging.DEBUG)
from gw_common import process_gausslog_file, SCAN_DICT

logger = logging.getLogger(__name__)
DISABLE_REMOVE = logger.isEnabledFor(logging.DEBUG)

__author__ = 'hmayes'

TEST_DIR = os.path.dirname(__file__)
MAIN_DIR = os.path.dirname(TEST_DIR)
DATA_DIR = os.path.join(TEST_DIR, 'test_data')
SUB_DATA_DIR = os.path.join(DATA_DIR, 'check_gauss')
ALT_DATA_DIR = os.path.join(DATA_DIR, 'gausslog_unique')
LOG_2_COM_DIR = os.path.join(DATA_DIR, 'gausslog2com')
SUB_SUB_DIR = os.path.join(SUB_DATA_DIR, 'temp_dir')

FOR_HARTREE_DIR = os.path.join(MAIN_DIR, 'for_hartree')
SINGLE_FILE = os.path.join(SUB_DATA_DIR, 'me2propprpnt_7.log')

LIST_FILE = os.path.join(SUB_DATA_DIR, 'list.txt')
CONV_239_PNG = os.path.join(ALT_DATA_DIR, 'hexyl_acrylate_239_conv_steps.png')
CONV_419_PNG = os.path.join(ALT_DATA_DIR, 'hexyl_acrylate_419_conv_steps.png')
CONV_239_OUT = os.path.join(ALT_DATA_DIR, 'hexyl_acrylate_239_conv_steps.csv')
CONV_419_OUT = os.path.join(ALT_DATA_DIR, 'hexyl_acrylate_419_conv_steps.csv')
GOOD_CONV_239_OUT = os.path.join(ALT_DATA_DIR, 'hexyl_acrylate_239_conv_steps_good.csv')
GOOD_CONV_419_OUT = os.path.join(ALT_DATA_DIR, 'hexyl_acrylate_419_conv_steps_good.csv')

DIOXOLAN_OUT = os.path.join(SUB_DATA_DIR, 'dioxolan4ol_ts4_ts_conv_steps.csv')
GOOD_DIOXOLAN_OUT = os.path.join(SUB_DATA_DIR, 'dioxolan4ol_ts4_ts_conv_steps_good.csv')
DIOXOLAN_PNG = os.path.join(SUB_DATA_DIR, 'dioxolan4ol_ts4_ts_conv_steps.png')

IN_PROGRESS_OUT = os.path.join(SUB_DATA_DIR, 'ti_eg5_dime_tpa_tse_conv_steps.csv')
IN_PROGRESS_PNG = os.path.join(SUB_DATA_DIR, 'ti_eg5_dime_tpa_tse_conv_steps.png')

NOT_FROM_CHK_FILE = os.path.join(SUB_DATA_DIR, 'acyl-min_ts5.out')


# data

SCAN_ARRAY = np.asarray([[0., 7.233], [10., 6.60829], [20., 5.29499], [30., 3.6031], [40., 1.9289], [50., 0.64086],
                         [60., 6.275e-06], [70., 0.0298], [80., 0.5979], [90., 1.497], [100., 2.4848], [110., 3.2781],
                         [120., 3.6707], [130., 3.577], [140., 3.07964], [150., 2.357], [160., 1.6372], [170., 1.1119],
                         [180., 0.926], [190., 1.1118], [200., 1.6371], [210., 2.357], [220., 3.0797], [230., 3.577],
                         [240., 3.6707], [250., 3.2781], [260., 2.4849], [270., 1.4971], [280., 0.5979], [290., 0.0298],
                         [300., 0.], [310., 0.641], [320., 1.93], [330., 3.6], [340., 5.3], [350., 6.61], [360., 7.23]])


def make_fill_sub_dir():
    # Start fresh with sub_sub_dir
    silent_remove(SUB_SUB_DIR, dir_with_files=True)
    os.makedirs(SUB_SUB_DIR)
    base_com = 'a579.com'
    orig_com = os.path.join(LOG_2_COM_DIR, base_com)
    temp_com = os.path.join(SUB_SUB_DIR, base_com)
    copyfile(orig_com, temp_com)
    base_log = 'a579.log'
    orig_log = os.path.join(LOG_2_COM_DIR, base_log)
    temp_log = os.path.join(SUB_SUB_DIR, base_log)
    copyfile(orig_log, temp_log)
    base_cp_log = 'ipvc_11_10_cp.log'
    orig_cp_log = os.path.join(LOG_2_COM_DIR, base_cp_log)
    temp_cp_log = os.path.join(SUB_SUB_DIR, base_cp_log)
    copyfile(orig_cp_log, temp_cp_log)
    temp_files_list = [temp_com, temp_log, temp_cp_log]
    return temp_files_list


def test_curve_func(x, a, b):
    return a*x*x + b


class TestFittingCurve(unittest.TestCase):
    def testFit(self):
        x_data = np.linspace(0, 4, 50)
        y = test_curve_func(x_data, 2.5, 1.3)
        np.random.seed(1729)
        y_noise = 0.2 * np.random.normal(size=x_data.size)
        y_data = y + y_noise
        b = 2.2
        png_out = os.path.join(SUB_DATA_DIR, "test.png")
        try:
            silent_remove(png_out)
            popt, pcov = curve_fit(lambda x, a: test_curve_func(x, a, b), x_data, y_data)
            a_fit = popt[0]
            plt.plot(x_data, y_data, 'b-', label='data')
            plt.plot(x_data, test_curve_func(x_data, a_fit, b), 'r-', label=f'fit: a={a_fit:5.3f}')
            plt.legend()
            plt.savefig(png_out, transparent=True, bbox_inches='tight',)
            plt.close()
            self.assertTrue(os.path.isfile(png_out))
        finally:
            silent_remove(png_out, disable=DISABLE_REMOVE)


class TestCheckGaussNoOut(unittest.TestCase):
    # These test failure cases that do not produce output files
    def testNoArgs(self):
        test_input = []
        # main(test_input)
        with capture_stderr(main, test_input) as output:
            self.assertTrue("Could not find files" in output)
        silent_remove(FOR_HARTREE_DIR, disable=DISABLE_REMOVE)

    def testHelp(self):
        test_input = ['-h']
        if logger.isEnabledFor(logging.DEBUG):
            main(test_input)
        with capture_stderr(main, test_input) as output:
            self.assertFalse(output)
        with capture_stdout(main, test_input) as output:
            self.assertTrue("optional arguments" in output)
        silent_remove(FOR_HARTREE_DIR, disable=DISABLE_REMOVE)

    def testWrongKey(self):
        test_input = ['-ghost']
        if logger.isEnabledFor(logging.DEBUG):
            main(test_input)
        with capture_stderr(main, test_input) as output:
            self.assertTrue("unrecognized arguments" in output)
        with capture_stdout(main, test_input) as output:
            self.assertTrue("optional arguments" in output)
        silent_remove(FOR_HARTREE_DIR, disable=DISABLE_REMOVE)

    def testWrongDir(self):
        test_input = ["-d", "ghost"]
        if logger.isEnabledFor(logging.DEBUG):
            main(test_input)
        with capture_stderr(main, test_input) as output:
            self.assertTrue("Could not find" in output)
        silent_remove(FOR_HARTREE_DIR, disable=DISABLE_REMOVE)

    def testConflictingOptions(self):
        test_input = ["-s", "-z"]
        # main(test_input)
        # if logger.isEnabledFor(logging.DEBUG):
        #     main(test_input)
        with capture_stderr(main, test_input) as output:
            self.assertTrue("Choose either" in output)
        silent_remove(FOR_HARTREE_DIR, disable=DISABLE_REMOVE)

    def testNoLastStepNum(self):
        test_input = ["-t"]
        # main(test_input)
        # if logger.isEnabledFor(logging.DEBUG):
        #     main(test_input)
        with capture_stderr(main, test_input) as output:
            self.assertTrue("expected one argument" in output)

    def testNonIntLastStepNum(self):
        test_input = ["-t", "ghost"]
        # main(test_input)
        # if logger.isEnabledFor(logging.DEBUG):
        #     main(test_input)
        with capture_stderr(main, test_input) as output:
            self.assertTrue("integer must be provided" in output)

    def testDirInsteadOfFile(self):
        test_input = ["-f", SUB_DATA_DIR, "-z"]
        with capture_stderr(main, test_input) as output:
            self.assertTrue("Problems reading file" in output)

    def testNoScanInfo(self):
        no_scan_log = os.path.join(SUB_DATA_DIR, 'me2propprpnt_7.log')
        test_input = ["-f", no_scan_log, "--scan", "test.png"]
        # main(test_input)
        with capture_stderr(main, test_input) as output:
            self.assertTrue("Did not find expected parameter scan info" in output)

    def testMismatchScanFiles(self):
        # no_scan_log = os.path.join(SUB_DATA_DIR, 'me2propprpnt_7.log')
        # test_input = ["-f", no_scan_log, "--scan", "test.png"]
        # # main(test_input)
        # with capture_stderr(main, test_input) as output:
        #     self.assertTrue("Did not find expected parameter scan info" in output)

        list_fname = os.path.join(SUB_DATA_DIR, "scan_list.txt")
        fnames = [os.path.join(SUB_DATA_DIR, "tieg4pdc1scan.log"),
                  os.path.join(SUB_DATA_DIR, "pet_dimer_scan_pos_tzvp.log")]
        list_to_file(fnames, list_fname)
        out_png_fname = os.path.join(SUB_DATA_DIR, "pet_dimer_scan.png")
        test_input = ["-l", list_fname, "--scan", out_png_fname]
        try:
            # main(test_input)
            with capture_stderr(main, test_input) as output:
                self.assertTrue("cannot" in output)
        finally:
            silent_remove(list_fname, disable=DISABLE_REMOVE)
            silent_remove(out_png_fname, disable=DISABLE_REMOVE)
            pass


class TestCheckGauss(unittest.TestCase):
    # These test/demonstrate different options
    def testBasicUse(self):
        silent_remove(SUB_SUB_DIR, dir_with_files=True)
        make_dir(SUB_SUB_DIR)
        fname_after_run = os.path.join(SUB_SUB_DIR, 'for_hartree', 'pet_mono_637_tzvp.log')
        silent_remove(fname_after_run)

        file1_to_copy = os.path.join(SUB_DATA_DIR, 'pet_mono_637_tzvp.tpl')
        temp_fname1 = os.path.join(SUB_SUB_DIR, 'pet_mono_637_tzvp.log')
        file2_to_copy = os.path.join(SUB_DATA_DIR, 'me2propprpnt_7.log')
        temp_fname2 = os.path.join(SUB_SUB_DIR, 'me2propprpnt_7.log')
        file3_to_copy = os.path.join(SUB_DATA_DIR, 'pet_mono_671_tzvp.log')
        temp_fname3 = os.path.join(SUB_SUB_DIR, 'pet_mono_671_tzvp.log')

        good_output = "The following files completed normally:\n" \
                      "    tests/test_data/check_gauss/temp_dir/pet_mono_637_tzvp.log\n" \
                      "The following files may have failed:\n" \
                      "    tests/test_data/check_gauss/temp_dir/me2propprpnt_7.log\n" \
                      "The following files may still be running:\n" \
                      "    tests/test_data/check_gauss/temp_dir/pet_mono_671_tzvp.log\n"

        try:
            copyfile(file1_to_copy, temp_fname1)
            copyfile(file2_to_copy, temp_fname2)
            copyfile(file3_to_copy, temp_fname3)
            test_input = ["-d", SUB_SUB_DIR]
            # main(test_input)
            with capture_stdout(main, test_input) as output:
                self.assertTrue(good_output in output)
            self.assertTrue(os.path.isfile(fname_after_run))
        finally:
            silent_remove(SUB_SUB_DIR, dir_with_files=True, disable=DISABLE_REMOVE)
            pass

    def testSingleFinalConvergence(self):
        test_input = ["-f", SINGLE_FILE, "-z"]
        good_out = 'File Name                            Convergence Convergence_Error\n'\
                   'me2propprpnt_7.log                      111.4981 True\n'
        with capture_stdout(main, test_input) as output:
            self.assertTrue(good_out in output)

    def testListFinalConvergence(self):
        test_input = ["-l", LIST_FILE, "-z"]
        # main(test_input)
        good_out = 'File Name                            Convergence Convergence_Error\n' \
                   'hexyl_acrylate_239.log                    1.1100 False\n' \
                   'hexyl_acrylate_419.log                    0.0706 False\n'
        with capture_stdout(main, test_input) as output:
            self.assertTrue(good_out in output)

    def testListEachStepConvergence(self):
        test_input = ["-l", LIST_FILE, "-s"]
        expected_f_names = [CONV_239_OUT, CONV_239_PNG, CONV_419_OUT, CONV_419_PNG]
        try:
            for f_name in expected_f_names:
                silent_remove(f_name)
            main(test_input)
            self.assertFalse(diff_lines(CONV_239_OUT, GOOD_CONV_239_OUT))
            self.assertFalse(diff_lines(CONV_419_OUT, GOOD_CONV_419_OUT))
            self.assertTrue(os.path.isfile(CONV_239_PNG))
            self.assertTrue(os.path.isfile(CONV_419_PNG))
        finally:
            for f_name in expected_f_names:
                silent_remove(f_name, disable=DISABLE_REMOVE)
            pass

    def testOneEachStepOut(self):
        # tests searching directory with checking convergence, plus using an alternate extension
        test_input = ["-s", "-d", SUB_DATA_DIR, "-e", "ts.out"]
        expected_f_names = [DIOXOLAN_OUT, DIOXOLAN_PNG]
        try:
            for f_name in expected_f_names:
                silent_remove(f_name)
            main(test_input)
            self.assertFalse(diff_lines(DIOXOLAN_OUT, GOOD_DIOXOLAN_OUT))
        finally:
            for f_name in expected_f_names:
                silent_remove(f_name, disable=DISABLE_REMOVE)
            pass

    def testStepConv(self):
        # script reported "Requested convergence data data not found" for this file that is still running
        # this test helped remove this warning
        in_file = os.path.join(SUB_DATA_DIR, "ti_eg5_dime_tpa_tse.log")
        test_input = ["-s", "-f", in_file]
        expected_f_names = [IN_PROGRESS_OUT, IN_PROGRESS_PNG]
        try:
            for f_name in expected_f_names:
                silent_remove(f_name)
            # main(test_input)
            with capture_stderr(main, test_input) as output:
                self.assertTrue("final convergence report" in output)
            for f_name in expected_f_names:
                self.assertTrue(os.path.isfile(f_name))
        finally:
            for f_name in expected_f_names:
                silent_remove(f_name, disable=DISABLE_REMOVE)
            pass

    def testOneEachStopAtStep(self):
        # tests searching directory with checking convergence, plus using an alternate extension
        test_input = ["-t", "37", "-d", SUB_DATA_DIR, "-e", "ts.out"]
        good_output = "Steps sorted by convergence to step number 37 for file: dioxolan4ol_ts4_ts.out\n" \
                      "    StepNum  Convergence\n" \
                      "         35    161.099\n" \
                      "         37    211.219\n" \
                      "          1    342.935\n" \
                      "         36    392.523\n" \
                      "         34    456.271\n"
        # main(test_input)
        with capture_stdout(main, test_input) as output:
            self.assertTrue(good_output in output)

    def testBest10Steps(self):
        # tests searching directory with checking convergence, plus using an alternate extension
        test_input = ["-b", "-d", SUB_DATA_DIR, "-e", "ts.out"]
        good_output = "Best (up to 10) steps sorted by convergence for file: dioxolan4ol_ts4_ts.out\n" \
                      "    StepNum  Convergence\n" \
                      "         43      1.031\n" \
                      "         44      1.146\n" \
                      "         42     42.515\n" \
                      "         41     42.785\n" \
                      "         40     56.686\n" \
                      "         39    102.103\n" \
                      "         35    161.099\n" \
                      "         37    211.219\n" \
                      "         38    252.513\n" \
                      "          1    342.935\n"
        # main(test_input)
        with capture_stdout(main, test_input) as output:
            self.assertTrue(good_output in output)

    def testBest10StepsList(self):
        test_input = ["-b", "-l", LIST_FILE]
        good_output = "Best (up to 10) steps sorted by convergence for file: hexyl_acrylate_239.log\n" \
                      "    StepNum  Convergence\n" \
                      "         17      1.110\n" \
                      "         16      1.138\n" \
                      "          1    749.248\n" \
                      "Best (up to 10) steps sorted by convergence for file: hexyl_acrylate_419.log\n" \
                      "    StepNum  Convergence\n" \
                      "          8      0.071\n" \
                      "          7      4.319\n" \
                      "          1    666.063\n"

        # main(test_input)
        with capture_stdout(main, test_input) as output:
            self.assertTrue(good_output in output)
            pass

    def testAllStepsStdOutList(self):
        test_input = ["-a", "-l", LIST_FILE]
        good_output = "Convergence of all steps for file: hexyl_acrylate_239.log\n" \
                      "    StepNum  Convergence\n" \
                      "          1    749.248\n" \
                      "         16      1.138\n" \
                      "         17      1.110\n" \
                      "Convergence of all steps for file: hexyl_acrylate_419.log\n" \
                      "    StepNum  Convergence\n" \
                      "          1    666.063\n" \
                      "          7      4.319\n" \
                      "          8      0.071\n"
        # main(test_input)
        with capture_stdout(main, test_input) as output:
            print(output)
            self.assertTrue(good_output in output)
            pass

    def testAllStepsNotFromChk(self):
        # if not reading from checkpoint, stoich comes before Coordinates; let's still catch that first step
        test_input = ["-a", "-f", NOT_FROM_CHK_FILE]
        good_output = "Convergence of all steps for file: acyl-min_ts5.out\n" \
                      "    StepNum  Convergence\n" \
                      "          1    120.064\n" \
                      "          2    248.641\n" \
                      "          3    310.534\n"
        # main(test_input)
        with capture_stdout(main, test_input) as output:
            self.assertTrue(good_output in output)
            pass

    def testDirSubdirsLastStep(self):
        temp_files_list = make_fill_sub_dir()
        test_input = ["-ds", SUB_SUB_DIR, "-z"]
        good_output = 'File Name                            Convergence Convergence_Error\n' \
                      'a579.log                                  0.1175 False\n' \
                      'ipvc_11_10_cp.log                            nan None'
        try:
            # main(test_input)
            with capture_stdout(main, test_input) as output:
                self.assertTrue(good_output in output)
        finally:
            for temp_name in temp_files_list:
                silent_remove(temp_name, disable=DISABLE_REMOVE)
            silent_remove(SUB_SUB_DIR, disable=DISABLE_REMOVE, dir_with_files=True)
            pass

    def testDirSubdirsBestSteps(self):
        temp_files_list = make_fill_sub_dir()
        test_input = ["-ds", SUB_SUB_DIR, "-b"]
        good_output = "Best (up to 10) steps sorted by convergence for file: a579.log\n" \
                      "    StepNum  Convergence\n" \
                      "         72      0.117\n" \
                      "          2      3.965\n" \
                      "          1    173.021\n" \
                      "No convergence data found for file: ipvc_11_10_cp.log\n"
        try:
            # main(test_input)
            with capture_stdout(main, test_input) as output:
                self.assertTrue(good_output in output)
        finally:
            for temp_name in temp_files_list:
                silent_remove(temp_name, disable=DISABLE_REMOVE)
            silent_remove(SUB_SUB_DIR, disable=DISABLE_REMOVE, dir_with_files=True)
            pass

    def testCheckFinalConverg(self):
        input_fname = os.path.join(SUB_DATA_DIR, "prop_acetate_8.log")
        test_input = ["-f", input_fname, "-z"]
        # main(test_input)
        with capture_stdout(main, test_input) as output:
            self.assertTrue("False (negligible forces)" in output)

    def testMakeScanPlot(self):
        list_fname = os.path.join(SUB_DATA_DIR, "scan_list.txt")
        fnames = [os.path.join(SUB_DATA_DIR, "pet_dimer_scan_neg_tzvp.log"),
                  os.path.join(SUB_DATA_DIR, "pet_dimer_scan_pos_tzvp.log")]
        list_to_file(fnames, list_fname)
        out_png_fname = os.path.join(SUB_DATA_DIR, "pet_dimer_scan.png")
        silent_remove(out_png_fname)
        test_input = ["-l", list_fname, "--scan", out_png_fname]
        try:
            # main(test_input)
            with capture_stdout(main, test_input) as output:
                self.assertTrue("Barriers" in output)
            self.assertTrue(os.path.isfile(out_png_fname))
        finally:
            silent_remove(list_fname, disable=DISABLE_REMOVE)
            silent_remove(out_png_fname, disable=DISABLE_REMOVE)
            pass


class TestCheckScanPlotParts(unittest.TestCase):
    def testCharmmDihedralFunction(self):
        params = np.asarray([0.8119, 0.2299, -1.7798, -0.8174, 1.53312,
                             88.429, 42.743, 52.9498, 34.9852, 36.352])
        x_pts = np.linspace(0., 360., 50)
        y_pts = charmm_dihedral(x_pts, *params, 1, 2, 3, 4, 6)
        good_y_pts = [3.4592541867020987, 1.9654961492578094, 0.5233786564400287, -0.3551669123072688,
                      -0.5712516976219202, -0.4274046900427555, -0.3832103505493323, -0.7300940072900284,
                      -1.3814831777489025, -1.905592035906774, -1.7854642803433234, -0.7562843579553601,
                      0.9793992506057351, 2.7929346855927433, 3.9182205690838714, 3.8373272531770986,
                      2.553609079056327, 0.6003793408251865, -1.2237247945539975, -2.2522955022449347,
                      -2.268893304902418, -1.5771685790378793, -0.7926596189639956, -0.4783366509178344,
                      -0.8294150631525361, -1.5811889610689298, -2.1854126077851705, -2.148758008951085,
                      -1.336657786568245, -0.06770003116817813, 1.0554479690081546, 1.4817602347536882,
                      1.0221425354975133, -0.04347485400443718, -1.1243475223459654, -1.653191216768191,
                      -1.4062421657029236, -0.6257331766119394, 0.12713040579745627, 0.3072919488465118,
                      -0.28082182345398543, -1.3434755030891519, -2.2363231661216227, -2.310698408850898,
                      -1.2787338410679685, 0.6048162744267835, 2.647055326856186, 4.054050226925558, 4.31979000424518,
                      3.4592541867021036]
        self.assertTrue(np.allclose(y_pts, good_y_pts))

    def testFindBestFit(self):
        # Note: method not currently used in check_gauss; CHARMM dihedral equation fitting would be used to fit a FF
        #    The polynomial fitting is to provide a pretty image with an analytical equation
        x_vals = SCAN_ARRAY[:, 0]
        y_vals = SCAN_ARRAY[:, 1]
        x_pts = np.linspace(0., 360., 50)
        png_out_fname = os.path.join(SUB_DATA_DIR, "test_fit.png")
        out_name_charmm = os.path.join(SUB_DATA_DIR, "test_fit_charmm.png")
        out_name_poly = os.path.join(SUB_DATA_DIR, "test_fit_poly.png")
        out_names = [out_name_charmm, out_name_poly]
        try:
            test_input = [x_vals, y_vals, x_pts, png_out_fname]
            for fname in out_names:
                silent_remove(fname)
            # find_good_fit(*test_input)
            with capture_stdout(find_good_fit, *test_input) as output:
                self.assertTrue("0.44" in output)
                self.assertTrue("0.07" in output)
            for fname in out_names:
                self.assertTrue(os.path.isfile(fname))
        finally:
            for fname in out_names:
                silent_remove(fname, disable=DISABLE_REMOVE)
            pass

    def testMakePlot(self):
        png_out_fname = os.path.join(SUB_DATA_DIR, "scan.png")
        try:
            silent_remove(png_out_fname)
            x_fit, y_fit = plot_scan(SCAN_ARRAY, png_out_fname)
            self.assertTrue(os.path.isfile(png_out_fname))

            # now testing second function:
            find_stable_points(x_fit, y_fit)
            with capture_stdout(find_stable_points, x_fit, y_fit) as output:
                self.assertTrue(" 7.3" in output)
                self.assertTrue(" 3.7" in output)
                self.assertTrue(" 2.7" in output)
        finally:
            silent_remove(png_out_fname, disable=DISABLE_REMOVE)
            pass

    def testProcessScanArray(self):
        scan_data = np.asarray([[180., -1602.04394884], [-170., -1602.04365181], [-160., -1602.04281471],
                                [-150., -1602.04166751], [-140., -1602.04051595], [-130., -1602.03972340],
                                [-120., -1602.03957399], [-110., -1602.04019960], [-100., -1602.04146379],
                                [-90., -1602.04303792], [-80., -1602.04447082], [-70., -1602.04537612],
                                [-60., -1602.04542365], [-50., -1602.04440242], [-40., -1602.04234981],
                                [-30., -1602.03968185], [-20., -1602.03698562], [-10., -1602.03489272],
                                [-0., -1602.03389706]])
        good_processed_data = np.asarray(
            [[180., -1602.04394884], [190., -1602.04365181], [200., -1602.04281471], [210., -1602.04166751],
             [220., -1602.04051595], [230., -1602.03972340], [240., -1602.03957399], [250., -1602.04019960],
             [260., -1602.04146379], [270., -1602.04303792], [280., -1602.04447082], [290., -1602.04537612],
             [300., -1602.04542365], [310., -1602.04440242], [320., -1602.04234981], [330., -1602.03968185],
             [340., -1602.03698562], [350., -1602.03489272], [360., -1602.03389706]])
        first_vals_diff = process_scan_array(scan_data)
        self.assertAlmostEqual(first_vals_diff, 10)
        self.assertTrue(np.allclose(scan_data, good_processed_data))

    def testFindScanStepsSVP(self):
        # this file has optimization at every step, including the first
        scan_file = os.path.join(DATA_DIR, "check_gauss", "pet_dimer_scan_neg.log")
        gausslog_contents = process_gausslog_file(scan_file, collect_scan_steps=True)
        scan_dict = gausslog_contents[SCAN_DICT]
        good_scan_dict = {'179.9999': -1602.04394884,
                          '169.9999': -1602.04365179,
                          '159.9999': -1602.04281467,
                          '149.9999': -1602.04166748,
                          '139.9999': -1602.04051593,
                          '129.9999': -1602.03972339,
                          '119.9998': -1602.03957399,
                          '109.9998': -1602.04019963,
                          '99.9998': -1602.04146383,
                          '89.9998': -1602.04303797,
                          '79.9998': -1602.04447086,
                          '69.9998': -1602.04537614,
                          '59.9999': -1602.04542364,
                          '49.9999': -1602.04440237,
                          '39.9999': -1602.04234975,
                          '29.9999': -1602.03968177,
                          '19.9999': -1602.03698555,
                          '9.9999': -1602.03489267,
                          '-0.0002': -1602.03389705}
        self.assertTrue(len(scan_dict) == len(good_scan_dict))
        for scan_idx, scan_val in scan_dict.items():
            self.assertAlmostEqual(scan_val, good_scan_dict[scan_idx])

    def testFindScanStepsTZVP(self):
        # this file starts with an optimized structure, so only one step for the first step
        scan_file = os.path.join(DATA_DIR, "check_gauss", "pet_dimer_scan_neg_tzvp.log")
        gausslog_contents = process_gausslog_file(scan_file, collect_scan_steps=True)
        scan_dict = gausslog_contents[SCAN_DICT]
        good_scan_dict = {'179.9999': -1603.87812457,
                          '169.9999': -1603.877835,
                          '159.9999': -1603.87700572,
                          '149.9999': -1603.87590461,
                          '139.9999': -1603.87484789,
                          '129.9999': -1603.87417888,
                          '119.9998': -1603.87415467,
                          '109.9998': -1603.87485973,
                          '99.9998': -1603.87614613,
                          '89.9998': -1603.87767359,
                          '79.9998': -1603.87898846,
                          '69.9998': -1603.87971082,
                          '59.9999': -1603.87952259,
                          '49.9999': -1603.87835992,
                          '39.9999': -1603.87628773,
                          '29.9999': -1603.87366776,
                          '19.9999': -1603.87107387,
                          '9.9999': -1603.86908504,
                          '-0.0001': -1603.86813326}
        self.assertTrue(len(scan_dict) == len(good_scan_dict))
        for scan_idx, scan_val in scan_dict.items():
            self.assertAlmostEqual(scan_val, good_scan_dict[scan_idx])

    def testCollectOutputScanSteps(self):
        fnames = [os.path.join(SUB_DATA_DIR, "pet_dimer_scan_neg_tzvp.log"),
                  os.path.join(SUB_DATA_DIR, "pet_dimer_scan_pos_tzvp.log")]
        scan_array = collect_output_scan_steps(fnames)
        good_scan_array = np.asarray([[0., 7.26502855], [9.9999, 6.66777759], [19.9999, 5.41976792],
                                      [29.9999, 3.79207738], [39.9999, 2.14802139], [49.9999, 0.847702544],
                                      [59.9999, 0.118116108], [69.9998, 0.], [79.9998, 0.453287741],
                                      [89.9998, 1.27838112], [99.9998, 2.23687674], [109.9998, 3.04410492],
                                      [119.9998, 3.48653675], [129.9999, 3.47134474], [139.9999, 3.05153463],
                                      [149.9999, 2.38843282], [159.9999, 1.69747587], [169.9999, 1.17709482],
                                      [179.9999, 0.995386898], [189.9998, 1.17708854], [199.9998, 1.69745704],
                                      [209.9998, 2.38841400], [219.9998, 3.05151581], [229.9998, 3.47133847],
                                      [239.9999, 3.48654302], [249.9999, 3.04412374], [259.9999, 2.23690184],
                                      [269.9999, 1.27840622], [279.9999, 0.453312842], [289.9999, 0.00000627502835],
                                      [299.9998, 0.118103557], [309.9998, 0.847671169], [319.9998, 2.14797747],
                                      [329.9998, 3.79202718], [339.9998, 5.41971772], [349.9998, 6.66774621],
                                      [359.9998, 7.26502227]])
        self.assertTrue(np.allclose(scan_array, good_scan_array))
