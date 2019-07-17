#!/usr/bin/env python
"""
makes and runs gaussian input
"""

from __future__ import print_function
import sys
import argparse
import subprocess
import re
import os
from nrel_tools.common import (InvalidDataError, warning, process_cfg, create_out_fname,
                               GOOD_RET, INPUT_ERROR, IO_ERROR, INVALID_DATA,
                               read_tpl)
from nrel_tools.fill_tpl import (OUT_DIR, MAIN_SEC, fill_save_tpl)

try:
    # noinspection PyCompatibility
    from ConfigParser import ConfigParser, MissingSectionHeaderError
except ImportError:
    # noinspection PyCompatibility
    from configparser import ConfigParser, MissingSectionHeaderError

__author__ = 'hmayes'


# Constants #

# Config keys
CONFIG_FILE = 'config_file_name'
# OPT_TPL = 'opt_tpl'
# STABLE_TPL = 'stable_tpl'
# FREQ_TPL = 'freq_tpl'
# PFB_TPL = 'pfb_tpl'
# FB_TPL = 'fb_tpl'
JOB_LIST = 'job_list'
TPL_DICT = 'dictionary of tpls for jobs'
FIRST_SUBMIT_TPL = "first_submit_tpl"
REMAINING_JOB_SUBMIT_TPL = "rest_submit_tpl"
FIRST_JOB_CHK = 'chk_for_first_job'

DEF_CFG_FILE = 'run_gauss_bde.ini'
DEF_SLURM_NO_CHK_TPL = 'run_gauss_no_old_chk.tpl'
DEF_SLURM_FROM_CHK_TPL = 'run_gauss_from_old_chk.tpl'
# DEF_OPT_TPL = 'opt.tpl'
# DEF_STABLE_TPL = 'stable.tpl'
# DEF_FREQ_TPL = 'freq.tpl'
# DEF_PFB_TPL = 'pfb.tpl'
# DEF_FB_TPL = 'fb.tpl'
DEF_JOB_LIST = ['', 'opt', 'stable', 'freq', 'pfb', 'fb']

# Set notation
DEF_CFG_VALS = {CONFIG_FILE: DEF_CFG_FILE,
                OUT_DIR: None,
                FIRST_SUBMIT_TPL: DEF_SLURM_NO_CHK_TPL,
                REMAINING_JOB_SUBMIT_TPL: DEF_SLURM_FROM_CHK_TPL,
                # OPT_TPL: DEF_OPT_TPL,
                # STABLE_TPL: DEF_STABLE_TPL,
                # FREQ_TPL: DEF_FREQ_TPL,
                # PFB_TPL: DEF_PFB_TPL,
                # FB_TPL: DEF_FB_TPL,
                JOB_LIST: DEF_JOB_LIST,
                FIRST_JOB_CHK: None,
                }
REQ_KEYS = {
            }

JOB_NAME = 'job_name'
OLD_JOB_NAME = 'old_job_name'
INPUT_FILE = 'input_file'
GAU_GOOD_PAT = re.compile(r"Normal termination of Gaussian.*")


def read_cfg(f_loc, cfg_proc=process_cfg):
    """
    Reads the given configuration file, returning a dict with the converted values supplemented by default values.

    :param f_loc: The location of the file to read.
    :param cfg_proc: The processor to use for the raw configuration values.  Uses default values when the raw
        value is missing.
    :return: A dict of the processed configuration file's data.
    """
    config = ConfigParser()
    good_files = config.read(f_loc)

    if not good_files:
        raise IOError('Could not read file {}'.format(f_loc))
    main_proc = cfg_proc(dict(config.items(MAIN_SEC)), DEF_CFG_VALS, REQ_KEYS, int_list=False)

    for job in main_proc[JOB_LIST]:
        if job == '':
            continue
        tpl_name = job + '.tpl'
        if os.path.isfile(tpl_name):
            main_proc[TPL_DICT][job] = tpl_name
        else:
            raise InvalidDataError("For job '{}', could not find a template file '{}'".format(job, tpl_name))

    for file in [main_proc[REMAINING_JOB_SUBMIT_TPL], main_proc[FIRST_SUBMIT_TPL]]:
        if not os.path.isfile(file):
            raise InvalidDataError("Could not find the submit template '{}'".format(file))

    return main_proc


def parse_cmdline(argv):
    """
    Returns the parsed argument list and return code.
    `argv` is a list of arguments, or `None` for ``sys.argv[1:]``.
    """
    if argv is None:
        argv = sys.argv[1:]

    # initialize the parser object:
    parser = argparse.ArgumentParser(description='Runs the series of Gaussian jobs for BDEs using Slurm.')
    parser.add_argument("job_name", help="The job name to run")
    parser.add_argument("-c", "--config", help="The location of the configuration file in ini format. "
                                               "The default file name is {}, located in the "
                                               "base directory where the program as run.".format(DEF_CFG_FILE),
                        default=DEF_CFG_FILE, type=read_cfg)
    parser.add_argument("-o", "--old_chk_file", help="The base name of the checkpoint file (do not include '.chk')"
                                                     "to be used for the first job (optional).",
                        default=None)
    args = None
    try:
        args = parser.parse_args(argv)
        if args.old_chk_file:
            args.config[FIRST_JOB_CHK] = args.old_chk_file
    except IOError as e:
        warning("Problems reading file:", e)
        parser.print_help()
        return args, IO_ERROR
    except (KeyError, InvalidDataError, MissingSectionHeaderError, SystemExit) as e:
        if hasattr(e, 'code') and e.code == 0:
            return args, GOOD_RET
        warning(e)
        parser.print_help()
        return args, INPUT_ERROR

    return args, GOOD_RET


def main(argv=None):
    # Read input
    args, ret = parse_cmdline(argv)
    if ret != GOOD_RET or args is None:
        return ret

    cfg = args.config
    job_name = args.job_name

    # Read template and data files
    try:
        tpl_dict = {JOB_NAME: job_name}
        # gau_tpl_files = {'_opt': cfg[OPT_TPL], '_stable': cfg[STABLE_TPL], '_freq': cfg[FREQ_TPL],
        #                  '_pfb': cfg[PFB_TPL], '_fb': cfg[FB_TPL]}

        # First job, svp
        for job in cfg[JOB_LIST]:
            if job == '':
                new_job_name = tpl_dict[JOB_NAME]
                tpl_dict[INPUT_FILE] = job_name + ".com"
                if cfg[FIRST_JOB_CHK]:
                    tpl_dict[OLD_JOB_NAME] = cfg[FIRST_JOB_CHK]
                    tpl_file = cfg[REMAINING_JOB_SUBMIT_TPL]
                else:
                    tpl_file = cfg[FIRST_SUBMIT_TPL]
            else:
                new_job_name = tpl_dict[JOB_NAME] + '_' + job
                tpl_file = cfg[REMAINING_JOB_SUBMIT_TPL]
                tpl_dict[OLD_JOB_NAME] = tpl_dict[JOB_NAME]
                tpl_dict[INPUT_FILE] = cfg[TPL_DICT][job]
            filled_tpl_name = create_out_fname(new_job_name, ext=".sh", base_dir=cfg[OUT_DIR])
            print("Running {}".format(new_job_name))

            tpl_dict[JOB_NAME] = new_job_name
            tpl_str = read_tpl(tpl_file)
            fill_save_tpl(cfg, tpl_str, tpl_dict, tpl_file, filled_tpl_name)
            subprocess.call(["chmod", "+x", filled_tpl_name])
            p1 = subprocess.Popen(filled_tpl_name)
            p1.wait()
            out_file = tpl_dict[JOB_NAME] + ".log"
            last_line = subprocess.check_output(["tail", "-1",  out_file]).strip().decode("utf-8")
            if GAU_GOOD_PAT.match(last_line):
                print("Successfully completed {}".format(out_file))
                # os.remove('infile_' + job_name)
                os.remove(filled_tpl_name)
            else:
                return INVALID_DATA

    except IOError as e:
        warning("Problems reading file:", e)
        return IO_ERROR
    except InvalidDataError as e:
        warning("Problems reading data:", e)
        return INVALID_DATA

    return GOOD_RET  # success


if __name__ == '__main__':
    status = main()
    sys.exit(status)
