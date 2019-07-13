#!/usr/bin/env python
"""
makes and runs gaussian input
"""

from __future__ import print_function
import os
import sys
import argparse
import subprocess
import re
from nrel_tools.common import (InvalidDataError, warning, process_cfg, create_out_fname, list_to_file,
                               GOOD_RET, INPUT_ERROR, IO_ERROR, INVALID_DATA,
                               read_tpl)
from nrel_tools.fill_tpl import (OUT_DIR, MAIN_SEC, TPL_VALS_SEC, TPL_EQS_SEC,
                                 TPL_VALS, TPL_EQ_PARAMS, NEW_FNAME, fill_save_tpl)

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
SLURM_NO_CHK_TPL = 'slurm_no_old_chk'
SLURM_FROM_CHK_TPL = 'slurm_from_chk'
OPT_TPL = 'opt_tpl'
STABLE_TPL = 'stable_tpl'

# Defaults
DEF_CFG_FILE = 'run_gauss_bde.ini'
DEF_SLURM_NO_CHK_TPL = 'run_gauss_no_old_chk.tpl'
DEF_SLURM_FROM_CHK_TPL = 'run_gauss_from_old_chk.tpl'
DEF_OPT_TPL = 'opt.tpl'
DEF_STABLE_TPL = 'stable.tpl'


# Set notation
DEF_CFG_VALS = {CONFIG_FILE: DEF_CFG_FILE,
                OUT_DIR: None,
                SLURM_NO_CHK_TPL: DEF_SLURM_NO_CHK_TPL,
                SLURM_FROM_CHK_TPL: DEF_SLURM_FROM_CHK_TPL,
                OPT_TPL: DEF_OPT_TPL,
                STABLE_TPL: DEF_STABLE_TPL,
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
    main_proc = cfg_proc(dict(config.items(MAIN_SEC)), DEF_CFG_VALS, REQ_KEYS)

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
    args = None
    try:
        args = parser.parse_args(argv)
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
        slurm_file_name = create_out_fname(job_name, ext=".sh", base_dir=cfg[OUT_DIR])
        # job_names = ['', '_opt', '_stable']
        job_names = ['', ]
        gau_tpl_files = {'_opt': cfg[OPT_TPL], '_stable': STABLE_TPL}

        # First job, svp
        for job in job_names:
            if job == '':
                tpl_file = cfg[SLURM_NO_CHK_TPL]
            else:
                tpl_file = read_tpl(cfg[SLURM_FROM_CHK_TPL])
                tpl_dict[OLD_JOB_NAME] = tpl_dict[JOB_NAME]
                tpl_dict[JOB_NAME] = job_name + job
                tpl_dict[INPUT_FILE] = gau_tpl_files[job]
            tpl_str = read_tpl(tpl_file)
            fill_save_tpl(cfg, tpl_str, tpl_dict, tpl_file, slurm_file_name)
            subprocess.call(["chmod", "+x", slurm_file_name])
            subprocess.call(slurm_file_name)
            out_file = tpl_dict[JOB_NAME] + ".log"
            last_line = subprocess.check_output(["tail", "-1",  out_file]).strip().decode("utf-8")
            if GAU_GOOD_PAT.match(last_line):
                print("Successfully completed {}".format(out_file))
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
