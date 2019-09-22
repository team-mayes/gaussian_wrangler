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
JOB_LIST = 'job_list'
FOLLOW_JOBS_LIST = 'follow_job_list'
TPL_DICT = 'dictionary of tpls for jobs'
JOB_RUN_TPL = "job_run_tpl"
FIRST_JOB_CHK = 'chk_for_first_job'
OLD_CHECK_ECHO = 'old_check_echo'
# config keys for spawning additional jobs
SBATCH_TPL = 'sbatch_tpl'
INI_TPL = 'ini_tpl'
PARTITION = 'partition'
RUN_TIME = 'run_time'
ACCOUNT = 'account'
EMAIL = 'email'
ALL_NEW = 'all_new'
OPT_OLD_JOB_NAME = 'opt_old_name'
RUN_GAUSS_INI = 'run_gauss_ini'
QOS = 'qos'
LIST_OF_JOBS = 'list_of_jobs'
SETUP_SUBMIT = 'setup_submit'
START_FROM_SAME_CHK = 'start_from_job_name_chk'
NO_SUBMIT = 'no_submit'

DEF_CFG_FILE = 'run_gauss.ini'
DEF_JOB_RUN_TPL = 'run_gauss_job.tpl'
DEF_SLURM_FROM_CHK_TPL = 'run_gauss_from_old_chk.tpl'
DEF_JOB_LIST = ['', 'opt', 'stable', 'freq', 'pfb', 'fb']
DEF_PARTITION = 'short'
DEF_RUN_TIME = '4:00:00'
DEF_ACCOUNT = 'bpms'
DATA_DIR = os.path.join(os.path.dirname(__file__), 'data')
DEF_SBATCH_TPL = os.path.join(DATA_DIR, 'sbatch.tpl')
DEF_INI_TPL = os.path.join(DATA_DIR, 'run_gauss_ini.tpl')
DEF_FOLLOW_JOBS_LIST = None

# Set notation
DEF_CFG_VALS = {OUT_DIR: None,
                JOB_RUN_TPL: DEF_JOB_RUN_TPL,
                JOB_LIST: DEF_JOB_LIST,
                FOLLOW_JOBS_LIST: DEF_FOLLOW_JOBS_LIST,
                FIRST_JOB_CHK: None,
                PARTITION: DEF_PARTITION,
                RUN_TIME: DEF_RUN_TIME,
                ACCOUNT: DEF_ACCOUNT,
                EMAIL: None,
                ALL_NEW: False,
                SBATCH_TPL: DEF_SBATCH_TPL,
                INI_TPL: DEF_INI_TPL,
                QOS: 'normal',
                START_FROM_SAME_CHK: False,
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
    main_proc = cfg_proc(dict(config.items(MAIN_SEC)), DEF_CFG_VALS, REQ_KEYS, int_list=False, store_extra_keys=True)

    main_proc[CONFIG_FILE] = f_loc
    main_proc[TPL_DICT] = {}

    all_job_types = main_proc[JOB_LIST].copy()
    follow_jobs_list = []
    if main_proc[FOLLOW_JOBS_LIST]:
        threads = [thread.split(',') for thread in main_proc[FOLLOW_JOBS_LIST].split(';')]
        for thread in threads:
            thread = [job.strip() for job in thread]
            follow_jobs_list.append(thread)
            all_job_types += thread
        main_proc[FOLLOW_JOBS_LIST] = follow_jobs_list
    else:
        main_proc[FOLLOW_JOBS_LIST] = []

    for job in all_job_types:
        if job == '':
            continue
        if job in main_proc:
            tpl_name = main_proc[job]
        else:
            tpl_name = job + '.tpl'
        if os.path.isfile(tpl_name):
            main_proc[TPL_DICT][job] = tpl_name
        else:
            raise InvalidDataError("For job '{}', could not find a template file '{}'\n"
                                   "You may specify the template to use (including path) using {} as a key in the "
                                   "config file.".format(job, tpl_name, job))

    if not os.path.isfile(main_proc[JOB_RUN_TPL]):
        raise InvalidDataError("Could not find the submit template '{}'".format(main_proc[JOB_RUN_TPL]))

    return main_proc


def parse_cmdline(argv):
    """
    Returns the parsed argument list and return code.
    `argv` is a list of arguments, or `None` for ``sys.argv[1:]``.
    """
    if argv is None:
        argv = sys.argv[1:]

    # initialize the parser object:
    parser = argparse.ArgumentParser(description='Sets up and runs series of Gaussian jobs, checking between jobs '
                                                 'for normal termination.')
    parser.add_argument("job_name", help="The job name to run. If the first job to run is '', a Gaussian input "
                                         "file (with extension '.com') is needed. Otherwise, a checkpoint file "
                                         "(with extension '.chk') is needed.")
    parser.add_argument("-c", "--config", help="The location of the configuration file in ini format. "
                                               "The default file name is {}, located in the base directory "
                                               "where the program as run.".format(DEF_CFG_FILE),
                        default=DEF_CFG_FILE, type=read_cfg)
    parser.add_argument("-n", "--no_submit", help="Set up jobs without submitting them.", action="store_true",
                        default=False)
    parser.add_argument("-o", "--old_chk_file", help="The base name of the checkpoint file (do not include '.chk')"
                                                     "to be used for the first job (optional).", default=None)
    parser.add_argument("-s", "--setup_submit", help="The script will setup and submit, rather than run, the provided "
                                                     "'job_name'. Any extension, or none, can be included in the job "
                                                     "name. If a 'single_job' or 'list_of_jobs' are not specified, "
                                                     "the script will instead attempt to run the 'job_name'. The "
                                                     "default is False.", action="store_true", default=False)
    parser.add_argument("-l", "--list_of_jobs", help="The input in the position of 'job_name' will be read as a file "
                                                     "name with a list of jobs to set up and submit. Each job name "
                                                     "should be on a separate line. Any extension, or none, can follow "
                                                     "the job name. If a 'setup_submit' or 'list_of_jobs' are not "
                                                     "specified, the script will instead attempt to run the 'job_name'."
                                                     " The default is False.", action="store_true", default=False)
    parser.add_argument("-t", "--testing", help="Run in testing mode, which will not check for normal Gaussian "
                                                "termination before continuing. Default is False.",
                        action="store_true", default=False)

    args = None
    try:
        args = parser.parse_args(argv)
        if args.setup_submit and args.list_of_jobs:
            raise InvalidDataError("Cannot choose both 'setup_submit' and 'list_of_jobs' options")
        if args.list_of_jobs:
            if not os.path.isfile(args.job_name):
                raise IOError("When using the 'list_of_jobs' option, the first positional argument ('job_name') must "
                              "be the name of the file with the list of jobs. Could not read: {}".format(args.job_name))
        if args.old_chk_file:
            args.config[FIRST_JOB_CHK] = os.path.splitext(args.old_chk_file)[0]

        # commenting below makes the current directory the default
        # if not args.config[OUT_DIR]:
        #     args.config[OUT_DIR] = os.path.dirname(args.config[CONFIG_FILE])
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


def run_job(job, job_name_perhaps_with_dir, tpl_dict, cfg, testing_flag):
    # Determine if it will run fresh or from an old checkpoint
    if job == '':
        new_job_name = tpl_dict[JOB_NAME]
        tpl_dict[INPUT_FILE] = job_name_perhaps_with_dir + ".com"
        if cfg[FIRST_JOB_CHK]:
            tpl_dict[OLD_CHECK_ECHO] = 'echo "%OldChk={}.chk" >> $INFILE'.format(cfg[FIRST_JOB_CHK])
        else:
            tpl_dict[OLD_CHECK_ECHO] = ''
    else:
        new_job_name = tpl_dict[JOB_NAME] + '_' + job
        tpl_dict[OLD_JOB_NAME] = tpl_dict[JOB_NAME]
        tpl_dict[OLD_CHECK_ECHO] = 'echo "%OldChk={}.chk" >> $INFILE'.format(tpl_dict[OLD_JOB_NAME])
        tpl_dict[INPUT_FILE] = cfg[TPL_DICT][job]
    if not os.path.isfile(tpl_dict[INPUT_FILE]):
        raise IOError(tpl_dict[INPUT_FILE])

    tpl_file = cfg[JOB_RUN_TPL]
    job_runner_fname = create_out_fname(new_job_name, ext=".sh", base_dir=cfg[OUT_DIR])
    print("Running {}".format(new_job_name))

    tpl_dict[JOB_NAME] = new_job_name
    tpl_str = read_tpl(tpl_file)
    fill_save_tpl(cfg, tpl_str, tpl_dict, tpl_file, job_runner_fname)
    subprocess.call(["chmod", "+x", job_runner_fname])
    p1 = subprocess.Popen(job_runner_fname)
    p1.wait()
    if testing_flag:
        print("Testing mode; did not check for normal Gaussian termination.")
    else:
        out_file = tpl_dict[JOB_NAME] + ".log"
        last_line = subprocess.check_output(["tail", "-1",  out_file]).strip().decode("utf-8")
        if GAU_GOOD_PAT.match(last_line):
            print("Successfully completed {}".format(out_file))
            os.remove(job_runner_fname)
        else:
            raise InvalidDataError('Job failed: {}'.format(out_file))


def create_sbatch_dict(cfg, tpl_dict, new_ini_fname, start_from_job_name_chk=True):
    sbatch_dict = {PARTITION: cfg[PARTITION], RUN_TIME: cfg[RUN_TIME], ACCOUNT: cfg[ACCOUNT],
                   JOB_NAME: tpl_dict[JOB_NAME], RUN_GAUSS_INI: new_ini_fname, QOS: cfg[QOS]
                   }

    if cfg[FIRST_JOB_CHK]:
        sbatch_dict[OPT_OLD_JOB_NAME] = '-o ' + cfg[FIRST_JOB_CHK]
    elif start_from_job_name_chk:
        sbatch_dict[OPT_OLD_JOB_NAME] = '-o ' + tpl_dict[JOB_NAME]
    else:
        sbatch_dict[OPT_OLD_JOB_NAME] = ''

    if cfg[EMAIL]:
        sbatch_dict[EMAIL] = '#SBATCH --mail-type=FAIL\n#SBATCH --mail-type=END\n' \
                             '#SBATCH --mail-user={}'.format(cfg[EMAIL])
    else:
        sbatch_dict[EMAIL] = ''

    return sbatch_dict


def add_to_ini(filled_tpl_name, thread, tpl_dict):
    with open(filled_tpl_name, 'a') as f:
        f.write('\n')
        for job in thread:
            if job == '':
                continue
            f.write('{} = {}\n'.format(job, tpl_dict[job]))


def create_ini_dict(cfg, thread):
    ini_dict = {JOB_RUN_TPL: cfg[JOB_RUN_TPL],
                JOB_LIST: ','.join(thread),
                # next are for setting up submits with spawning; doesn't hurt if not needed
                PARTITION: cfg[PARTITION], QOS: cfg[QOS], RUN_TIME: cfg[RUN_TIME], ACCOUNT: cfg[ACCOUNT],
                SBATCH_TPL: cfg[SBATCH_TPL],
                ALL_NEW: cfg[ALL_NEW],
                INI_TPL: cfg[INI_TPL],
                OUT_DIR: cfg[OUT_DIR],
                FIRST_JOB_CHK: cfg[FIRST_JOB_CHK],
                EMAIL: cfg[EMAIL],
                }
    jobs_string = None
    if len(cfg[FOLLOW_JOBS_LIST]) > 0 and (cfg[SETUP_SUBMIT] or cfg[LIST_OF_JOBS]):
        thread_list = []
        for thread in cfg[FOLLOW_JOBS_LIST]:
            thread_list.append(','.join(thread))
        jobs_string = ';'.join(thread_list)
    ini_dict[FOLLOW_JOBS_LIST] = jobs_string

    return ini_dict


def setup_and_submit(cfg, index, thread, tpl_dict):
    ini_dict = create_ini_dict(cfg, thread)
    tpl_str = read_tpl(cfg[INI_TPL])
    if cfg[SETUP_SUBMIT]:
        base_name = tpl_dict[JOB_NAME]
    elif cfg[LIST_OF_JOBS]:
        base_name = tpl_dict[JOB_NAME]
    else:
        base_name = cfg[CONFIG_FILE]
    new_ini_fname = create_out_fname(base_name, suffix=str(index), ext='.ini', base_dir=cfg[OUT_DIR])
    fill_save_tpl(cfg, tpl_str, ini_dict, cfg[INI_TPL], new_ini_fname)
    add_to_ini(new_ini_fname, thread, cfg[TPL_DICT])
    tpl_str = read_tpl(cfg[SBATCH_TPL])
    sbatch_dict = create_sbatch_dict(cfg, tpl_dict, new_ini_fname, start_from_job_name_chk=cfg[START_FROM_SAME_CHK])
    new_sbatch_fname = create_out_fname(base_name, suffix=str(index), ext='.slurm', base_dir=cfg[OUT_DIR])
    fill_save_tpl(cfg, tpl_str, sbatch_dict, cfg[SBATCH_TPL], new_sbatch_fname)
    if not cfg[NO_SUBMIT]:
        try:
            sbatch_result = subprocess.check_output(["sbatch", new_sbatch_fname]).strip().decode("utf-8")
            print(sbatch_result)
        except IOError as e:
            print(e)


def main(argv=None):
    # Read input
    args, ret = parse_cmdline(argv)
    if ret != GOOD_RET or args is None:
        return ret

    cfg = args.config
    # to keep these values handy
    cfg[LIST_OF_JOBS] = args.list_of_jobs
    cfg[SETUP_SUBMIT] = args.setup_submit
    cfg[NO_SUBMIT] = args.no_submit

    # Read template and data files
    try:
        if args.list_of_jobs:
            with open(args.job_name) as f:
                for line in f:
                    base_name = os.path.splitext(os.path.basename(line.strip()))[0]
                    tpl_dict = {JOB_NAME: base_name}
                    setup_and_submit(cfg, '', cfg[JOB_LIST], tpl_dict)
            return GOOD_RET

        job_name_perhaps_with_dir = args.job_name
        job_name = os.path.splitext(os.path.basename(args.job_name))[0]
        tpl_dict = {JOB_NAME: job_name}

        if args.setup_submit:
            setup_and_submit(cfg, '', cfg[JOB_LIST], tpl_dict)
            return GOOD_RET

        for job in cfg[JOB_LIST]:
            run_job(job, job_name_perhaps_with_dir, tpl_dict, cfg, args.testing)

        if len(cfg[FOLLOW_JOBS_LIST]) > 1:
            for index, thread in enumerate(cfg[FOLLOW_JOBS_LIST]):
                if index == 0 and not cfg[ALL_NEW]:
                    continue
                setup_and_submit(cfg, index, thread, tpl_dict)

        if len(cfg[FOLLOW_JOBS_LIST]) > 0 and not cfg[ALL_NEW]:
            for job in cfg[FOLLOW_JOBS_LIST][0]:
                run_job(job, job_name_perhaps_with_dir, tpl_dict, cfg, args.testing)

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
