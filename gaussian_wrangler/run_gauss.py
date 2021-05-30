#!/usr/bin/env python
"""
makes and runs gaussian input
"""

import sys
import argparse
import subprocess
import re
import os
from configparser import ConfigParser, MissingSectionHeaderError
from common_wrangler.common import (GOOD_RET, INPUT_ERROR, IO_ERROR, INVALID_DATA, OUT_DIR, MAIN_SEC,
                                    InvalidInputError, InvalidDataError, warning,
                                    create_out_fname, get_fname_root, list_to_file, process_cfg, read_tpl, str_to_file)
from common_wrangler.fill_tpl import fill_save_tpl
from gaussian_wrangler.gw_common import GAU_HEADER_PAT
from gaussian_wrangler import __version__

__author__ = 'hmayes'


# Constants #

# For checking unix node memory and cpu, specifically for this script
SOURCE_DIR = os.path.dirname(__file__)
DATA_DIR = os.path.join(SOURCE_DIR, 'data')

MEMINFO_FILE = os.path.join(DATA_DIR, 'meminfo')
CPUINFO_FILE = os.path.join(DATA_DIR, 'cpuinfo')
DF_H = os.path.join(DATA_DIR, 'df_h')

MEM_TOT_PAT = re.compile(r"MemTotal.*")
MEM_FREE_PAT = re.compile(r"MemFree.*")
CACHE_PAT = re.compile(r"cache size.*")
PROC_PAT = re.compile(r"processor.*")

SCRATCH_DIR = 'scratch_dir'
DEF_ROUTE = 'default_route'

# Config keys
CONFIG_FILE = 'config_file_name'
JOB_LIST = 'job_list'
FOLLOW_JOBS_LIST = 'follow_job_list'
GAUSS_IN_EXT = 'gaussian_input_ext'
TPL_DICT = 'dictionary of tpls for jobs'
JOB_RUN_TPL = "job_run_tpl"
FIRST_JOB_CHK = 'chk_for_first_job'
OLD_CHECK_ECHO = 'old_check_echo'
# config keys for spawning additional jobs
SBATCH_TPL = 'sbatch_tpl'
JOB_DESCRIP = 'job_descrip'
PARTITION = 'partition'
RUN_TIME = 'run_time'
ACCOUNT = 'account'
EMAIL = 'email'
USER = 'user'
MEM = 'mem'
PROC_LIST = 'proc_list'
ALL_NEW = 'all_new'
CHK_BASENAME = 'chk_base_name'
READ_CHK_OPTION = 'read_chk_option'
RUN_GAUSS_INI = 'run_gauss_ini'
QOS = 'qos'
LIST_OF_JOBS = 'list_of_jobs'
SETUP_SUBMIT = 'setup_submit'
START_FROM_SAME_CHK = 'start_from_job_name_chk'
NO_SUBMIT = 'no_submit'
CHECK_FOR_CHK = "check_for_chk"
CHK_EXT = '.chk'
KEYS_FOR_SPAWNING_SBATCH = [JOB_RUN_TPL, PARTITION, QOS, RUN_TIME, ACCOUNT, SBATCH_TPL, EMAIL, ALL_NEW,
                            USER, PROC_LIST, MEM]
KEYS_FOR_SPAWNING_INIS = [USER, PROC_LIST, MEM, FIRST_JOB_CHK, OLD_CHECK_ECHO]

DEF_CFG_FILE = 'run_gauss.ini'
DEF_JOB_RUN_TPL = 'run_gauss_job.tpl'
DEF_GAUSS_IN_EXT = '.com'
DEF_JOB_LIST = None
DEF_PARTITION = 'short'
DEF_RUN_TIME = '4:00:00'
DEF_ACCOUNT = 'bpms'
DEF_SBATCH_TPL = os.path.join(DATA_DIR, 'sbatch.tpl')
DEF_FOLLOW_JOBS_LIST = None
DEF_OLD_CHK_STR = 'echo "%OldChk={}.chk" >> ${{INFILE}}'

# Set notation
DEF_CFG_VALS = {OUT_DIR: None,
                GAUSS_IN_EXT: DEF_GAUSS_IN_EXT,
                JOB_RUN_TPL: DEF_JOB_RUN_TPL,
                JOB_LIST: None,
                FOLLOW_JOBS_LIST: DEF_FOLLOW_JOBS_LIST,
                FIRST_JOB_CHK: None,
                PARTITION: DEF_PARTITION,
                RUN_TIME: DEF_RUN_TIME,
                ACCOUNT: DEF_ACCOUNT,
                EMAIL: None,
                ALL_NEW: False,
                SBATCH_TPL: DEF_SBATCH_TPL,
                QOS: 'normal',
                START_FROM_SAME_CHK: False,
                OLD_CHECK_ECHO: DEF_OLD_CHK_STR,
                CHECK_FOR_CHK: True,
                USER: None,
                PROC_LIST: None,
                MEM: None,
                SETUP_SUBMIT: False,
                LIST_OF_JOBS: False,
                SCRATCH_DIR: None,
                }
REQ_KEYS = {
            }

JOB_NAME = 'job_name'
OLD_JOB_NAME = 'old_job_name'
INPUT_FILE = 'input_file'
GAU_GOOD_PAT = re.compile(r"Normal termination of Gaussian.*")
GUESS_READ_OR_GEOM_CHK_PAT = re.compile(r"^.*\b(guess.*read|geom.*check)\b.*$", re.I)


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

    all_job_types = []
    for job_list_key in [JOB_LIST, FOLLOW_JOBS_LIST]:
        job_list = main_proc[job_list_key]
        jobs_list = []
        # need to distinguish between NoneType and ''
        if job_list is not None:
            threads = [thread.split(',') for thread in job_list.split(';')]
            for thread in threads:
                thread = [job.strip() for job in thread]
                jobs_list.append(thread)
                all_job_types += thread
            main_proc[job_list_key] = jobs_list
        else:
            main_proc[job_list_key] = []

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
    parser.add_argument("job_name", help="The job name to run. If the first job to run is '', a Gaussian input file "
                                         "(with extension '{}' or specified with '{}' argument in the config file) is "
                                         "needed. Otherwise, a checkpoint file (with extension '.chk') is "
                                         "needed.".format(DEF_GAUSS_IN_EXT, GAUSS_IN_EXT))
    parser.add_argument("-c", "--config", help="The location of the configuration file in ini format. "
                                               "The default file name is {}, located in the base directory "
                                               "where the program as run.".format(DEF_CFG_FILE),
                        default=DEF_CFG_FILE, type=read_cfg)
    parser.add_argument("-i", "--ignore_chk_warning", help="Ignore warning that a chk file cannot be found in the "
                                                           "current directory for a job that will attempt to read it. "
                                                           "Default is False.",
                        action="store_true", default=False)
    parser.add_argument("-l", "--list_of_jobs", help="The input in the position of 'job_name' will be read as a file "
                                                     "name with a list of jobs to set up and submit. Each job name "
                                                     "should be on a separate line. Any extension, or none, can follow "
                                                     "the job name. If a 'setup_submit' or 'list_of_jobs' are not "
                                                     "specified, the script will instead attempt to run the 'job_name'."
                                                     " The default is False.", action="store_true", default=False)
    parser.add_argument("-n", "--no_submit", help="Set up jobs without submitting them. This flag only effects the "
                                                  "'-s' and '-l' options.", action="store_true", default=False)
    parser.add_argument("-o", "--old_chk_fname", help="The name of the checkpoint file (will use base name plus "
                                                      "'.chk' whether or not an extension of any type is provided) "
                                                      "to be used for the first job (optional).", default=None)
    parser.add_argument("-s", "--setup_submit", help="The script will setup and submit, rather than run, the provided "
                                                     "'job_name'. Any extension, or none, can be included in the job "
                                                     "name. If a 'single_job' or 'list_of_jobs' are not specified, "
                                                     "the script will instead attempt to run the 'job_name'. The "
                                                     "default is False.", action="store_true", default=False)
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
                raise IOError("When using the 'list_of_jobs' option, the first positional argument \n    ('job_name') "
                              "must be the name of the file with the list of jobs. "
                              "Could not read: {}".format(args.job_name))

        if not (args.list_of_jobs or args.setup_submit):
            if len(args.config[JOB_LIST]) > 1:
                raise InvalidDataError("Found ';' in the '{}'. This option (setting up multiple job threads) is "
                                       "currently only supported for setting up (and optionally submitting) jobs "
                                       "(using the '-s' or '-l' options).".format(JOB_LIST))
            elif len(args.config[JOB_LIST]) == 1:
                args.config[JOB_LIST] = args.config[JOB_LIST][0]

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


def get_proc_info(testing_mode):
    """
    Read standard linux /proc/cpuinfo file, and return list of procs in format Gaussian wants
    # :param hostname: string, used to report where the information was retrieved
    :param testing_mode: flag to not actually look for the running machine's /proc/cpuinfo during testing,
              as this would fail for non-unix machines, and give different results for different machines.
              Instead, it will read a file that is part of the package, which is a copy of a file from a unix machine
    :return: string: list of procs in format Gaussian wants, and max_cache according to Gaussian's algorithm
    """
    if testing_mode:
        cpu_info_file = CPUINFO_FILE
    else:
        # as noted above, do not want to run this in testing mode, so will not be covered
        cpu_info_file = "/proc/cpuinfo"
    # minimize subprocess calls specific processor info, and being extra careful to check obtained data
    raw_proc_info = subprocess.check_output(["cat", cpu_info_file]).decode("utf-8").split('\n')

    raw_proc_list = []
    # Gaussian default (conservative) is 1024 * 1024
    cache_size = 1024
    for line in raw_proc_info:
        if PROC_PAT.match(line) or CACHE_PAT.match(line):
            split_line = line.split()
            if CACHE_PAT.match(line):
                # fine that it overwrites every time; would be weird if different
                assert split_line[-1] == 'KB'
                cache_size = int(split_line[-2])
            else:
                raw_proc_list.append(int(split_line[-1]))

    num_procs = len(raw_proc_list)
    first_proc = raw_proc_list[0]
    last_proc = raw_proc_list[-1]

    assert (last_proc - first_proc + 1) == num_procs
    proc_list = "{}-{}".format(first_proc, last_proc)

    max_cache = (cache_size * 1024) / num_procs

    return num_procs, proc_list, max_cache


def get_node_mem(testing_mode):
    """
    Read standard linux /proc/meminfo file, and return list of procs in format Gaussian wants
    # :param hostname: string, used to report where the information was retrieved
    :param testing_mode: flag to not actually look for the running machine's /proc/cpuinfo during testing,
              as this would fail for non-unix machines, and give different results for different machines.
              Instead, it will read a file that is part of the package, which is a copy of a file from a unix machine
    :return: string: maximum memory that Gaussian may allocation, in form Gaussian wants
    """
    if testing_mode:
        mem_info_file = MEMINFO_FILE
    else:
        # as noted above, do not want to run this in testing mode, so will not be covered
        mem_info_file = "/proc/meminfo"

    # To make IDE happy, and check that both pieces of info have been obtained
    mem_tot, mem_free = 0, 0

    raw_mem_info = subprocess.check_output(["cat", mem_info_file]).decode("utf-8").split('\n')
    # not relying on expected order, even though likely okay
    for line in raw_mem_info:
        if MEM_TOT_PAT.match(line) or MEM_FREE_PAT.match(line):
            split_line = line.split()
            assert split_line[-1] == 'kB'
            if MEM_TOT_PAT.match(line):
                mem_tot = int(split_line[1])
            else:
                mem_free = int(split_line[1])
            if mem_tot and mem_free:
                # no other info needed, so stop iteration
                break
    mem_alloc = int(min(mem_tot * .75, mem_free * .85))
    print("    Found {} kB total memory, and {} kB available memory.\n    Will allocate up to {} kB (the lessor of 75% "
          "of MemTotal or 85% of MemFree).\n".format(mem_tot, mem_free, mem_alloc))
    return "{}KB".format(mem_alloc)


def get_max_disk(testing_mode):
    """
    Read standard linux /proc/meminfo file, and return list of procs in format Gaussian wants
    :param testing_mode: flag to not actually look for the running machine's /proc/cpuinfo during testing,
              as this would fail for non-unix machines, and give different results for different machines.
              Instead, it will read a file that is part of the package, which is a copy of a file from a unix machine
    :return: string: maximum disk space that Gaussian may use, in form Gaussian wants
    """
    if testing_mode:
        raw_disk_info = subprocess.check_output(["cat", DF_H]).decode("utf-8").split('\n')
    else:
        # as noted above, do not want to run this in testing mode, so will not be covered
        raw_disk_info = subprocess.check_output(["df", "-h"]).decode("utf-8").split('\n')

    # this is conservative for Eagle and Comet
    raw_avail = '6G'
    for line_num, line in enumerate(raw_disk_info):
        split_line = line.split()
        if line_num == 0:
            assert split_line[-4] == "Avail"
            assert split_line[-2] == "Mounted"
        else:
            if split_line[-1] == '/':
                raw_avail = split_line[-3]
                break
    avail_list = re.split(r'([\d.]+)(\w+)', raw_avail)
    max_disk_num = 0.9 * float(avail_list[1])
    max_disk_unit = avail_list[2]
    max_disk = "{:.2f}{}".format(max_disk_num, max_disk_unit)

    return max_disk


def run_job(job, job_name_perhaps_with_dir, tpl_dict, cfg, testing_mode):
    # Determine if it will run fresh or from an old checkpoint
    if job == '':
        new_job_name = tpl_dict[JOB_NAME]
        tpl_dict[INPUT_FILE] = job_name_perhaps_with_dir + cfg[GAUSS_IN_EXT]
        if cfg[FIRST_JOB_CHK]:
            tpl_dict[OLD_CHECK_ECHO] = cfg[OLD_CHECK_ECHO].format(cfg[FIRST_JOB_CHK])
        else:
            tpl_dict[OLD_CHECK_ECHO] = ''
    else:
        new_job_name = tpl_dict[JOB_NAME] + '_' + job
        tpl_dict[OLD_JOB_NAME] = tpl_dict[JOB_NAME]
        tpl_dict[OLD_CHECK_ECHO] = cfg[OLD_CHECK_ECHO].format(tpl_dict[OLD_JOB_NAME])
        tpl_dict[INPUT_FILE] = cfg[TPL_DICT][job]

    tpl_file = cfg[JOB_RUN_TPL]
    job_runner_fname = create_out_fname(new_job_name, ext=".sh", base_dir=cfg[OUT_DIR])
    print("Running {}".format(new_job_name))

    tpl_dict[JOB_NAME] = new_job_name
    for key_name in [USER, MEM, PROC_LIST, ]:
        if key_name in cfg:
            tpl_dict[key_name] = cfg[key_name]

    tpl_str = read_tpl(tpl_file)
    # if either MEM or PROC_LIST is the default (Nonetype), and is used to run the job, get info from the node before
    #    creating the job script
    mem_required = '{' + MEM + '}' in tpl_str
    get_mem = mem_required and not tpl_dict[MEM]

    proc_required = '{' + PROC_LIST + '}' in tpl_str
    get_proc = proc_required and not tpl_dict[PROC_LIST]

    default_gauss_required = '{' + DEF_ROUTE + '}' in tpl_str

    num_procs = 1  # to make IDE happy
    proc_list = '0'  # to make IDE happy

    if get_mem or get_proc or default_gauss_required:
        # explicitly check each possible required info flag, because any or all can be requested
        if testing_mode:
            hostname = subprocess.check_output(["echo", "r1i7n35"]).decode("utf-8").strip()
        else:
            #  Will not be covered in testing mode, as is not part of written code to be tested
            hostname = subprocess.check_output(["hostname"]).decode("utf-8").strip()
        print("Obtaining available memory and/or number of processors on node {}.\n    "
              "Note: this program assumes the whole node will be allocated to Gaussian.\n".format(hostname))
        if get_mem:
            tpl_dict[MEM] = get_node_mem(testing_mode)

        max_cache = 1024 * 1024  # to make IDE happy; Gaussian default (conservative) is 1024 * 1024
        if get_proc or default_gauss_required:
            num_procs, proc_list, max_cache = get_proc_info(testing_mode)
        if get_proc:
            tpl_dict[PROC_LIST] = proc_list
            print("    Found {} processors. Will allow use of cpus {}.\n".format(num_procs, proc_list))

        if get_mem or get_proc:
            print("    The user may override these values by specifying the '{}' and/or '{}' keywords in the "
                  "configuration file.\n    Be sure to use the formatting Gaussian expects.\n".format(MEM, PROC_LIST))

        if default_gauss_required:
            max_disk = get_max_disk(testing_mode)
            max_cache = int(max_cache)
            print("Since '{}' found in the {}, read machine specs to determine CacheSize={} and "
                  "MaxDisk={}".format(DEF_ROUTE, JOB_RUN_TPL, max_cache, max_disk))
            default_route_list = ["-#- CacheSize={}".format(max_cache), "-#- MaxDisk={}".format(max_disk)]
            fname = create_out_fname('Default.Route', base_dir=cfg[SCRATCH_DIR])
            list_to_file(default_route_list, fname)
            tpl_dict[DEF_ROUTE] = ''   # there is an action triggered, not a value needed, so replaced with blank space

    move_on = False
    while not move_on:
        try:
            fill_save_tpl(tpl_str, tpl_dict, tpl_file, job_runner_fname)
            move_on = True
        except KeyError as e:
            missing_key = e.args[0].split("\'")[1]
            if missing_key in cfg:
                tpl_dict[missing_key] = cfg[missing_key]
            else:
                raise e
    subprocess.call(["chmod", "+x", job_runner_fname])
    if testing_mode:
        print("Testing mode; did not run job script or check Gaussian output for normal termination.\n")
    else:
        # do not want this tested, as actually running Gaussian would take too long, and not what should be tested
        p1 = subprocess.Popen(job_runner_fname)
        p1.wait()
        out_file = tpl_dict[JOB_NAME] + ".log"
        last_line = subprocess.check_output(["tail", "-1",  out_file]).strip().decode("utf-8")
        if GAU_GOOD_PAT.match(last_line):
            print("Successfully completed {}\n".format(out_file))
            os.remove(job_runner_fname)
        else:
            raise InvalidDataError('Job failed: {}'.format(out_file))


def create_sbatch_dict(cfg, tpl_dict, new_ini_fname, current_job_list, start_from_job_name_chk=True,
                       ignore_chk_warning=False):
    sbatch_dict = {PARTITION: cfg[PARTITION], RUN_TIME: cfg[RUN_TIME], ACCOUNT: cfg[ACCOUNT],
                   JOB_NAME: tpl_dict[JOB_NAME], RUN_GAUSS_INI: new_ini_fname, QOS: cfg[QOS],
                   JOB_DESCRIP: tpl_dict[JOB_DESCRIP],
                   }

    if cfg[FIRST_JOB_CHK]:
        if not os.path.isfile(cfg[FIRST_JOB_CHK] + CHK_EXT):
            raise InvalidInputError("Could not find specified '{}': {}".format(FIRST_JOB_CHK,
                                                                               cfg[FIRST_JOB_CHK] + CHK_EXT))
        sbatch_dict[OLD_CHECK_ECHO] = '-o ' + cfg[FIRST_JOB_CHK]
    elif start_from_job_name_chk:
        fname_to_check = tpl_dict[JOB_NAME] + CHK_EXT
        if not os.path.isfile(fname_to_check):
            raise InvalidDataError("Could not find required checkpoint file: {}".format(fname_to_check))
        sbatch_dict[OLD_CHECK_ECHO] = '-o ' + tpl_dict[JOB_NAME]
    else:
        sbatch_dict[OLD_CHECK_ECHO] = ''
        if current_job_list[0] == '' and cfg[CHECK_FOR_CHK]:
            # in the case when there is no old_check_file, make sure the first input file does not try to read from chk
            # IOError is already caught; no don't need to add a try loop
            with open(tpl_dict[INPUT_FILE]) as f:
                try:
                    read_route = False
                    for line in f:
                        line = line.strip()
                        # route can be multiple lines, so first fine the line, then continue until a blank is reached
                        if GAU_HEADER_PAT.match(line):
                            read_route = True
                            while line != '':
                                if GUESS_READ_OR_GEOM_CHK_PAT.match(line) and not ignore_chk_warning:
                                    raise InvalidDataError("Did not find an old checkpoint file to read, but the "
                                                           "Gaussian input header indicates that Gaussian will attempt "
                                                           "and fail to read from a checkpoint:\n   file:  {}\n"
                                                           "  route:  {} ".format(tpl_dict[INPUT_FILE], line))
                                line = next(f).strip()
                    if not read_route:
                        raise StopIteration
                except StopIteration:
                    raise InvalidDataError('The specified input file does not appear valid: {}'
                                           ''.format(tpl_dict[INPUT_FILE]))

    if cfg[EMAIL]:
        sbatch_dict[EMAIL] = '#SBATCH --mail-type=FAIL\n#SBATCH --mail-type=END\n' \
                             '#SBATCH --mail-user={}'.format(cfg[EMAIL])
    else:
        sbatch_dict[EMAIL] = ''

    return sbatch_dict


def create_ini_with_req_keys(thread, tpl_dict, cfg, new_ini_fname):
    """
    Adds to the ini_tpl if a non-default parameter needs to be specified, and is not already included
    :param thread: the jobs to be run, to check that if a location of the tpl is specified, it gets included
                   in the new ini file
    :param tpl_dict: dictionary of locations of job tpl files
    :param new_ini_fname: name and location of the ini_tpl created
    :param cfg: configuration dict
    :return: tpl_str: str with values to be filled
    """
    # This is the minimum needed; more added later if needed
    job_list_string = ','.join(thread)
    # we always need the header, job_list, and job_run_tpl
    tpl_str = '[main]\njob_run_tpl = {}\njob_list = {}'.format(cfg[JOB_RUN_TPL], job_list_string, )

    if len(cfg[FOLLOW_JOBS_LIST]) > 0 and (cfg[SETUP_SUBMIT] or cfg[LIST_OF_JOBS]):
        thread_list = []
        for thread in cfg[FOLLOW_JOBS_LIST]:
            thread_list.append(','.join(thread))
        jobs_follow_string = ';'.join(thread_list)
        # make sure to catch all relevant ini into needed in the newly created ini
        if cfg[FOLLOW_JOBS_LIST] is not None:
            tpl_str += '\n{} = {}'.format(FOLLOW_JOBS_LIST, jobs_follow_string)
            for key_word in KEYS_FOR_SPAWNING_SBATCH:
                # if the value is not what is default, check if it will be printed in the created ini
                # (by checking the tpl_str). If not, add to the tpl_str.
                if cfg[key_word] != DEF_CFG_VALS[key_word]:
                    if not (key_word in tpl_str):
                        # might as well just fill in the value here; easier than setting it up to be filled later
                        tpl_str += '\n{} = {}'.format(key_word, cfg[key_word])
    # Add the following keywords only if they are not the default values
    for key_word in KEYS_FOR_SPAWNING_INIS:
        if key_word not in tpl_str:
            if cfg[key_word] != DEF_CFG_VALS[key_word]:
                if key_word == OLD_CHECK_ECHO:
                    tpl_str += '\n{} = {}'.format(key_word, cfg[key_word].replace('%', '%%'))
                else:
                    tpl_str += '\n{} = {}'.format(key_word, cfg[key_word])

    # job tpl locations will be needed
    tpl_loc_list_tuples = sorted(tpl_dict.items(), key=lambda x: x[1])
    for job, tpl_loc in tpl_loc_list_tuples:
        if job != '' and job in tpl_str:
            tpl_str += '\n{} = {}'.format(job, tpl_loc)

    str_to_file(tpl_str, new_ini_fname, print_info=True)


def setup_and_submit(cfg, current_job_list, tpl_dict, testing_mode, chk_warn):
    if len(current_job_list) == 1 and current_job_list[0] == '':
        suffix = ''
    else:
        if current_job_list[0] == '':
            suffix = '_' + '_'.join(current_job_list[1:])
        else:
            suffix = '_' + '_'.join(current_job_list)
    tpl_dict[JOB_DESCRIP] = tpl_dict[JOB_NAME] + suffix

    new_ini_fname = create_out_fname(tpl_dict[JOB_DESCRIP], ext='.ini', base_dir=cfg[OUT_DIR])
    new_sbatch_fname = create_out_fname(tpl_dict[JOB_DESCRIP], ext='.slurm', base_dir=cfg[OUT_DIR])

    sbatch_dict = create_sbatch_dict(cfg, tpl_dict, os.path.relpath(new_ini_fname), current_job_list,
                                     start_from_job_name_chk=cfg[START_FROM_SAME_CHK], ignore_chk_warning=chk_warn)
    tpl_str = read_tpl(cfg[SBATCH_TPL])
    fill_save_tpl(tpl_str, sbatch_dict, cfg[SBATCH_TPL], new_sbatch_fname)

    # read ini_tpl and check if it has fields for submitting spawned jobs, if needed
    create_ini_with_req_keys(current_job_list, cfg[TPL_DICT], cfg, new_ini_fname)

    if not cfg[NO_SUBMIT]:
        # Do not want to actually (attempt to) submit a job during testing; this way, do not have to specify both
        #   testing mode and NO_SUBMIT (could make NO_SUBMIT if in testing mode, but no real advantage to that
        if testing_mode:
            sbatch_result = subprocess.check_output(["echo", "Running in testing mode: "
                                                             "'sbatch' not called"]).decode("utf-8").strip()
        else:
            #  Will not be covered in testing mode, as is not part of written code to be tested
            sbatch_result = subprocess.check_output(["sbatch", new_sbatch_fname]).decode("utf-8").strip()
        print(sbatch_result)


def main(argv=None):
    print(f"Running GaussianWrangler script run_gauss version {__version__}")
    # Read input
    args, ret = parse_cmdline(argv)
    if ret != GOOD_RET or args is None:
        return ret

    cfg = args.config

    try:
        # overwrite default values from reading config if they were specified in command line
        args_key_map = [(args.list_of_jobs, False, LIST_OF_JOBS),
                        (args.old_chk_fname, None, FIRST_JOB_CHK),
                        (args.setup_submit, False, SETUP_SUBMIT)]
        for arg_val, arg_default, cfg_key in args_key_map:
            if arg_val != arg_default:
                cfg[cfg_key] = arg_val
        # The following do not have default config options, so overwrite
        cfg[NO_SUBMIT] = args.no_submit
        if cfg[FIRST_JOB_CHK]:
            # remove extension (if any) from cfg[FIRST_JOB_CHK]
            cfg[FIRST_JOB_CHK] = os.path.splitext(cfg[FIRST_JOB_CHK])[0]

        # for the "list_of_jobs" option, "job_name" is actually the name of the name of file with the list of jobs
        if args.list_of_jobs:
            with open(args.job_name) as f:
                for line in f:
                    s_line = line.strip()
                    if len(s_line) == 0:
                        continue
                    input_job_file = os.path.splitext(s_line)[0] + cfg[GAUSS_IN_EXT]
                    base_name = get_fname_root(s_line)
                    tpl_dict = {JOB_NAME: base_name, INPUT_FILE: input_job_file}
                    for thread_index, thread in enumerate(cfg[JOB_LIST]):
                        setup_and_submit(cfg, thread, tpl_dict, args.testing, args.ignore_chk_warning)
            return GOOD_RET

        # otherwise, job_name is actually the job name. We can to ignore any extension on it
        job_name_perhaps_with_dir = os.path.splitext(args.job_name)[0]
        job_name = os.path.basename(job_name_perhaps_with_dir)
        tpl_dict = {JOB_NAME: job_name, INPUT_FILE: job_name_perhaps_with_dir + cfg[GAUSS_IN_EXT]}
        if cfg[JOB_LIST][0] == '' and not os.path.isfile(tpl_dict[INPUT_FILE]):
            raise IOError("Could not find input file: {}".format(tpl_dict[INPUT_FILE]))

        if args.setup_submit:
            for thread_index, thread in enumerate(cfg[JOB_LIST]):
                setup_and_submit(cfg, thread, tpl_dict, args.testing, args.ignore_chk_warning)
            return GOOD_RET

        for job in cfg[JOB_LIST]:
            run_job(job, job_name_perhaps_with_dir, tpl_dict, cfg, args.testing)

        if len(cfg[FOLLOW_JOBS_LIST]) > 1:
            for thread_index, thread in enumerate(cfg[FOLLOW_JOBS_LIST]):
                if thread_index == 0 and not cfg[ALL_NEW]:
                    continue
                setup_and_submit(cfg, thread, tpl_dict, args.testing, args.ignore_chk_warning)

        if len(cfg[FOLLOW_JOBS_LIST]) > 0 and not cfg[ALL_NEW]:
            for job in cfg[FOLLOW_JOBS_LIST][0]:
                run_job(job, job_name_perhaps_with_dir, tpl_dict, cfg, args.testing)

    except IOError as e:
        warning("Problems reading file:", e)
        return IO_ERROR
    except (subprocess.CalledProcessError, KeyError) as e:
        warning("", e)
    except InvalidInputError as e:
        warning("Check input:", e)
        return INVALID_DATA
    except InvalidDataError as e:
        warning("Invalid data:", e)
        return INVALID_DATA

    return GOOD_RET  # success


if __name__ == '__main__':
    status = main()
    sys.exit(status)
