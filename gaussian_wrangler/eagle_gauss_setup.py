#!/usr/bin/env python
"""
Script to add environment-specific information to a Gaussian job
Set up for slurm systems only

Assumptions:
   - Running on one node, and the whole node
"""
import argparse
import os
import re
import subprocess
import sys
from common_wrangler.common import (GOOD_RET, INPUT_ERROR, warning, InvalidDataError, IO_ERROR)

SOURCE_DIR = os.path.dirname(__file__)
DATA_DIR = os.path.join(SOURCE_DIR, 'data')

MEMINFO_FILE = os.path.join(DATA_DIR, 'meminfo')
CPUINFO_FILE = os.path.join(DATA_DIR, 'cpuinfo')

MEM_TOT_PAT = re.compile(r"MemTotal.*")
MEM_FREE_PAT = re.compile(r"MemFree.*")


def get_node_info(node_name, testing_mode):
    if testing_mode:
        mem_info_file = MEMINFO_FILE
        cpu_info_file = CPUINFO_FILE
    else:
        mem_info_file = "/proc/meminfo"
        cpu_info_file = "/proc/cpuinfo"

    # To make IDE happy
    mem_tot, mem_free = 0, 0

    # avoiding multiple subprocess calls to get detailed processor info, and being extra careful to check obtained data
    raw_proc_info = subprocess.check_output(["grep", "^processor", cpu_info_file]).decode("utf-8").strip().split('\n')
    print(raw_proc_info)
    num_procs = len(raw_proc_info)
    first_proc = int(raw_proc_info[0].split()[-1])
    last_proc = int(raw_proc_info[-1].split()[-1])
    assert (last_proc - first_proc + 1) == num_procs
    proc_list = "{}-{}".format(first_proc, last_proc)

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
                # no other info needed, so stop iteration
                break
    mem_alloc = int(min(mem_tot * .75, mem_free * .85))
    print('On node {}, Found {} processors, {} kB total memory and {} kB free memory.\n'
          'Will instruct Gaussian to use up to {} processors and {} kB of memory.'.format(node_name, num_procs,
                                                                                          mem_tot, mem_free,
                                                                                          num_procs, mem_alloc))
    return proc_list, mem_alloc


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
    parser.add_argument("node_name", help="The node name where the job is being run. This script assumes that exactly "
                                          "one node is used, with all processors "
                                          "needed. Otherwise, a checkpoint file (with extension '.chk') is "
                                          "needed.")
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
    parser.add_argument("-o", "--old_chk_file", help="The base name of the checkpoint file (do not include '.chk')"
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
    except (KeyError, InvalidDataError, SystemExit) as e:
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

    try:
        num_procs, mem_alloc = get_node_info(args.node_name, args.testing)
        print(num_procs, mem_alloc)

    except IOError as e:
        warning("Problems reading file:", e)
        return IO_ERROR
    except (subprocess.CalledProcessError, ValueError) as e:
        warning("", e)
    # except InvalidInputError as e:
    #     warning("Check input:", e)
    #     return INVALID_DATA
    # except InvalidDataError as e:
    #     warning("Problems reading data:", e)
    #     return INVALID_DATA

    return GOOD_RET  # success


if __name__ == '__main__':
    status = main()
    sys.exit(status)
