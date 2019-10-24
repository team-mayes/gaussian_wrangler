#!/usr/bin/env python
"""
Script to add environment-specific information to a Gaussian job
Set up for slurm systems only

Assumptions:
   - Running on one node, and the whole node
"""
import argparse
import subprocess
import sys
from common_wrangler.common import (GOOD_RET, INPUT_ERROR, warning, InvalidDataError)


def get_node_info(node_name):
    print(node_name)
    num_procs = subprocess.check_output(["grep", "-c", "^processor", "/proc/meminfo"])
    # 97311328 kB
    # r1i6n14
    # cat  /proc/meminfo
    # num_processors = grep -c ^processor /proc/cpuinfo
    # sinfo -n r1i6n14 --format="%e %m"
    max_mem = 1000
    free_mem = 900
    alloc_mem = min(max_mem * .75, free_mem * .85)
    return num_procs, max_mem, free_mem, alloc_mem


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
    #
    # cfg = args.config

    num_procs, max_mem, free_mem, alloc_mem = get_node_info(args.node_name)
    print("NumProcs: {}".format(num_procs))
    print("MemTotal: {:14} kB\nMemFree:  {:14} kB\n"
          "MemToGauss: {:12} kB".format(num_procs, max_mem, free_mem, alloc_mem))

    # except IOError as e:
    #     warning("Problems reading file:", e)
    #     return IO_ERROR
    # except subprocess.CalledProcessError as e:
    #     warning("", e)
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
