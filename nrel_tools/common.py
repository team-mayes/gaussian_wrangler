# coding=utf-8

"""
Common methods for this project.
"""

from __future__ import print_function, division
# Util Methods #
import argparse
import collections
import csv
import difflib
import glob
from datetime import datetime
import re
import shutil
import errno
import fnmatch
from itertools import chain, islice
import math
import numpy as np
import os
from shutil import copy2, Error, copystat
import six
import sys
from contextlib import contextmanager
# from matplotlib import rc
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
from matplotlib.patches import Rectangle

# Constants #

TPL_IO_ERR_MSG = "Couldn't read template at: '{}'"
MISSING_SEC_HEADER_ERR_MSG = "Configuration files must start with a section header such as '[main]'. Check file: {}"
BACKUP_TS_FMT = "_%Y-%m-%d_%H-%M-%S_%f"

# Boltzmann's Constant in kcal/mol Kelvin
BOLTZ_CONST = 0.0019872041
KB = 1.38064852e-23  # [J/K]

# Planck's Constant in kcal s / mol
PLANCK_CONST = 9.53707E-14
H = 6.626070e-34  # [Js]

# Universal gas constant in kcal/mol K
RG = 0.001985877534

EHPART_TO_KCAL_MOL = 627.5094709  # [kcal/mol/(Eh/part)]

XYZ_ORIGIN = np.zeros(3)

# for figures
DEF_FIG_WIDTH = 10
DEF_FIG_HEIGHT = 6
DEF_AXIS_SIZE = 20
DEF_TICK_SIZE = 15
DEF_FIG_DIR = './figs/'

# Tolerance initially based on double standard machine precision of 5 × 10−16 for float64 (decimal64)
# found to be too stringent
TOL = 0.00000000001
# similarly, use this to round away the insignificant digits!
SIG_DECIMALS = 12

# For converting atomic number to species
ATOM_NUM_DICT = {1: 'H', 2: 'He', 3: 'Li', 4: 'Be', 5: 'B', 6: 'C', 7: 'N', 8: 'O', 9: 'F', 10: 'Ne',
                 11: 'Na', 12: 'Mg', 13: 'Al', 14: 'Si', 15: 'P', 16: 'S', 17: 'Cl', 18: 'Ar',
                 19: 'K', 20: 'Ca', 21: 'Sc', 22: 'Ti', 23: 'V', 24: 'Cr', 25: 'Mn', 26: 'Fe', 27: 'Co', 28: 'Ni',
                 29: 'Cu',
                 30: 'Zn', 31: 'Ga', 32: 'Ge', 33: 'As', 34: 'Se', 35: 'Br', 36: 'Kr',
                 }

# Sections for reading files
SEC_TIMESTEP = 'timestep'
SEC_NUM_ATOMS = 'dump_num_atoms'
SEC_BOX_SIZE = 'dump_box_size'
SEC_ATOMS = 'atoms_section'
SEC_HEAD = 'head_section'
SEC_TAIL = 'tail_section'
ATOM_TYPE = 'atom_type'
ATOM_COORDS = 'atom_coords'
GAU_HEADER_PAT = re.compile(r"#.*")
# noinspection RegExpRepeatedSpace
GAU_COORD_PAT = re.compile(r"Center     Atomic      Atomic             Coordinates.*")
GAU_SEP_PAT = re.compile(r"---------------------------------------------------------------------.*")
GAU_E_PAT = re.compile(r"SCF Done:.*")
GAU_CHARGE_PAT = re.compile(r"Charge =.*")
CHARGE = 'Charge'
MULT = 'Mult'
MULTIPLICITY = 'Multiplicity'

# From template files
BASE_NAME = 'base_name'
NUM_ATOMS = 'num_atoms'
HEAD_CONTENT = 'head_content'
ATOMS_CONTENT = 'atoms_content'
TAIL_CONTENT = 'tail_content'

# Lammps-specific sections
MASSES = 'Masses'
PAIR_COEFFS = 'Pair Coeffs'
ATOMS = 'Atoms'
BOND_COEFFS = 'Bond Coeffs'
BONDS = 'Bonds'
ANGLE_COEFFS = 'Angle Coeffs'
ANGLES = 'Angles'
DIHE_COEFFS = 'Dihedral Coeffs'
DIHES = 'Dihedrals'
IMPR_COEFFS = 'Improper Coeffs'
IMPRS = 'Impropers'
LAMMPS_SECTION_NAMES = [MASSES, PAIR_COEFFS, ATOMS, BOND_COEFFS, BONDS, ANGLE_COEFFS, ANGLES,
                        DIHE_COEFFS, DIHES, IMPR_COEFFS, IMPRS]

# PDB file info
PDB_FORMAT = '{:s}{:s}{:s}{:s}{:4d}    {:8.3f}{:8.3f}{:8.3f}{:s}'
PDB_LINE_TYPE_LAST_CHAR = 6
PDB_ATOM_NUM_LAST_CHAR = 11
PDB_ATOM_TYPE_LAST_CHAR = 17
PDB_RES_TYPE_LAST_CHAR = 22
PDB_MOL_NUM_LAST_CHAR = 28
PDB_X_LAST_CHAR = 38
PDB_Y_LAST_CHAR = 46
PDB_Z_LAST_CHAR = 54
PDB_BEFORE_ELE_LAST_CHAR = 76
PDB_ELE_LAST_CHAR = 78

# Error Codes
# The good status code
GOOD_RET = 0
INPUT_ERROR = 1
IO_ERROR = 2
INVALID_DATA = 3

PY2 = sys.version_info[0] == 2
# PY3 = sys.version_info[0] == 3


# Exceptions #

class MdError(Exception):
    pass


class InvalidInputError(MdError):
    pass


class InvalidDataError(MdError):
    pass


class NotFoundError(MdError):
    pass


class ArgumentParserError(Exception):
    pass


class TemplateNotReadableError(Exception):
    pass


class ThrowingArgumentParser(argparse.ArgumentParser):
    def error(self, message):
        raise ArgumentParserError(message)


def warning(*objs):
    """Writes a message to stderr."""
    print("WARNING: ", *objs, file=sys.stderr)


# Test utilities


# From http://schinckel.net/2013/04/15/capture-and-test-sys.stdout-sys.stderr-in-unittest.testcase/
@contextmanager
def capture_stdout(command, *args, **kwargs):
    # pycharm doesn't know six very well, so ignore the false warning
    # noinspection PyCallingNonCallable
    out, sys.stdout = sys.stdout, six.StringIO()
    command(*args, **kwargs)
    sys.stdout.seek(0)
    yield sys.stdout.read()
    sys.stdout = out


@contextmanager
def capture_stderr(command, *args, **kwargs):
    # pycharm doesn't know six very well, so ignore the false warning
    # noinspection PyCallingNonCallable
    err, sys.stderr = sys.stderr, six.StringIO()
    command(*args, **kwargs)
    sys.stderr.seek(0)
    yield sys.stderr.read()
    sys.stderr = err


# Calculations #


def calc_kbt(temp_k):
    """
    Returns the given temperature in Kelvin multiplied by Boltzmann's Constant.

    @param temp_k: A temperature in Kelvin.
    @return: The given temperature in Kelvin multiplied by Boltzmann's Constant.
    """
    return BOLTZ_CONST * temp_k


def calc_k(temp, delta_gibbs):
    """
    Returns the rate coefficient calculated from Transition State Theory in inverse seconds
    @param temp: the temperature in Kelvin
    @param delta_gibbs: the change in Gibbs free energy in kcal/mol
    @return: rate coefficient in inverse seconds
    """
    return BOLTZ_CONST * temp / PLANCK_CONST * math.exp(-delta_gibbs / (RG * temp))


def pbc_calc_vector(a, b, box):
    """
    Finds the vectors between two points
    @param a: xyz coords 1
    @param b: xyz coords 2
    @param box: vector with PBC box dimensions
    @return: returns the vector a - b
    """
    vec = np.subtract(a, b)
    return vec - np.multiply(box, np.asarray(list(map(round, vec / box))))


def first_pbc_image(xyz_coords, box):
    """
    Moves xyz coords to the first PBC image, centered at the origin
    @param xyz_coords: coordinates to center (move to first image)
    @param box: PBC box dimensions
    @return: xyz coords (np array) moved to the first image
    """
    return pbc_calc_vector(xyz_coords, XYZ_ORIGIN, box)


def pbc_vector_avg(a, b, box):
    diff = pbc_calc_vector(a, b, box)
    mid_pt = np.add(b, np.divide(diff, 2.0))
    # mid-point may not be in the first periodic image. Make it so by getting its difference from the origin
    return pbc_calc_vector(mid_pt, np.zeros(len(mid_pt)), box)


def unit_vector(vector):
    """ Returns the unit vector of the vector.
    http://stackoverflow.com/questions/2827393/angles-between-two-n-dimensional-vectors-in-python
    """
    return vector / np.linalg.norm(vector)


def vec_angle(vec_1, vec_2):
    """
    Calculates the angle between the vectors (p2 - p1) and (p0 - p1)
    Note: assumes the vector calculation accounted for the PBC
    @param vec_1: xyz coordinate for the first pt
    @param vec_2: xyz for 2nd pt
    @return: the angle in between the vectors
    """
    unit_vec_1 = unit_vector(vec_1)
    unit_vec_2 = unit_vector(vec_2)

    return np.rad2deg(np.arccos(np.clip(np.dot(unit_vec_1, unit_vec_2), -1.0, 1.0)))


def vec_dihedral(vec_ba, vec_bc, vec_cd):
    """
    calculates the dihedral angle from the vectors b --> a, b --> c, c --> d
    where a, b, c, and d are the four points
    From:
    http://stackoverflow.com/questions/20305272/
      dihedral-torsion-angle-from-four-points-in-cartesian-coordinates-in-python
    Khouli formula
    1 sqrt, 1 cross product
    @param vec_ba: the vector connecting points b --> a, accounting for pbc
    @param vec_bc: b --> c
    @param vec_cd: c --> d
    @return: dihedral angle in degrees
    """
    # normalize b1 so that it does not influence magnitude of vector
    # rejections that come next
    vec_bc = unit_vector(vec_bc)

    # vector rejections
    # v = projection of b0 onto plane perpendicular to b1
    #   = b0 minus component that aligns with b1
    # w = projection of b2 onto plane perpendicular to b1
    #   = b2 minus component that aligns with b1
    v = vec_ba - np.dot(vec_ba, vec_bc) * vec_bc
    w = vec_cd - np.dot(vec_cd, vec_bc) * vec_bc

    # angle between v and w in a plane is the torsion angle
    # v and w may not be normalized but that's fine since tan is y/x
    x = np.dot(v, w)
    y = np.dot(np.cross(vec_bc, v), w)
    return np.degrees(np.arctan2(y, x))


# Other #

def chunk(seq, chunk_size, process=iter):
    """Yields items from an iterator in iterable chunks.
    From https://gist.github.com/ksamuel/1275417

    @param seq: The sequence to chunk.
    @param chunk_size: The size of the returned chunks.
    @param process: The function to use for creating the iterator.  Useful for iterating over different
    data structures.
    @return: Chunks of the given size from the given sequence.
    """
    it = iter(seq)
    while True:
        yield process(chain([six.next(it)], islice(it, chunk_size - 1)))


# I/O #

def read_tpl(tpl_loc):
    """Attempts to read the given template location and throws A
    TemplateNotReadableError if it can't read the given location.

    :param tpl_loc: The template location to read.
    :raise TemplateNotReadableError: If there is an IOError reading the location.
    """
    try:
        return file_to_str(tpl_loc)
    except IOError:
        raise TemplateNotReadableError(TPL_IO_ERR_MSG.format(tpl_loc))


def make_dir(tgt_dir):
    """
    Creates the given directory and its parent directories if it
    does not already exist.

    Keyword arguments:
    tgt_dir -- The directory to create

    """
    if not os.path.exists(tgt_dir):
        os.makedirs(tgt_dir)
    elif not os.path.isdir(tgt_dir):
        raise NotFoundError("Resource {} exists and is not a dir".format(tgt_dir))


def file_to_str(f_name):
    """
    Reads and returns the contents of the given file.

    @param f_name: The location of the file to read.
    @return: The contents of the given file.
    :raises: IOError if the file can't be opened for reading.
    """
    with open(f_name) as f:
        return f.read()


def file_rows_to_list(c_file):
    """
    Given the name of a file, returns a list of its rows, after filtering out empty rows
    @param c_file: file location
    @return: list of non-empty rows
    """
    with open(c_file) as f:
        row_list = [row.strip() for row in f.readlines()]
        return list(filter(None, row_list))


def str_to_file(str_val, f_name, mode='w', print_info=False):
    """
    Writes the string to the given file.
    @param str_val: The string to write.
    @param f_name: The location of the file to write
    @param mode: default mode is to overwrite file
    @param print_info: boolean to specify whether to print action to stdout
    """
    with open(f_name, mode) as f:
        f.write(str_val)
    if print_info:
        print("Wrote file: {}".format(f_name))


def round_to_print(val):
    """
    To remove floating point digits that are imprecise due to machine precision
    @param val: a float
    @return: a float without insignificant digits
    """
    return round(val, SIG_DECIMALS)


def np_float_array_from_file(data_file, delimiter=" ", header=False, gather_hist=False):
    """
    Adds to the basic np.loadtxt by performing data checks.
    @param data_file: file expected to have space-separated values, with the same number of entries per row
    @param delimiter: default is a space-separated file
    @param header: default is no header; alternately, specify number of header lines
    @param gather_hist: default is false; gather data to make histogram of non-numerical data
    @return: a numpy array or InvalidDataError if unsuccessful, followed by the header_row (None if none specified)
    """
    header_row = None
    hist_data = {}
    with open(data_file) as csv_file:
        csv_list = list(csv.reader(csv_file, delimiter=delimiter))
    if header:
        header_row = csv_list[0]

    try:
        data_array = np.genfromtxt(data_file, dtype=np.float64, delimiter=delimiter, skip_header=header)
    except ValueError:
        data_array = None
        line_len = None
        if header:
            first_line = 1
        else:
            first_line = 0
        for row in csv_list[first_line:]:
            if len(row) == 0:
                continue
            s_len = len(row)
            if line_len is None:
                line_len = s_len
            elif s_len != line_len:
                raise InvalidDataError('File could not be read as an array of floats: {}\n  Expected '
                                       'values separated by "{}" with an equal number of columns per row.\n'
                                       '  However, found {} values on the first data row'
                                       '  and {} values on the later row: "{}")'
                                       ''.format(data_file, delimiter, line_len, s_len, row))
            data_vector = np.empty([line_len], dtype=np.float64)
            for col in range(line_len):
                try:
                    data_vector[col] = float(row[col])
                except ValueError:
                    data_vector[col] = np.nan
                    if gather_hist:
                        col_key = str(row[col])
                        if col in hist_data:
                            if col_key in hist_data[col]:
                                hist_data[col][col_key] += 1
                            else:
                                hist_data[col][col_key] = 1
                        else:
                            hist_data[col] = {col_key: 1}
            if data_array is None:
                data_array = np.copy(data_vector)
            else:
                data_array = np.vstack((data_array, data_vector))
    if len(data_array.shape) == 1:
        raise InvalidDataError("File contains a vector, not an array of floats: {}\n".format(data_file))
    if np.isnan(data_array).any():
        warning("Encountered entry (or entries) which could not be converted to a float. "
                "'nan' will be returned for the stats for that column.")
    return data_array, header_row, hist_data


def read_csv_to_list(data_file, delimiter=',', header=False):
    """
    Reads file of values; did not use np.loadtxt because can have floats and strings
    @param data_file: name of delimiter-separated file with the same number of entries per row
    @param delimiter: string: delimiter between column values
    @param header: boolean to denote if file contains a header
    @return: a list containing the data (removing header row, if one is specified) and a list containing the
             header row (empty if no header row specified)
    """
    with open(data_file) as csv_file:
        csv_list = list(csv.reader(csv_file, delimiter=delimiter, quoting=csv.QUOTE_NONNUMERIC))

    header_row = []

    if header:
        first_line = 1
        header_row = csv_list[0]
    else:
        first_line = 0

    return csv_list[first_line:], header_row


def create_backup_filename(orig):
    base, ext = os.path.splitext(orig)
    now = datetime.now()
    return "".join((base, now.strftime(BACKUP_TS_FMT), ext))


def find_backup_filenames(orig):
    base, ext = os.path.splitext(orig)
    found = glob.glob(base + "*" + ext)
    try:
        found.remove(orig)
    except ValueError:
        # Original not present; ignore.
        pass
    return found


def silent_remove(filename, disable=False):
    """
    Removes the target file name, catching and ignoring errors that indicate that the
    file does not exist.

    @param filename: The file to remove.
    @param disable: boolean to flag if want to disable removal
    """
    if not disable:
        try:
            if os.path.isdir(filename):
                os.rmdir(filename)
            else:
                os.remove(filename)
        except OSError as e:
            if e.errno != errno.ENOENT:
                raise


def allow_write(f_loc, overwrite=False):
    """
    Returns whether to allow writing to the given location.

    @param f_loc: The location to check.
    @param overwrite: Whether to allow overwriting an existing location.
    @return: Whether to allow writing to the given location.
    """
    if os.path.exists(f_loc) and not overwrite:
        warning("Not overwriting existing file '{}'".format(f_loc))
        return False
    return True


def move_existing_file(f_loc):
    """
    Renames an existing file using a timestamp based on the move time.

    @param f_loc: The location to check.
    """
    if os.path.exists(f_loc):
        shutil.move(f_loc, create_backup_filename(f_loc))


def get_fname_root(src_file):
    """

    @param src_file:
    @return: the file root name (no directory, no extension)
    """
    return os.path.splitext(os.path.basename(src_file))[0]


def create_out_fname(src_file, prefix='', suffix='', remove_prefix=None, base_dir=None, ext=None):
    """Creates an outfile name for the given source file.

    @param remove_prefix: string to remove at the beginning of file name
    @param src_file: The file to process.
    @param prefix: The file prefix to add, if specified.
    @param suffix: The file suffix to append, if specified.
    @param base_dir: The base directory to use; defaults to `src_file`'s directory.
    @param ext: The extension to use instead of the source file's extension;
        defaults to the `scr_file`'s extension.
    @return: The output file name.
    """

    if base_dir is None:
        base_dir = os.path.dirname(src_file)

    file_name = os.path.basename(src_file)
    if remove_prefix is not None and file_name.startswith(remove_prefix):
        base_name = file_name[len(remove_prefix):]
    else:
        base_name = os.path.splitext(file_name)[0]

    if ext is None:
        ext = os.path.splitext(file_name)[1]

    return os.path.abspath(os.path.join(base_dir, prefix + base_name + suffix + ext))


def find_files_by_dir(tgt_dir, pat):
    """Recursively searches the target directory tree for files matching the given pattern.
    The results are returned as a dict with a list of found files keyed by the absolute
    directory name.
    @param tgt_dir: The target base directory.
    @param pat: The file pattern to search for.
    @return: A dict where absolute directory names are keys for lists of found file names
        that match the given pattern.
    """
    match_dirs = {}
    for root, dirs, files in os.walk(tgt_dir):
        matches = [match for match in files if fnmatch.fnmatch(match, pat)]
        if matches:
            match_dirs[os.path.abspath(root)] = matches
    return match_dirs


def copytree(src, dst, symlinks=False, ignore=None):
    """This is a copy of the standard Python shutil.copytree, but it
    allows for an existing destination directory.

    Recursively copy a directory tree using copy2().

    If exception(s) occur, an Error is raised with a list of reasons.

    If the optional symlinks flag is true, symbolic links in the
    source tree result in symbolic links in the destination tree; if
    it is false, the contents of the files pointed to by symbolic
    links are copied.

    The optional ignore argument is a callable. If given, it
    is called with the `src` parameter, which is the directory
    being visited by copytree(), and `names` which is the list of
    `src` contents, as returned by os.listdir():

        callable(src, names) -> ignored_names

    Since copytree() is called recursively, the callable will be
    called once for each directory that is copied. It returns a
    list of names relative to the `src` directory that should
    not be copied.

    XXX Consider this example code rather than the ultimate tool.

    @param src: The source directory.
    @param dst: The destination directory.
    @param symlinks: Whether to follow symbolic links.
    @param ignore: A callable for items to ignore at a given level.
    """
    names = os.listdir(src)
    if ignore is not None:
        ignored_names = ignore(src, names)
    else:
        ignored_names = set()

    if not os.path.exists(dst):
        os.makedirs(dst)

    errors = []
    for name in names:
        if name in ignored_names:
            continue
        src_name = os.path.join(src, name)
        dst_name = os.path.join(dst, name)
        try:
            if symlinks and os.path.islink(src_name):
                link_to = os.readlink(src_name)
                os.symlink(link_to, dst_name)
            elif os.path.isdir(src_name):
                copytree(src_name, dst_name, symlinks, ignore)
            else:
                # Will raise a SpecialFileError for unsupported file types
                copy2(src_name, dst_name)
        # catch the Error from the recursive copytree so that we can
        # continue with other files
        except Error as err:
            errors.extend(err.args[0])
        except EnvironmentError as why:
            errors.append((src_name, dst_name, str(why)))
    try:
        copystat(src, dst)
    except OSError as why:
        # can't copy file access times on Windows
        # noinspection PyUnresolvedReferences
        if why.winerror is None:
            errors.extend((src, dst, str(why)))
    if errors:
        raise Error(errors)


# CSV #

def read_csv_header(src_file):
    """Returns a list containing the values from the first row of the given CSV
    file or None if the file is empty.

    @param src_file: The CSV file to read.
    @return: The first row or None if empty.
    """
    with open(src_file) as csv_file:
        for row in csv.reader(csv_file):
            return list(row)


def convert_dict_line(all_conv, data_conv, line):
    s_dict = {}
    for s_key, s_val in line.items():
        if data_conv and s_key in data_conv:
            try:
                s_dict[s_key] = data_conv[s_key](s_val)
            except ValueError as e:
                warning("Could not convert value '{}' from column '{}': '{}'.  Leaving as str".format(s_val, s_key, e))
                s_dict[s_key] = s_val
        elif all_conv:
            try:
                s_dict[s_key] = all_conv(s_val)
            except ValueError as e:
                warning("Could not convert value '{}' from column '{}': '{}'.  Leaving as str".format(s_val, s_key, e))
                s_dict[s_key] = s_val
        else:
            s_dict[s_key] = s_val
    return s_dict


def read_csv(src_file, data_conv=None, all_conv=None, quote_style=csv.QUOTE_MINIMAL):
    """
    Reads the given CSV (comma-separated with a first-line header row) and returns a list of
    dicts where each dict contains a row's data keyed by the header row.

    @param src_file: The CSV to read.
    @param data_conv: A map of header keys to conversion functions.  Note that values
        that throw a TypeError from an attempted conversion are left as strings in the result.
    @param all_conv: A function to apply to all values in the CSV.  A specified data_conv value
        takes precedence.
    @param quote_style: how to read the dictionary
    @return: A list of dicts containing the file's data.
    """
    result = []
    with open(src_file) as csv_file:
        csv_reader = csv.DictReader(csv_file, quoting=quote_style)
        for line in csv_reader:
            result.append(convert_dict_line(all_conv, data_conv, line))
    return result


def read_csv_to_dict(src_file, col_name, data_conv=None, all_conv=None):
    """
    Reads the given CSV (comma-separated with a first-line header row) and returns a
    dict of dicts indexed on the given col_name. Each dict contains a row's data keyed by the header row.

    @param src_file: The CSV to read.
    @param col_name: the name of the column to index on
    @param data_conv: A map of header keys to conversion functions.  Note that values
        that throw a TypeError from an attempted conversion are left as strings in the result.
    @param all_conv: A function to apply to all values in the CSV.  A specified data_conv value
        takes precedence.
    @return: A list of dicts containing the file's data.
    """
    result = {}
    with open(src_file) as csv_file:
        try:
            csv_reader = csv.DictReader(csv_file, quoting=csv.QUOTE_NONNUMERIC)
            create_dict(all_conv, col_name, csv_reader, data_conv, result, src_file)
        except ValueError:
            csv_reader = csv.DictReader(csv_file)
            create_dict(all_conv, col_name, csv_reader, data_conv, result, src_file)
    return result


def create_dict(all_conv, col_name, csv_reader, data_conv, result, src_file):
    for line in csv_reader:
        val = convert_dict_line(all_conv, data_conv, line)
        if col_name in val:
            try:
                col_val = int(val[col_name])
            except ValueError:
                col_val = val[col_name]
            if col_val in result:
                warning("Duplicate values found for {}. Value for key will be overwritten.".format(col_val))
            result[col_val] = convert_dict_line(all_conv, data_conv, line)
        else:
            raise InvalidDataError("Could not find value for {} in file {} on line {}."
                                   "".format(col_name, src_file, line))


def write_csv(data, out_fname, fieldnames, extrasaction="raise", mode='w', quote_style=csv.QUOTE_NONNUMERIC,
              print_message=True, round_digits=False):
    """
    Writes the given data to the given file location.

    @param round_digits: if desired, provide decimal number for rounding
    @param data: The data to write (list of dicts).
    @param out_fname: The name of the file to write to.
    @param fieldnames: The sequence of field names to use for the header.
    @param extrasaction: What to do when there are extra keys.  Acceptable
        values are "raise" or "ignore".
    @param mode: default mode is to overwrite file
    @param print_message: boolean to flag whether to note that file written or appended
    @param quote_style: dictates csv output style
    """
    with open(out_fname, mode) as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames, extrasaction=extrasaction, quoting=quote_style)
        if mode == 'w':
            writer.writeheader()
        if round_digits:
            for row_id in range(len(data)):
                new_dict = {}
                for key, val in data[row_id].items():
                    if isinstance(val, float):
                        new_dict[key] = round(val, round_digits)
                    else:
                        new_dict[key] = val
                data[row_id] = new_dict
        writer.writerows(data)
    if print_message:
        if mode == 'a':
            print("  Appended: {}".format(out_fname))
        elif mode == 'w':
            print("Wrote file: {}".format(out_fname))


def list_to_csv(data, out_fname, delimiter=',', mode='w', quote_style=csv.QUOTE_NONNUMERIC,
                print_message=True, round_digits=False):
    """
    Writes the given data to the given file location.
    @param data: The data to write (list of lists).
    @param out_fname: The name of the file to write to.
    @param delimiter: string
    @param mode: default mode is to overwrite file
    @param quote_style: csv quoting style
    @param print_message: boolean to allow update
    @param round_digits: boolean to affect printing output; supply an integer to round to that number of decimals
    """
    with open(out_fname, mode) as csv_file:
        writer = csv.writer(csv_file, delimiter=delimiter, quoting=quote_style)
        if round_digits:
            for row_id in range(len(data)):
                new_row = []
                for val in data[row_id]:
                    if isinstance(val, float):
                        new_row.append(round(val, round_digits))
                    else:
                        new_row.append(val)
                data[row_id] = new_row
        writer.writerows(data)
    if print_message:
        print("Wrote file: {}".format(out_fname))


# Other input/output files


def read_csv_dict(d_file, ints=True, one_to_one=True, pdb_dict=False, str_float=False):
    """
    If an dictionary file is given, read it and return the dict[old]=new.
    Checks that all keys are unique.
    If one_to_one=True, checks that there 1:1 mapping of keys and values.

    @param d_file: the file with csv of old_id,new_id
    @param ints: boolean to indicate if the values are to be read as integers
    @param one_to_one: flag to check for one-to-one mapping in the dict
    @param pdb_dict: flag to format as required for the PDB output
    @param str_float: indicates dictionary is a string followed by a float
    @return: new_dict
    """
    new_dict = {}
    if pdb_dict:
        ints = False
        one_to_one = False
    elif str_float:
        ints = False
        one_to_one = False
        pdb_dict = False
    # If d_file is None, return the empty dictionary, as no dictionary file was specified
    if d_file is not None:
        with open(d_file) as csv_file:
            reader = csv.reader(csv_file)
            key_count = 0
            for row in reader:
                if len(row) == 0:
                    continue
                if len(row) == 2:
                    if pdb_dict:
                        atom_type = row[0].strip()
                        type_len = len(atom_type)
                        element_type = row[1].strip()
                        if len(element_type) > 2 or type_len > 4:
                            raise InvalidDataError("Error reading line '{}' in file: {}\n  "
                                                   "Expected to read atom_type,element_type, with atom type no more "
                                                   "than 4 characters and element_type no more than 2."
                                                   "".format(row, d_file))
                        if type_len == 4:
                            atom_type = ' {:s} '.format(atom_type)
                        else:
                            atom_type = '  {:4s}'.format(atom_type)
                        new_dict[atom_type] = '{:>2s}'.format(element_type)
                    elif ints:
                        new_dict[int(row[0])] = int(row[1])
                    elif str_float:
                        new_dict[row[0]] = float(row[1])
                    else:
                        new_dict[row[0]] = row[1]
                    key_count += 1
                else:
                    raise InvalidDataError("Error reading line '{}' in file: {}\n"
                                           "  Expected exactly two comma-separated values per row."
                                           "".format(row, d_file))
        if key_count == len(new_dict):
            if one_to_one:
                for key in new_dict:
                    if not (key in new_dict.values()):
                        raise InvalidDataError('Did not find a 1:1 mapping of key,val ids in {}'.format(d_file))
        else:
            raise InvalidDataError('A non-unique key value (first column) found in file: {}\n'.format(d_file))
    return new_dict


def create_element_dict(dict_file, pdb_dict=True, one_to_one=False):
    # This is used when need to add atom types to PDB file
    element_dict = {}
    if dict_file is not None:
        return read_csv_dict(dict_file, pdb_dict=pdb_dict, ints=False, one_to_one=one_to_one)
    return element_dict


def list_to_file(list_to_print, fname, list_format=None, delimiter=' ', mode='w', print_message=True):
    """
    Writes the list of sequences to the given file in the specified format for a PDB.

    @param list_to_print: A list of lines to print. The list may be a list of lists, list of strings, or a mixture.
    @param fname: The location of the file to write.
    @param list_format: Specified formatting for the line if the line is  list.
    @param delimiter: If no format is given and the list contains lists, the delimiter will join items in the list.
    @param print_message: boolean to determine whether to write to output if the file is printed or appended
    @param mode: write by default; can be changed to allow appending to file.
    """
    with open(fname, mode) as w_file:
        for line in list_to_print:
            if isinstance(line, six.string_types):
                w_file.write(line + '\n')
            elif isinstance(line, collections.Iterable):
                if list_format is None:
                    w_file.write(delimiter.join(map(str, line)) + "\n")
                else:
                    w_file.write(list_format.format(*line) + '\n')
    if print_message:
        rel_path_fname = os.path.relpath(fname)
        if mode == 'w':
            print("Wrote file: {}".format(rel_path_fname))
        elif mode == 'a':
            print("  Appended: {}".format(rel_path_fname))


def print_qm_kind(int_list, element_name, fname, mode='w'):
    """
    Writes the list to the given file, formatted for CP2K to read as qm atom indices.

    @param int_list: The list to write.
    @param element_name: element type to designate
    @param fname: The location of the file to write.
    @param mode: default is to write to a new file. Use option to designate to append to existing file.
    """
    with open(fname, mode) as m_file:
        m_file.write('    &QM_KIND {}\n'.format(element_name))
        m_file.write('        MM_INDEX {}\n'.format(' '.join(map(str, int_list))))
        m_file.write('    &END QM_KIND\n')
    if mode == 'w':
        print("Wrote file: {}".format(fname))


def print_mm_kind(atom_type, radius, fname, mode='w'):
    """
    Writes the list to the given file, formatted for CP2K to read as qm atom indices.

    @param atom_type: (str) MM atom type
    @param radius: radius to list for covalent radius (smoothing point charge)
    @param fname: The location of the file to write.
    @param mode: default is to write to a new file. Use option to designate to append to existing file.
    """
    with open(fname, mode) as m_file:
        m_file.write('    &MM_KIND {}\n'.format(atom_type))
        m_file.write('        RADIUS {}\n'.format(radius))
        m_file.write('    &END MM_KIND\n')
    if mode == 'w':
        print("Wrote file: {}".format(fname))


def print_qm_links(c_alpha_dict, c_beta_dict, f_name, mode="w"):
    """
    Note: this needs to be tested. Only ran once to get the protein residues set up correctly.
    @param c_alpha_dict: dict of protein residue to be broken to c_alpha atom id
    @param c_beta_dict: as above, but for c_beta
    @param f_name: The location of the file to write.
    @param mode: default is to write to a new file. Use option to designate to append to existing file.
    """
    with open(f_name, mode) as m_file:
        for resid in c_beta_dict:
            m_file.write('    !! Break resid {} between CA and CB, and cap CB with hydrogen\n'
                         '    &LINK\n       MM_INDEX  {}  !! CA\n       QM_INDEX  {}  !! CB\n'
                         '       LINK_TYPE  IMOMM\n       ALPHA_IMOMM  1.5\n'
                         '    &END LINK\n'.format(resid, c_alpha_dict[resid], c_beta_dict[resid]))
    if mode == 'w':
        print("Wrote file: {}".format(f_name))


# Conversions #

def to_int_list(raw_val):
    return_vals = []
    for val in raw_val.split(','):
        return_vals.append(int(val.strip()))
    return return_vals


def to_list(raw_val):
    return_vals = []
    for val in raw_val.split(','):
        return_vals.append(val.strip())
    return return_vals


def str_to_bool(s):
    """
    Basic converter for Python boolean values written as a str.
    @param s: The value to convert.
    @return: The boolean value of the given string.
    @raises: ValueError if the string value cannot be converted.
    """
    if s == 'True':
        return True
    elif s == 'False':
        return False
    else:
        raise ValueError("Cannot covert {} to a bool".format(s))


def fmt_row_data(raw_data, fmt_str):
    """ Formats the values in the dicts in the given list of raw data using
    the given format string.

    *This may not be needed at all*
    Now that I'm using csv.QUOTE_NONNUMERIC, generally don't want to format floats to strings

    @param raw_data: The list of dicts to format.
    @param fmt_str: The format string to use when formatting.
    @return: The formatted list of dicts.
    """
    fmt_rows = []
    for row in raw_data:
        fmt_row = {}
        for key, raw_val in row.items():
            fmt_row[key] = fmt_str.format(raw_val)
        fmt_rows.append(fmt_row)
    return fmt_rows


def conv_raw_val(param, def_val, int_list=True):
    """
    Converts the given parameter into the given type (default returns the raw value).  Returns the default value
    if the param is None.
    @param param: The value to convert.
    @param def_val: The value that determines the type to target.
    @param int_list: flag to specify if lists should converted to a list of integers
    @return: The converted parameter value.
    """
    if param is None:
        return def_val
    if isinstance(def_val, bool):
        if param in ['T', 't', 'true', 'TRUE', 'True']:
            return True
        else:
            return False
    if isinstance(def_val, int):
        return int(param)
    if isinstance(def_val, float):
        return float(param)
    if isinstance(def_val, list):
        if int_list:
            return to_int_list(param)
        else:
            return to_list(param)
    return param


def process_cfg(raw_cfg, def_cfg_vals=None, req_keys=None, int_list=True, store_extra_keys=False):
    """
    Converts the given raw configuration, filling in defaults and converting the specified value (if any) to the
    default value's type.
    @param raw_cfg: The configuration map.
    @param def_cfg_vals: dictionary of default values
    @param req_keys: dictionary of required types
    @param int_list: flag to specify if lists should converted to a list of integers
    @param store_extra_keys: boolean to skip error if there are unexpected keys
    @return: The processed configuration.

    """
    proc_cfg = {}
    extra_keys = []
    for key in raw_cfg:
        if not (key in def_cfg_vals or key in req_keys):
            if store_extra_keys:
                extra_keys.append(key)
            else:
                raise InvalidDataError("Unexpected key '{}' in configuration ('ini') file.".format(key))
    key = None
    try:
        for key, def_val in def_cfg_vals.items():
            proc_cfg[key] = conv_raw_val(raw_cfg.get(key), def_val, int_list)
        for key, type_func in req_keys.items():
            proc_cfg[key] = type_func(raw_cfg[key])
        for key in extra_keys:
            proc_cfg[key] = raw_cfg[key]
    except KeyError as e:
        raise KeyError("Missing config val for key '{}'".format(key, e))
    except Exception as e:
        raise InvalidDataError('Problem with config vals on key {}: {}'.format(key, e))

    return proc_cfg


def dequote(s):
    """
    from: http://stackoverflow.com/questions/3085382/python-how-can-i-strip-first-and-last-double-quotes
    If a string has single or double quotes around it, remove them.
    Make sure the pair of quotes match.
    If a matching pair of quotes is not found, return the string unchanged.
    """
    if isinstance(s, str) and len(s) > 0:
        if (s[0] == s[-1]) and s.startswith(("'", '"')):
            return s[1:-1]
    return s


def quote(s):
    """
    Converts a variable into a quoted string
    """
    if (s[0] == s[-1]) and s.startswith(("'", '"')):
        return str(s)
    return '"' + str(s) + '"'


def single_quote(s):
    """
    Converts a variable into a quoted string
    """
    if s[0] == s[-1]:
        if s.startswith("'"):
            return str(s)
        elif s.startswith('"'):
            s = dequote(s)
    return "'" + str(s) + "'"


# Comparisons #

def conv_num(s):
    try:
        return int(s)
    except ValueError:
        try:
            return float(s)
        except ValueError:
            return s


def diff_lines(floc1, floc2, delimiter=","):
    """
    Determine all lines in a file are equal.
    This function became complicated because of edge cases:
        Do not want to flag files as different if the only difference is due to machine precision diffs of floats
    Thus, if the files are not immediately found to be the same:
        If not, test if the line is a csv that has floats and the difference is due to machine precision.
        Be careful if one value is a np.nan, but not the other (the diff evaluates to zero)
        If not, return all lines with differences.
    @param floc1: file location 1
    @param floc2: file location 1
    @param delimiter: defaults to CSV
    @return: a list of the lines with differences
    """
    diff_lines_list = []
    # Save diffs to strings to be converted to use csv parser
    output_plus = ""
    output_neg = ""
    with open(floc1, 'r') as file1:
        with open(floc2, 'r') as file2:
            diff = list(difflib.ndiff(file1.read().splitlines(), file2.read().splitlines()))

    for line in diff:
        if line.startswith('-') or line.startswith('+'):
            diff_lines_list.append(line)
            if line.startswith('-'):
                output_neg += line[2:] + '\n'
            elif line.startswith('+'):
                output_plus += line[2:] + '\n'

    if len(diff_lines_list) == 0:
        return diff_lines_list

    warning("Checking for differences between files {} {}".format(floc1, floc2))
    try:
        # take care of parentheses
        for char in ('(', ')', '[', ']'):
            output_plus = output_plus.replace(char, delimiter)
            output_neg = output_neg.replace(char, delimiter)
        # pycharm doesn't know six very well
        # noinspection PyCallingNonCallable
        diff_plus_lines = list(csv.reader(six.StringIO(output_plus), delimiter=delimiter, quoting=csv.QUOTE_NONNUMERIC))
        # noinspection PyCallingNonCallable
        diff_neg_lines = list(csv.reader(six.StringIO(output_neg), delimiter=delimiter, quoting=csv.QUOTE_NONNUMERIC))
    except ValueError:
        diff_plus_lines = output_plus.split('\n')
        diff_neg_lines = output_neg.split('\n')
        for diff_list in [diff_plus_lines, diff_neg_lines]:
            for line_id in range(len(diff_list)):
                # noinspection PyTypeChecker
                diff_list[line_id] = [x.strip() for x in diff_list[line_id].split(delimiter)]

    if len(diff_plus_lines) == len(diff_neg_lines):
        # if the same number of lines, there is a chance that the difference is only due to difference in
        # floating point precision. Check each value of the line, split on whitespace or comma
        diff_lines_list = []
        for line_plus, line_neg in zip(diff_plus_lines, diff_neg_lines):
            if len(line_plus) == len(line_neg):
                # print("Checking for differences between: ", line_neg, line_plus)
                for item_plus, item_neg in zip(line_plus, line_neg):
                    try:
                        item_plus = float(item_plus)
                        item_neg = float(item_neg)
                        # if difference greater than the tolerance, the difference is not just precision
                        # Note: if only one value is nan, the float diff is zero!
                        #  Thus, check for diffs only if neither are nan; show different if only one is nan
                        diff_vals = False
                        if np.isnan(item_neg) != np.isnan(item_plus):
                            diff_vals = True
                            warning("Comparing '{}' to '{}'.".format(item_plus, item_neg))
                        elif not (np.isnan(item_neg) and np.isnan(item_plus)):
                            # noinspection PyTypeChecker
                            if not np.isclose(item_neg, item_plus, TOL):
                                diff_vals = True
                                warning("Values {} and {} differ.".format(item_plus, item_neg))
                        if diff_vals:
                            diff_lines_list.append("- " + " ".join(map(str, line_neg)))
                            diff_lines_list.append("+ " + " ".join(map(str, line_plus)))
                            break
                    except ValueError:
                        # not floats, so the difference is not just precision
                        if item_plus != item_neg:
                            diff_lines_list.append("- " + " ".join(map(str, line_neg)))
                            diff_lines_list.append("+ " + " ".join(map(str, line_plus)))
                            break
            # Not the same number of items in the lines
            else:
                diff_lines_list.append("- " + " ".join(map(str, line_neg)))
                diff_lines_list.append("+ " + " ".join(map(str, line_plus)))
    return diff_lines_list


# Data Structures #

def unique_list(a_list):
    """ Creates an ordered list from a list of tuples or other hashable items.
    From https://code.activestate.com/recipes/576694/#c6
    """
    m_map = {}
    o_set = []
    for item in a_list:
        if item not in m_map:
            m_map[item] = 1
            o_set.append(item)
    return o_set


def conv_str_to_func(func_name):
    """
    Convert a name of a function into a function, if possible
    @param func_name: string to be converted (if possible)
    @return: either the function or error
    """
    name_func_dict = {"None": None,
                      "str": str,
                      "int": int,
                      "float": float,
                      "bool": bool,
                      }
    if func_name is None:
        return func_name
    elif func_name in name_func_dict:
        return name_func_dict[func_name]
    else:
        raise InvalidDataError("Invalid type entry '{}'. Valid options are ")


# Processing LAMMPS files #

def find_dump_section_state(line, sec_timestep=SEC_TIMESTEP, sec_num_atoms=SEC_NUM_ATOMS, sec_box_size=SEC_BOX_SIZE,
                            sec_atoms=SEC_ATOMS):
    atoms_pat = re.compile(r"^ITEM: ATOMS id mol type q x y z.*")
    if line == 'ITEM: TIMESTEP':
        return sec_timestep
    elif line == 'ITEM: NUMBER OF ATOMS':
        return sec_num_atoms
    elif line == 'ITEM: BOX BOUNDS pp pp pp':
        return sec_box_size
    elif atoms_pat.match(line):
        return sec_atoms


def process_pdb_file(pdb_file, atom_info_only=False):
    if atom_info_only:
        pdb_data = {NUM_ATOMS: 0, SEC_HEAD: [], SEC_ATOMS: {}, SEC_TAIL: []}
    else:
        pdb_data = {NUM_ATOMS: 0, SEC_HEAD: [], SEC_ATOMS: [], SEC_TAIL: []}
    atom_id = 0

    with open(pdb_file) as f:
        for line in f:
            line = line.strip()
            if len(line) == 0:
                continue
            line_head = line[:PDB_LINE_TYPE_LAST_CHAR]
            # head_content to contain Everything before 'Atoms' section
            # also capture the number of atoms
            # match 5 letters so don't need to set up regex for the ones that have numbers following the letters
            # noinspection SpellCheckingInspection
            if line_head[:-1] in ['HEADE', 'TITLE', 'REMAR', 'CRYST', 'MODEL', 'COMPN',
                                  'NUMMD', 'ORIGX', 'SCALE', 'SOURC', 'AUTHO', 'CAVEA',
                                  'EXPDT', 'MDLTY', 'KEYWD', 'OBSLT', 'SPLIT', 'SPRSD',
                                  'REVDA', 'JRNL ', 'DBREF', 'SEQRE', 'HET  ', 'HETNA',
                                  'HETSY', 'FORMU', 'HELIX', 'SHEET', 'SSBON', 'LINK ',
                                  'CISPE', 'SITE ', ]:
                # noinspection PyTypeChecker
                pdb_data[SEC_HEAD].append(line)

            # atoms_content to contain everything but the xyz
            elif line_head == 'ATOM  ' or line_head == 'HETATM':

                # By renumbering, handles the case when a PDB template has ***** after atom_id 99999.
                # For renumbering, making sure prints in the correct format, including num of characters:
                atom_id += 1
                if atom_id > 99999:
                    atom_num = format(atom_id, 'x')
                else:
                    atom_num = '{:5d}'.format(atom_id)
                # Alternately, use this:
                # atom_num = line[cfg[PDB_LINE_TYPE_LAST_CHAR]:cfg[PDB_ATOM_NUM_LAST_CHAR]]

                atom_type = line[PDB_ATOM_NUM_LAST_CHAR:PDB_ATOM_TYPE_LAST_CHAR]
                res_type = line[PDB_ATOM_TYPE_LAST_CHAR:PDB_RES_TYPE_LAST_CHAR]
                mol_num = int(line[PDB_RES_TYPE_LAST_CHAR:PDB_MOL_NUM_LAST_CHAR])
                pdb_x = float(line[PDB_MOL_NUM_LAST_CHAR:PDB_X_LAST_CHAR])
                pdb_y = float(line[PDB_X_LAST_CHAR:PDB_Y_LAST_CHAR])
                pdb_z = float(line[PDB_Y_LAST_CHAR:PDB_Z_LAST_CHAR])
                last_cols = line[PDB_Z_LAST_CHAR:]
                element_type = line[PDB_BEFORE_ELE_LAST_CHAR:PDB_ELE_LAST_CHAR]

                if atom_info_only:
                    atom_xyz = np.array([pdb_x, pdb_y, pdb_z])
                    pdb_data[SEC_ATOMS][atom_id] = {ATOM_TYPE: element_type, ATOM_COORDS: atom_xyz}
                else:
                    line_struct = [line_head, atom_num, atom_type, res_type, mol_num, pdb_x, pdb_y, pdb_z, last_cols]
                    # noinspection PyTypeChecker
                    pdb_data[SEC_ATOMS].append(line_struct)
            elif line_head == 'END':
                pdb_data[SEC_TAIL].append(line)
                break
            # tail_content to contain everything after the 'Atoms' section
            else:
                # noinspection PyTypeChecker
                pdb_data[SEC_TAIL].append(line)
    pdb_data[NUM_ATOMS] = len(pdb_data[SEC_ATOMS])
    return pdb_data


def process_gausscom_file(gausscom_file):
    # Grabs and stores in gausscom_content as a dictionary with the keys:
    #    SEC_HEAD: header (route section, blank lines, comments, and full charge and multiplicity line)
    #    CHARGE: overall charge (only) as int
    #    MULT: overall multiplicity (only) as int
    #    SEC_ATOMS: atoms as a dict of dicts, with atom_id as key to dict with
    #        ATOM_TYPE: atom_type (str), ATOM_COORDS: (np array)
    #    SEC_TAIL: everything including and after the blank line following SEC_ATOMS
    with open(gausscom_file) as d:
        gausscom_content = {SEC_HEAD: [], SEC_ATOMS: {}, SEC_TAIL: [],
                            BASE_NAME: os.path.splitext(os.path.basename(gausscom_file))[0]}
        section = SEC_HEAD
        atom_id = 1
        blank_header_lines = 0

        for line in d:
            line = line.strip()

            if section == SEC_HEAD:
                gausscom_content[SEC_HEAD].append(line)
                if GAU_HEADER_PAT.match(line):
                    continue
                elif len(line) == 0:
                    blank_header_lines += 1
                    if blank_header_lines == 2:
                        section = SEC_ATOMS
                        line = next(d).strip()
                        gausscom_content[SEC_HEAD].append(line)
                        split_line = line.split()
                        try:
                            gausscom_content[CHARGE] = int(split_line[0])
                            gausscom_content[MULT] = int(split_line[1])
                        except (IndexError, ValueError):
                            raise InvalidDataError("Error in reading file {}\n  as a Gaussian input file. On the line "
                                                   "where charge and multiplicity are expected, "
                                                   "found: '{}'".format(gausscom_file, line))
                    continue

            elif section == SEC_ATOMS:
                if len(line) == 0:
                    section = SEC_TAIL
                    gausscom_content[SEC_TAIL].append(line)
                    continue
                split_line = line.split()

                atom_type = split_line[0]
                atom_xyz = np.array(list(map(float, split_line[1:4])))
                gausscom_content[SEC_ATOMS][atom_id] = {ATOM_TYPE: atom_type, ATOM_COORDS: atom_xyz}
                atom_id += 1
            elif section == SEC_TAIL:
                gausscom_content[SEC_TAIL].append(line)

    return gausscom_content


def process_gausslog_file(gausslog_file):
    # Grabs and stores in gausslog_content as a dictionary with the keys:
    #    (fyi: unlike process_gausscom_file, no SEC_HEAD is collected)
    #    CHARGE: overall charge (only) as int
    #    MULT: overall multiplicity (only) as int
    #    SEC_ATOMS: atoms as a dict of dicts, with atom_id as key to dict with
    #        ATOM_TYPE: atom_type (str), ATOM_COORDS: (np array)
    #    SEC_TAIL: everything including and after the blank line following SEC_ATOMS
    with open(gausslog_file) as d:
        gausslog_content = {SEC_ATOMS: {}, BASE_NAME: os.path.splitext(os.path.basename(gausslog_file))[0]}
        section = SEC_HEAD
        atom_id = 1

        for line in d:
            line = line.strip()

            if section == SEC_HEAD:
                # only get overall charge and mult
                if GAU_CHARGE_PAT.match(line):
                    split_line = line.split('=')
                    gausslog_content[CHARGE] = int(split_line[1].split()[0])
                    gausslog_content[MULT] = int(split_line[2].split()[0])
                    section = SEC_TAIL

            elif section == SEC_TAIL:
                if GAU_COORD_PAT.match(line):
                    next(d)
                    next(d)
                    section = SEC_ATOMS

            elif section == SEC_ATOMS:
                # will keep overwriting coordinates until it gets to the end
                while not GAU_SEP_PAT.match(line):
                    split_line = line.split()
                    atom_type = ATOM_NUM_DICT[int(split_line[1])]
                    atom_xyz = np.array(list(map(float, split_line[3:6])))
                    gausslog_content[SEC_ATOMS][atom_id] = {ATOM_TYPE: atom_type, ATOM_COORDS: atom_xyz}
                    atom_id += 1
                    line = next(d).strip()
                section = SEC_TAIL
                atom_id = 1

    return gausslog_content


def longest_common_substring(s1, s2):
    """
    From https://en.wikibooks.org/wiki/Algorithm_Implementation/Strings/Longest_common_substring#Python
    @param s1: string 1
    @param s2: string 2
    @return: string: the longest common string!
    """
    # noinspection PyUnusedLocal
    m = [[0] * (1 + len(s2)) for i in range(1 + len(s1))]
    longest, x_longest = 0, 0
    for x in range(1, 1 + len(s1)):
        for y in range(1, 1 + len(s2)):
            if s1[x - 1] == s2[y - 1]:
                m[x][y] = m[x - 1][y - 1] + 1
                if m[x][y] > longest:
                    longest = m[x][y]
                    x_longest = x
            else:
                m[x][y] = 0
    return s1[x_longest - longest: x_longest]


# FIGURES

def save_figure(name, save_fig=True):
    """
    Specifies where and if to save a created figure
    :param name: Name for the file
    :param save_fig: boolean as to whether to save fig; defaults to true (specify False if not desired)
    :return: n/a
    """
    if save_fig:
        plt.savefig(name, bbox_inches='tight', transparent=True)


def make_fig(name, x_array, y1_array, y1_label="", ls1="-", color1="blue",
             x2_array=None, y2_array=None, y2_label="", ls2='--', color2='orange',
             x3_array=None, y3_array=None, y3_label="", ls3=':',
             x4_array=None, y4_array=None, y4_label="", ls4='-.',
             x5_array=None, y5_array=None, y5_label="", ls5='-', color4='red',
             x_fill=None, y_fill=None, x2_fill=None, y2_fill=None,
             fill1_label=None, fill2_label=None,
             fill_color_1="green", fill_color_2="blue",
             x_label="", y_label="", x_lima=None, x_limb=None, y_lima=None, y_limb=None, loc=0,
             fig_width=DEF_FIG_WIDTH, fig_height=DEF_FIG_HEIGHT, axis_font_size=DEF_AXIS_SIZE,
             tick_font_size=DEF_TICK_SIZE, hide_x=False):
    """
    Many defaults to it is easy to adjust
    """
    # rc('text', usetex=True)
    # a general purpose plotting routine; can plot between 1 and 5 curves
    fig, ax = plt.subplots(figsize=(fig_width, fig_height))
    ax.plot(x_array, y1_array, ls1, label=y1_label, linewidth=2, color=color1)
    if y2_array is not None:
        if x2_array is None:
            x2_array = x_array
        ax.plot(x2_array, y2_array, label=y2_label, ls=ls2, linewidth=2, color=color2)
    if y3_array is not None:
        if x3_array is None:
            x3_array = x_array
        ax.plot(x3_array, y3_array, label=y3_label, ls=ls3, linewidth=3, color='green')
    if y4_array is not None:
        if x4_array is None:
            x4_array = x_array
        ax.plot(x4_array, y4_array, label=y4_label, ls=ls4, linewidth=3, color=color4)
    if y5_array is not None:
        if x5_array is None:
            x5_array = x_array
        ax.plot(x5_array, y5_array, label=y5_label, ls=ls5, linewidth=3, color='purple')
    ax.set_xlabel(x_label, fontsize=axis_font_size)
    ax.set_ylabel(y_label, fontsize=axis_font_size)
    if x_limb is not None:
        if x_lima is None:
            x_lima = 0.0
        ax.set_xlim([x_lima, x_limb])

    if y_limb is not None:
        if y_lima is None:
            y_lima = 0.0
        ax.set_ylim([y_lima, y_limb])

    if x_fill is not None:
        plt.fill_between(x_fill, y_fill, 0, color=fill_color_1, alpha='0.75')

    if x2_fill is not None:
        plt.fill_between(x2_fill, y2_fill, 0, color=fill_color_2, alpha='0.5')

    ax.tick_params(labelsize=tick_font_size)
    ax.xaxis.set_minor_locator(AutoMinorLocator(5))
    ax.yaxis.set_minor_locator(AutoMinorLocator(5))
    if len(y1_label) > 0:
        ax.legend(loc=loc, fontsize=tick_font_size, )
    if fill1_label and fill2_label:
        p1 = Rectangle((0, 0), 1, 1, fc=fill_color_1, alpha=0.75)
        p2 = Rectangle((0, 0), 1, 1, fc=fill_color_2, alpha=0.5)
        ax.legend([p1, p2], [fill1_label, fill2_label], loc=loc, fontsize=tick_font_size, )
    if hide_x:
        ax.xaxis.set_visible(False)
    else:
        ax.xaxis.grid(True, 'minor')
        ax.xaxis.grid(True, 'major')
    # ax.yaxis.grid(True, 'minor')
    # ax.yaxis.grid(True, 'major', linewidth=1)
    save_figure(name)
