#!/usr/bin/env python
"""

"""

import os
import sys
import argparse
import numpy as np
from common_wrangler.common import (InvalidDataError, warning,
                                    GOOD_RET, INPUT_ERROR, IO_ERROR, INVALID_DATA,
                                    EHPART_TO_KCAL_MOL, DEF_FIG_HEIGHT, DEF_FIG_WIDTH,
                                    create_out_fname, make_fig)
from gaussian_wrangler import __version__


__author__ = 'hmayes'

# Constants #


# Config keys
DEF_Y_LABEL = '\u0394G at {} K (kcal/mol)'


def parse_cmdline(argv):
    """
    Returns the parsed argument list and return code.
    `argv` is a list of arguments, or `None` for ``sys.argv[1:]``.
    """
    if argv is None:
        argv = sys.argv[1:]

    # initialize the parser object:
    parser = argparse.ArgumentParser(description='Creates a plot showing given values as horizontal lines, connected '
                                                 'by diagonal lines, as in a \u0394G or \u0394H plot.')
    parser.add_argument("-c", "--conv", help="Flag to convert values from a.u. to kcal/mol. The default is False.",
                        action='store_true')
    parser.add_argument("-d", "--out_dir", help="A directory where output files should be saved. The default location "
                                                "is the current working directory.", default=None)
    parser.add_argument("-l", "--list", help="The location of the list of values (with labels) to plot.",
                        default=None)
    parser.add_argument("-t", "--temp", help="Temperature in K for the plot of \u0394G.", default=None)
    parser.add_argument("-o", "--output_fname", help="The name of the output file to be created. The default is the "
                                                     "same base name as the list, with the '.png' extension.",
                        default=None)
    parser.add_argument("-fh", "--fig_height", help="Figure height in inches. "
                                                    "The default is {} in.".format(DEF_FIG_HEIGHT), default=None)
    parser.add_argument("-fw", "--fig_width", help="Figure width in inches. "
                                                   "The default is {} in.".format(DEF_FIG_WIDTH), default=None)
    parser.add_argument("-y", "--y_axis_label", help="Text for the y-axis label. The default is: {}.\n"
                                                     "Be sure to include braces ({{}}) for the temperature to be "
                                                     "filled in.".format(DEF_Y_LABEL.format("(input 'temp')")),
                        default=None)

    args = None
    try:
        args = parser.parse_args(argv)
        if not args.list:
            raise InvalidDataError("A list of data must be supplied.")
        if not args.out_dir:
            args.out_dir = os.getcwd()
        # user can define a new directory as the output directory
        if not os.path.exists(args.out_dir):
            os.makedirs(args.out_dir)

    except (SystemExit, InvalidDataError) as e:
        if hasattr(e, 'code') and e.code == 0:
            return args, GOOD_RET
        warning(e)
        parser.print_help()
        return args, INPUT_ERROR

    return args, GOOD_RET


def plot_delta_g(fname, g_temp, data_list, convert_flag, fig_width, fig_height, y_label):
    """
    Makes a plot of delta G at the specified temp
    :param fname: string, to save plot
    :param g_temp: float, temp at which delta Gs were calculated
    :param data_list: list of data, starting with the label
    :param convert_flag: Boolean on whether to convert from a.u. to kcal/mol
    :param fig_width: None or string; if none use default, otherwise make is a float
    :param fig_height: None or string; if none use default, otherwise make is a float
    :param y_label: None or string; if none use default
    :return: nothing, just save
    """
    max_y_lines = 5
    x_axis = []
    for index in range(len(data_list[0]) - 1):
        x_axis += [index * 3, index * 3+1]
    x_axis = np.array(x_axis)
    y_axis = []
    y_labels = []
    for index in range(max_y_lines):
        try:
            current_row = data_list[index]
            y_labels.append(current_row[0])
            plot_data = np.array([float(x) for x in current_row[1:]])
            if convert_flag:
                plot_data *= EHPART_TO_KCAL_MOL
            plot_data -= plot_data[0]
            y_data = []
            for value_index in range(len(plot_data)):
                y_data += [plot_data[value_index], plot_data[value_index]]
            y_axis.append(np.array(y_data))
        except IndexError:
            y_labels.append(None)
            y_axis.append(None)

    if fig_width:
        fig_width = float(fig_width)
    else:
        fig_width = DEF_FIG_WIDTH

    if fig_height:
        fig_height = float(fig_height)
    else:
        fig_height = DEF_FIG_HEIGHT

    if not y_label:
        y_label = DEF_Y_LABEL

    make_fig(fname, x_axis, y_axis[0],
             x_label='reaction coordinate', y_label=y_label.format(g_temp),
             y1_label=y_labels[0], y2_label=y_labels[1], y3_label=y_labels[2], y4_label=y_labels[3],
             y5_label=y_labels[4], y2_array=y_axis[1], y3_array=y_axis[2], y4_array=y_axis[3], y5_array=y_axis[4],
             ls2='-', ls3='-', ls4='-', ls5='-', fig_width=fig_width, fig_height=fig_height,
             # y_lima=y_min, y_limb=y_max,
             hide_x=True,
             )


def main(argv=None):
    print(f"Running GaussianWrangler script plot_steps version {__version__}")

    # Read input
    args, ret = parse_cmdline(argv)
    if ret != GOOD_RET or args is None:
        return ret

    try:
        # Make a list of lists from the input file list
        with open(args.list) as f:
            row_list = [row.strip().split() for row in f.readlines()]
            row_list = list(filter(None, row_list))

        if args.output_fname:
            plot_fname = create_out_fname(args.output_fname, base_dir=args.out_dir, ext='.png')
        else:
            plot_fname = create_out_fname(args.list, base_dir=args.out_dir, ext='.png')
        plot_delta_g(plot_fname, args.temp, row_list, args.conv, args.fig_width, args.fig_height, args.y_axis_label)
        print("Wrote file: {}".format(plot_fname))

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
