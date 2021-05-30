# coding=utf-8

"""
Moved methods here to make no quite so long
"""
import os
import sys
import matplotlib.pyplot as plt
import matplotlib.path as m_path
import matplotlib.patches as m_patches
import numpy as np
from glob import glob
from collections import namedtuple
from gaussian_wrangler.vib_scale_factors import (get_free_space, SCALING_DATA, SCALING_REFS)
from common_wrangler.common import (warning,
                                    AU_TO_J, EHPART_TO_KCAL_MOL, GAS_CONSTANT,
                                    InvalidDataError)

# # To make the output exactly match Gaussian's, use the values below instead importing them from common_wrangler.common
# GAS_CONSTANT = 8.31441  # J / K / mol; in common, GAS_CONSTANT = 8.314462618
# EHPART_TO_KCAL_MOL = 627.5095  # kcal/mol/(Eh/part); in common, the value is 627.5094709
# AU_TO_J = 4.184 * EHPART_TO_KCAL_MOL * 1000.0  # This value changes based on which EHPART_TO_KCAL_MOL is used


# Some literature references
GRIMME_REF = "Grimme, S. Chem. Eur. J. 2012, 18, 9955-9964"
TRUHLAR_REF = "Ribeiro, R. F.; Marenich, A. V.; Cramer, C. J.; Truhlar, D. G. J. Phys. Chem. B 2011, 115, 14556-14562"
HEAD_GORDON_REF = "Li, Y.; Gomes, J.; Sharada, S. M.; Bell, A. T.; Head-Gordon, M. J. Phys. Chem. C " \
                  "2015, 119, 1840-1850"
D3_REF = "Grimme, S.; Atony, J.; Ehrlich S.; Krieg, H. J. Chem. Phys. 2010, 132, 154104"
D3BJ_REF = "Grimme S.; Ehrlich, S.; Goerigk, L. J. Comput. Chem. 2011, 32, 1456-1465"
ATM_REF = "Axilrod, B. M.; Teller, E. J. Chem. Phys. 1943, 11, 299 \n Muto, Y. Proc. Phys. np. Soc. Jpn. " \
          "1944, 17, 629"
ONIOM_SCALE_REF = "Simon, L.; Paton, R. S. J. Am. Chem. Soc. 2018, 140, 5412-5420"
CSD_REF = ("C. R. Groom, I. J. Bruno, M. P. Lightfoot and S. C. Ward, Acta Cryst. 2016, B72, 171-179\n   "
           "Cordero, B.; Gomez V.; Platero-Prats, A. E.; Reves, M.; Echeverria, J.; Cremades, E.; Barragan, F.; "
           "Alvarez, S. Dalton Trans. 2008, 2832-2838")

ALPHABET = 'abcdefghijklmnopqrstuvwxyz'


def print_check_fails(check_attribute, file, attribute, option2=None):
    # Function for printing unique checks
    unique_attr = {}
    for i, attr in enumerate(check_attribute):
        if option2:
            attr = (attr, option2[i])
        if attr not in unique_attr:
            unique_attr[attr] = [file[i]]
        else:
            unique_attr[attr].append(file[i])
    warning("Different {} found: ".format(attribute))
    for attr in unique_attr:
        if option2:
            if float(attr[0]) < 0:
                atr_str = '       {} {}:'.format(attr[0], attr[1])
            else:
                atr_str = '        {} {}:'.format(attr[0], attr[1])
        else:
            atr_str = '        -{}: '.format(attr)
        for filename in unique_attr[attr]:
            rel_path = os.path.relpath(filename)
            if filename == unique_attr[attr][-1]:
                warning('{} {}'.format(atr_str, rel_path))
            else:
                warning('{} {}, '.format(atr_str, rel_path))


def all_same(items):
    return all(x == items[0] for x in items)


def create_plot(options, thermo_data):
    for key in thermo_data:
        if not hasattr(thermo_data[key], "qh_gibbs_free_energy"):
            raise InvalidDataError("\nWarning! Could not find thermodynamic data for " + key + "\n")
        if not hasattr(thermo_data[key], "sp_energy") and options.spc:
            raise InvalidDataError("\nWarning! Could not find thermodynamic data for " + key + "\n")
    graph_data = GetPES(options.graph, thermo_data, options.temperature, options.gconf, options.qh)
    graph_reaction_profile(graph_data, options)


# Scatter points that may overlap when graphing
def jitter(data_sets, color, ax, nx, marker, edge_color='black'):
    for i, p in enumerate(data_sets):
        y = [p]
        x = np.random.normal(nx, 0.015, size=len(y))
        ax.plot(x, y, alpha=0.5, markersize=7, color=color, marker=marker, markeredgecolor=edge_color,
                markeredgewidth=1, linestyle='None')


def graph_reaction_profile(graph_data, options):
    # Graph a reaction profile
    print("\n   Graphing Reaction Profile\n")
    data = {}
    # Get PES data
    for i, d_path in enumerate(graph_data.path):
        g_data = []
        zero_val = graph_data.qhg_zero[i][0]
        for j, e_abs in enumerate(graph_data.e_abs[i]):
            species = graph_data.qhg_abs[i][j]
            relative = species - zero_val
            if graph_data.units == 'kJ/mol':
                formatted_g = AU_TO_J / 1000.0 * relative
            else:
                formatted_g = EHPART_TO_KCAL_MOL * relative  # Defaults to kcal/mol
            g_data.append(formatted_g)
        data[d_path] = g_data

    # Grab any additional formatting for graph
    with open(options.graph) as f:
        yaml = f.readlines()
    ylim, color, show_conf, show_gconf, show_title = None, None, True, False, True
    label_point, label_xaxis, dpi, dec, legend, colors, gridlines, title = True, True, False, 2, True, None, False, None
    for i, y_line in enumerate(yaml):
        if y_line.strip().find('FORMAT') > -1:
            for j, s_line in enumerate(yaml[i + 1:]):
                s_line = s_line.strip()
                if s_line.find('ylim') > -1:
                    try:
                        ylim = s_line.replace(':', '=').split("=")[1].replace(' ', '').strip().split(',')
                    except IndexError:
                        pass
                if s_line.find('color') > -1:
                    try:
                        colors = s_line.replace(':', '=').split("=")[1].replace(' ', '').strip().split(',')
                    except IndexError:
                        pass
                if s_line.find('title') > -1:
                    try:
                        title_input = s_line.replace(':', '=').split("=")[1].strip().split(',')[0]
                        if title_input == 'false' or title_input == 'False':
                            show_title = False
                        else:
                            title = title_input
                    except IndexError:
                        pass
                if s_line.find('dec') > -1:
                    try:
                        dec = int(s_line.replace(':', '=').split("=")[1].strip().split(',')[0])
                    except IndexError:
                        pass
                if s_line.find('pointlabel') > -1:
                    try:
                        label_input = s_line.replace(':', '=').split("=")[1].strip().split(',')[0].lower()
                        if label_input == 'false':
                            label_point = False
                    except IndexError:
                        pass
                if s_line.find('show_conformers') > -1:
                    try:
                        conformers = s_line.replace(':', '=').split("=")[1].strip().split(',')[0].lower()
                        if conformers == 'false':
                            show_conf = False
                    except IndexError:
                        pass
                if s_line.find('show_gconf') > -1:
                    try:
                        gconf_input = s_line.replace(':', '=').split("=")[1].strip().split(',')[0].lower()
                        if gconf_input == 'true':
                            show_gconf = True
                    except IndexError:
                        pass
                if s_line.find('xlabel') > -1:
                    try:
                        label_input = s_line.replace(':', '=').split("=")[1].strip().split(',')[0].lower()
                        if label_input == 'false':
                            label_xaxis = False
                    except IndexError:
                        pass
                if s_line.find('dpi') > -1:
                    try:
                        dpi = int(s_line.replace(':', '=').split("=")[1].strip().split(',')[0])
                    except IndexError:
                        pass
                if s_line.find('legend') > -1:
                    try:
                        legend_input = s_line.replace(':', '=').split("=")[1].strip().split(',')[0].lower()
                        if legend_input == 'false':
                            legend = False
                    except IndexError:
                        pass
                if s_line.find('gridlines') > -1:
                    try:
                        gridline_input = s_line.replace(':', '=').split("=")[1].strip().split(',')[0].lower()
                        if gridline_input == 'true':
                            gridlines = True
                    except IndexError:
                        pass
    # Do some graphing
    g_path = m_path.Path
    fig, ax = plt.subplots()
    for i, d_path in enumerate(graph_data.path):
        for j in range(len(data[d_path]) - 1):
            if colors is None:
                color = 'k'
                colors = ['k']
            else:
                if len(colors) > 1:
                    color = colors[i]
                else:
                    color = colors[0]
            if j == 0:
                path_patch = m_patches.PathPatch(g_path([(j, data[d_path][j]), (j + 0.5, data[d_path][j]),
                                                         (j + 0.5, data[d_path][j + 1]), (j + 1, data[d_path][j + 1])],
                                                        [g_path.MOVETO, g_path.CURVE4, g_path.CURVE4, g_path.CURVE4]),
                                                 label=d_path, fc="none", transform=ax.transData, color=color)
            else:
                path_patch = m_patches.PathPatch(
                    g_path([(j, data[d_path][j]), (j + 0.5, data[d_path][j]), (j + 0.5, data[d_path][j + 1]),
                            (j + 1, data[d_path][j + 1])],
                           [g_path.MOVETO, g_path.CURVE4, g_path.CURVE4, g_path.CURVE4]),
                    fc="none", transform=ax.transData, color=color)
            ax.add_patch(path_patch)
            plt.hlines(data[d_path][j], j - 0.15, j + 0.15)
        plt.hlines(data[d_path][-1], len(data[d_path]) - 1.15, len(data[d_path]) - 0.85)

    if show_conf:
        markers = ['o', 's', 'x', 'P', 'D']
        for i in range(len(graph_data.g_qhgvals)):  # i = reaction pathways
            for j in range(len(graph_data.g_qhgvals[i])):  # j = reaction steps
                for k in range(len(graph_data.g_qhgvals[i][j])):  # k = species
                    zero_val = graph_data.g_species_qhgzero[i][j][k]
                    points = graph_data.g_qhgvals[i][j][k]
                    points[:] = [((x - zero_val) + (graph_data.qhg_abs[i][j] - graph_data.qhg_zero[i][0]) + (
                            graph_data.g_rel_val[i][j] - graph_data.qhg_abs[i][j])) * EHPART_TO_KCAL_MOL for
                                 x in points]
                    if len(colors) > 1:
                        jitter(points, colors[i], ax, j, markers[k])
                    else:
                        jitter(points, color, ax, j, markers[k])
                    if show_gconf:
                        plt.hlines((graph_data.g_rel_val[i][j] - graph_data.qhg_zero[i][0]) * EHPART_TO_KCAL_MOL,
                                   j - 0.15, j + 0.15, linestyles='dashed')

    # Annotate points with energy level
    if label_point:
        for d_path in graph_data.path:
            for i, point in enumerate(data[d_path]):
                if dec == 1:
                    ax.annotate("{:.1f}".format(point), (i, point - fig.get_figheight() * fig.dpi * 0.025),
                                horizontalalignment='center')
                else:
                    ax.annotate("{:.2f}".format(point), (i, point - fig.get_figheight() * fig.dpi * 0.025),
                                horizontalalignment='center')
    if ylim is not None:
        ax.set_ylim(float(ylim[0]), float(ylim[1]))
    if show_title:
        if title is None:
            ax.set_title("Reaction Profile")
        else:
            ax.set_title(title)

    ax.set_ylabel(r"$G_{rel}$ (kcal / mol)")
    plt.minorticks_on()
    ax.tick_params(axis='x', which='minor', bottom=False)
    ax.tick_params(which='minor', labelright=True, right=True)
    ax.tick_params(labelright=True, right=True)
    if gridlines:
        ax.yaxis.grid(linestyle='--', linewidth=0.5)
        ax.xaxis.grid(linewidth=0)
    ax_label = []
    xaxis_text = []
    new_ax_text_list = []
    for i, d_path in enumerate(graph_data.path):
        new_ax_text = []
        ax_label.append(d_path)
        for j, e_abs in enumerate(graph_data.e_abs[i]):
            if i == 0:
                xaxis_text.append(graph_data.species[i][j])
            else:
                new_ax_text.append(graph_data.species[i][j])
        new_ax_text_list.append(new_ax_text)
    # Label rxn steps
    if label_xaxis:
        if colors is None:
            plt.xticks(range(len(xaxis_text)), xaxis_text, color='k')
        else:
            plt.xticks(range(len(xaxis_text)), xaxis_text, color=colors[0])
        locations, labels = plt.xticks()
        new_ax = []
        for i in range(len(ax_label)):
            if i > 0:
                y = ax.twiny()
                new_ax.append(y)
        for i in range(len(new_ax)):
            new_ax[i].set_xticks(locations)
            new_ax[i].set_xlim(ax.get_xlim())
            if len(colors) > 1:
                new_ax[i].tick_params(axis='x', colors=colors[i + 1])
            else:
                new_ax[i].tick_params(axis='x', colors='k')
            new_ax[i].set_xticklabels(new_ax_text_list[i + 1])
            new_ax[i].xaxis.set_ticks_position('bottom')
            new_ax[i].xaxis.set_label_position('bottom')
            new_ax[i].xaxis.set_ticks_position('none')
            new_ax[i].spines['bottom'].set_position(('outward', 15 * (i + 1)))
            new_ax[i].spines['bottom'].set_visible(False)
    else:
        plt.xticks(range(len(xaxis_text)))
        ax.xaxis.set_ticklabels([])
    if legend:
        plt.legend()
    if dpi:
        plt.savefig('Rxn_profile_' + options.graph.split('.')[0] + '.png', dpi=dpi)
    plt.show()


class GetPES:
    # Obtain relative thermochemistry between species and for reactions
    # noinspection PyTypeChecker
    def __init__(self, file, thermo_data, temperature, gconf, qh, cosmo=None):
        # Default values
        self.dec, self.units, self.boltz = 2, 'kcal/mol', False

        g_corr, qg_corr = None, None  # Make IDE happy

        with open(file) as f:
            data = f.readlines()
        folder, program, names, files, zeros, pes_list = None, None, [], [], [], []
        for i, d_line in enumerate(data):
            if d_line.strip().find('PES') > -1:
                for j, p_line in enumerate(data[i + 1:]):
                    if p_line.strip().startswith('#'):
                        pass
                    elif len(p_line) <= 2:
                        pass
                    elif p_line.strip().startswith('---'):
                        break
                    elif p_line.strip() != '':
                        pathway, pes = p_line.strip().replace(':', '=').split("=")
                        # Auto-grab first species as zero unless specified
                        pes_list.append(pes)
                        zeros.append(pes.strip().lstrip('[').rstrip(']').split(',')[0])
                        # Look at SPECIES block to determine filenames
            if d_line.strip().find('SPECIES') > -1:
                for j, s_line in enumerate(data[i + 1:]):
                    if s_line.strip().startswith('---'):
                        break
                    else:
                        if s_line.lower().strip().find('folder') < 0:
                            try:
                                n, f = (s_line.strip().replace(':', '=').split("="))
                                # Check the specified filename is also one that GoodVibes has thermochemistry for:
                                if f.find('*') == -1 and f not in pes_list:
                                    match = None
                                    for key in thermo_data:
                                        if os.path.splitext(os.path.basename(key))[0] in f.replace('[', ''). \
                                                replace(']', '').replace('+', ',').replace(' ', '').split(','):
                                            match = key
                                    if match:
                                        names.append(n.strip())
                                        files.append(match)
                                    else:
                                        warning("   Warning! " + f.strip() + ' is specified in ' + file +
                                                ' but no thermochemistry data found\n')
                                elif f not in pes_list:
                                    match = []
                                    for key in thermo_data:
                                        if os.path.splitext(os.path.basename(key))[0].find(f.strip().strip('*')) == 0:
                                            match.append(key)
                                    if len(match) > 0:
                                        names.append(n.strip())
                                        files.append(match)
                                    else:
                                        warning("   Warning! " + f.strip() + ' is specified in ' + file +
                                                ' but no thermochemistry data found\n')
                            except ValueError:
                                if s_line.isspace():
                                    pass
                                elif s_line.strip().find('#') > -1:
                                    pass
                                elif len(s_line) > 2:
                                    warning("   Warning! {} input is incorrectly formatted for s_line:\n    {}".
                                            format(os.path.relpath(file), s_line))
            # Look at FORMAT block to see if user has specified any formatting rules
            if d_line.strip().find('FORMAT') > -1:
                for j, f_line in enumerate(data[i + 1:]):
                    f_line = f_line.strip()
                    if f_line.find('dec') > -1:
                        try:
                            self.dec = int(f_line.replace(':', '=').split("=")[1].strip())
                        except IndexError:
                            pass
                    if f_line.find('zero') > -1:
                        zeros = []
                        try:
                            zeros.append(f_line.replace(':', '=').split("=")[1].strip())
                        except IndexError:
                            pass
                    if f_line.find('units') > -1:
                        try:
                            self.units = f_line.replace(':', '=').split("=")[1].strip()
                        except IndexError:
                            pass
                    if f_line.find('boltz') > -1:
                        try:
                            self.boltz = f_line.replace(':', '=').split("=")[1].strip()
                        except IndexError:
                            pass

        for i in range(len(files)):
            if len(files[i]) == 1:
                files[i] = files[i][0]
        species = dict(zip(names, files))
        self.path, self.species = [], []
        self.spc_abs, self.e_abs, self.zpe_abs, self.h_abs, self.qh_abs, self.s_abs, self.qs_abs, self.g_abs, \
            self.qhg_abs, self.cosmo_qhg_abs = [], [], [], [], [], [], [], [], [], []
        self.spc_zero, self.e_zero, self.zpe_zero, self.h_zero, self.qh_zero, self.ts_zero, self.qhts_zero, \
            self.g_zero, self.qhg_zero, self.cosmo_qhg_zero = [], [], [], [], [], [], [], [], [], []
        self.g_qhgvals, self.g_species_qhgzero, self.g_rel_val = [], [], []
        structure = None
        # Loop over .yaml file, grab energies, populate arrays and compute Boltzmann factors
        with open(file) as f:
            data = f.readlines()
        for i, d_line in enumerate(data):
            if d_line.strip().find('PES') > -1:
                n = 0
                for j, p_line in enumerate(data[i + 1:]):
                    if p_line.strip().startswith('#'):
                        pass
                    elif len(p_line) <= 2:
                        pass
                    elif p_line.strip().startswith('---'):
                        break
                    elif p_line.strip() != '':
                        try:
                            self.e_zero.append([])
                            self.spc_zero.append([])
                            self.zpe_zero.append([])
                            self.h_zero.append([])
                            self.qh_zero.append([])
                            self.ts_zero.append([])
                            self.qhts_zero.append([])
                            self.g_zero.append([])
                            self.qhg_zero.append([])
                            self.cosmo_qhg_zero.append([])
                            min_conf = False
                            spc_zero, e_zero, zpe_zero, h_zero, qh_zero, s_zero, qs_zero, g_zero, qhg_zero = \
                                0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
                            h_conf, h_tot, s_conf, s_tot, qh_conf, qh_tot, qs_conf, qs_tot, cosmo_qhg_zero = \
                                0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
                            zero_structures = zeros[n].replace(' ', '').split('+')
                            # Routine for 'zero' values
                            for structure in zero_structures:
                                try:
                                    if not isinstance(species[structure], list):
                                        if hasattr(thermo_data[species[structure]], "sp_energy"):
                                            spc_zero += thermo_data[species[structure]].sp_energy
                                        e_zero += thermo_data[species[structure]].scf_energy
                                        zpe_zero += thermo_data[species[structure]].zpe
                                        h_zero += thermo_data[species[structure]].enthalpy
                                        qh_zero += thermo_data[species[structure]].qh_enthalpy
                                        s_zero += thermo_data[species[structure]].entropy
                                        qs_zero += thermo_data[species[structure]].qh_entropy
                                        g_zero += thermo_data[species[structure]].gibbs_free_energy
                                        qhg_zero += thermo_data[species[structure]].qh_gibbs_free_energy
                                        cosmo_qhg_zero += thermo_data[species[structure]].cosmo_qhg
                                    else:  # If we have a list of different kinds of structures: loop over conformers
                                        g_min, boltz_sum = sys.float_info.max, 0.0
                                        # Find minimum G, along with associated enthalpy and entropy
                                        for conformer in species[structure]:
                                            if cosmo:
                                                if thermo_data[conformer].cosmo_qhg <= g_min:
                                                    min_conf = thermo_data[conformer]
                                                    g_min = thermo_data[conformer].cosmo_qhg
                                            else:
                                                if thermo_data[conformer].qh_gibbs_free_energy <= g_min:
                                                    min_conf = thermo_data[conformer]
                                                    g_min = thermo_data[conformer].qh_gibbs_free_energy
                                        for conformer in species[structure]:  # Get a Boltzmann sum for conformers
                                            if cosmo:
                                                g_rel = thermo_data[conformer].cosmo_qhg - g_min
                                            else:
                                                g_rel = thermo_data[conformer].qh_gibbs_free_energy - g_min
                                            boltz_fac = np.exp(-g_rel * AU_TO_J / GAS_CONSTANT / temperature)
                                            boltz_sum += boltz_fac
                                        for conformer in species[structure]:
                                            # Calculate relative data based on Gmin and the Boltzmann sum
                                            if cosmo:
                                                g_rel = thermo_data[conformer].cosmo_qhg - g_min
                                            else:
                                                g_rel = thermo_data[conformer].qh_gibbs_free_energy - g_min
                                            boltz_fac = np.exp(-g_rel * AU_TO_J / GAS_CONSTANT / temperature)
                                            boltz_prob = boltz_fac / boltz_sum
                                            if hasattr(thermo_data[conformer], "sp_energy"):
                                                if thermo_data[conformer].sp_energy == '!':
                                                    raise InvalidDataError("Not all files contain a SPC value; "
                                                                           "relative values will not be calculated.")
                                                else:
                                                    spc_zero += thermo_data[conformer].sp_energy * boltz_prob
                                            e_zero += thermo_data[conformer].scf_energy * boltz_prob
                                            zpe_zero += thermo_data[conformer].zpe * boltz_prob
                                            if gconf:  # Default calculate gconf correction for conformers
                                                h_conf += thermo_data[conformer].enthalpy * boltz_prob
                                                s_conf += thermo_data[conformer].entropy * boltz_prob
                                                s_conf += -GAS_CONSTANT / AU_TO_J * boltz_prob * np.log(boltz_prob)

                                                qh_conf += thermo_data[conformer].qh_enthalpy * boltz_prob
                                                qs_conf += thermo_data[conformer].qh_entropy * boltz_prob
                                                qs_conf += -GAS_CONSTANT / AU_TO_J * boltz_prob * np.log(boltz_prob)
                                            else:
                                                h_zero += thermo_data[conformer].enthalpy * boltz_prob
                                                s_zero += thermo_data[conformer].entropy * boltz_prob
                                                g_zero += thermo_data[conformer].gibbs_free_energy * boltz_prob

                                                qh_zero += thermo_data[conformer].qh_enthalpy * boltz_prob
                                                qs_zero += thermo_data[conformer].qh_entropy * boltz_prob
                                                qhg_zero += thermo_data[conformer].qh_gibbs_free_energy * boltz_prob
                                                cosmo_qhg_zero += thermo_data[conformer].cosmo_qhg * boltz_prob

                                        if gconf:
                                            h_adj = h_conf - min_conf.enthalpy
                                            h_tot = min_conf.enthalpy + h_adj
                                            s_adj = s_conf - min_conf.entropy
                                            s_tot = min_conf.entropy + s_adj
                                            g_corr = h_tot - temperature * s_tot
                                            qh_adj = qh_conf - min_conf.qh_enthalpy
                                            qh_tot = min_conf.qh_enthalpy + qh_adj
                                            qs_adj = qs_conf - min_conf.qh_entropy
                                            qs_tot = min_conf.qh_entropy + qs_adj
                                            if qh:
                                                qg_corr = qh_tot - temperature * qs_tot
                                            else:
                                                qg_corr = h_tot - temperature * qs_tot
                                except KeyError:
                                    warning("Structure {structure} has not been defined correctly as "
                                            "energy-zero in file\n   Make sure this structure matches one of the "
                                            "SPECIES defined in the same file\n")
                                    raise InvalidDataError("   Please edit " + file + " and try again\n")
                            # Set zero vals here
                            conformers, single_structure, mix = False, False, False
                            for structure in zero_structures:
                                if not isinstance(species[structure], list):
                                    single_structure = True
                                else:
                                    conformers = True
                            if conformers and single_structure:
                                mix = True
                            if gconf and min_conf:
                                if mix:
                                    h_mix = h_tot + h_zero
                                    s_mix = s_tot + s_zero
                                    g_mix = g_corr + g_zero
                                    qh_mix = qh_tot + qh_zero
                                    qs_mix = qs_tot + qs_zero
                                    qg_mix = qg_corr + qhg_zero
                                    cosmo_qhg_mix = qg_corr + cosmo_qhg_zero
                                    self.h_zero[n].append(h_mix)
                                    self.ts_zero[n].append(s_mix)
                                    self.g_zero[n].append(g_mix)
                                    self.qh_zero[n].append(qh_mix)
                                    self.qhts_zero[n].append(qs_mix)
                                    self.qhg_zero[n].append(qg_mix)
                                    self.cosmo_qhg_zero[n].append(cosmo_qhg_mix)
                                elif conformers:
                                    self.h_zero[n].append(h_tot)
                                    self.ts_zero[n].append(s_tot)
                                    self.g_zero[n].append(g_corr)
                                    self.qh_zero[n].append(qh_tot)
                                    self.qhts_zero[n].append(qs_tot)
                                    self.qhg_zero[n].append(qg_corr)
                                    self.cosmo_qhg_zero[n].append(qg_corr)
                            else:
                                self.h_zero[n].append(h_zero)
                                self.ts_zero[n].append(s_zero)
                                self.g_zero[n].append(g_zero)

                                self.qh_zero[n].append(qh_zero)
                                self.qhts_zero[n].append(qs_zero)
                                self.qhg_zero[n].append(qhg_zero)
                                self.cosmo_qhg_zero[n].append(cosmo_qhg_zero)

                            self.spc_zero[n].append(spc_zero)
                            self.e_zero[n].append(e_zero)
                            self.zpe_zero[n].append(zpe_zero)

                            self.species.append([])
                            self.e_abs.append([])
                            self.spc_abs.append([])
                            self.zpe_abs.append([])
                            self.h_abs.append([])
                            self.qh_abs.append([])
                            self.s_abs.append([])
                            self.g_abs.append([])
                            self.qs_abs.append([])
                            self.qhg_abs.append([])
                            self.cosmo_qhg_abs.append([])
                            self.g_qhgvals.append([])
                            self.g_species_qhgzero.append([])
                            self.g_rel_val.append([])  # graphing

                            pathway, pes = p_line.strip().replace(':', '=').split("=")
                            pes = pes.strip()
                            points = [entry.strip() for entry in pes.lstrip('[').rstrip(']').split(',')]
                            self.path.append(pathway.strip())
                            # Obtain relative values for each species
                            for m, point in enumerate(points):
                                if point != '':
                                    # Create values to populate
                                    point_structures = point.replace(' ', '').split('+')
                                    e_abs, spc_abs, zpe_abs, h_abs, qh_abs, s_abs, g_abs, qs_abs, qhg_abs, \
                                        cosmo_qhg_abs = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
                                    qh_conf, qh_tot, qs_conf, qs_tot, h_conf, h_tot, s_conf, s_tot, g_corr, \
                                        qg_corr = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
                                    min_conf = False
                                    rel_val = 0.0
                                    self.g_qhgvals[n].append([])
                                    self.g_species_qhgzero[n].append([])
                                    try:
                                        # Loop over structures, structures are species specified
                                        for k, structure in enumerate(point_structures):
                                            zero_conf = 0.0
                                            self.g_qhgvals[n][m].append([])
                                            if not isinstance(species[structure], list):  # Only one conf in structures
                                                e_abs += thermo_data[species[structure]].scf_energy
                                                if hasattr(thermo_data[species[structure]], "sp_energy"):
                                                    spc_abs += thermo_data[species[structure]].sp_energy
                                                zpe_abs += thermo_data[species[structure]].zpe
                                                h_abs += thermo_data[species[structure]].enthalpy
                                                qh_abs += thermo_data[species[structure]].qh_enthalpy
                                                s_abs += thermo_data[species[structure]].entropy
                                                g_abs += thermo_data[species[structure]].gibbs_free_energy
                                                qs_abs += thermo_data[species[structure]].qh_entropy
                                                qhg_abs += thermo_data[species[structure]].qh_gibbs_free_energy
                                                cosmo_qhg_abs += thermo_data[species[structure]].cosmo_qhg
                                                zero_conf += thermo_data[species[structure]].qh_gibbs_free_energy
                                                self.g_qhgvals[n][m][k].append(
                                                    thermo_data[species[structure]].qh_gibbs_free_energy)
                                                rel_val += thermo_data[species[structure]].qh_gibbs_free_energy
                                            else:
                                                # If a list of different kinds of structures: loop over conformers
                                                g_min, boltz_sum = sys.float_info.max, 0.0
                                                # Find minimum G, along with associated enthalpy and entropy
                                                for conformer in species[structure]:
                                                    if cosmo:
                                                        if thermo_data[conformer].cosmo_qhg <= g_min:
                                                            min_conf = thermo_data[conformer]
                                                            g_min = thermo_data[conformer].cosmo_qhg
                                                    else:
                                                        if thermo_data[conformer].qh_gibbs_free_energy <= g_min:
                                                            min_conf = thermo_data[conformer]
                                                            g_min = thermo_data[conformer].qh_gibbs_free_energy
                                                # Get a Boltzmann sum for conformers
                                                for conformer in species[structure]:
                                                    if cosmo:
                                                        g_rel = thermo_data[conformer].cosmo_qhg - g_min
                                                    else:
                                                        g_rel = thermo_data[conformer].qh_gibbs_free_energy - g_min
                                                    boltz_fac = np.exp(-g_rel * AU_TO_J / GAS_CONSTANT / temperature)
                                                    boltz_sum += boltz_fac
                                                # Calculate relative data based on Gmin and the Boltzmann sum
                                                for conformer in species[structure]:
                                                    if cosmo:
                                                        g_rel = thermo_data[conformer].cosmo_qhg - g_min
                                                    else:
                                                        g_rel = thermo_data[conformer].qh_gibbs_free_energy - g_min
                                                    boltz_fac = np.exp(-g_rel * AU_TO_J / GAS_CONSTANT / temperature)
                                                    boltz_prob = boltz_fac / boltz_sum
                                                    if hasattr(thermo_data[conformer], "sp_energy"):
                                                        if thermo_data[conformer].sp_energy == '!':
                                                            raise InvalidDataError("\n   Not all files contain a SPC "
                                                                                   "value; relative values will not "
                                                                                   " be calculated.\n")
                                                        else:
                                                            spc_abs += thermo_data[conformer].sp_energy * boltz_prob
                                                    e_abs += thermo_data[conformer].scf_energy * boltz_prob
                                                    zpe_abs += thermo_data[conformer].zpe * boltz_prob
                                                    if cosmo:
                                                        zero_conf += thermo_data[conformer].cosmo_qhg * boltz_prob
                                                        rel_val += thermo_data[conformer].cosmo_qhg * boltz_prob
                                                    else:
                                                        zero_conf += thermo_data[
                                                                         conformer].qh_gibbs_free_energy * boltz_prob
                                                        rel_val += thermo_data[
                                                                       conformer].qh_gibbs_free_energy * boltz_prob
                                                    if gconf:  # Default calculate gconf correction for conformers
                                                        h_conf += thermo_data[conformer].enthalpy * boltz_prob
                                                        s_conf += thermo_data[conformer].entropy * boltz_prob
                                                        s_conf += -(GAS_CONSTANT / AU_TO_J) * (boltz_prob *
                                                                                               np.log(boltz_prob))
                                                        qh_conf += thermo_data[conformer].qh_enthalpy * boltz_prob
                                                        qs_conf += thermo_data[conformer].qh_entropy * boltz_prob
                                                        qs_conf += -GAS_CONSTANT / AU_TO_J * (boltz_prob *
                                                                                              np.log(boltz_prob))
                                                    else:
                                                        h_abs += thermo_data[conformer].enthalpy * boltz_prob
                                                        s_abs += thermo_data[conformer].entropy * boltz_prob
                                                        g_abs += thermo_data[conformer].gibbs_free_energy * boltz_prob

                                                        qh_abs += thermo_data[conformer].qh_enthalpy * boltz_prob
                                                        qs_abs += thermo_data[conformer].qh_entropy * boltz_prob
                                                        qhg_abs += thermo_data[
                                                                       conformer].qh_gibbs_free_energy * boltz_prob
                                                        cosmo_qhg_abs += thermo_data[conformer].cosmo_qhg * boltz_prob
                                                    if cosmo:
                                                        self.g_qhgvals[n][m][k].append(thermo_data[conformer].cosmo_qhg)
                                                    else:
                                                        self.g_qhgvals[n][m][k].append(
                                                            thermo_data[conformer].qh_gibbs_free_energy)
                                                if gconf:
                                                    h_adj = h_conf - min_conf.enthalpy
                                                    h_tot = min_conf.enthalpy + h_adj
                                                    s_adj = s_conf - min_conf.entropy
                                                    s_tot = min_conf.entropy + s_adj
                                                    g_corr = h_tot - temperature * s_tot
                                                    qh_adj = qh_conf - min_conf.qh_enthalpy
                                                    qh_tot = min_conf.qh_enthalpy + qh_adj
                                                    qs_adj = qs_conf - min_conf.qh_entropy
                                                    qs_tot = min_conf.qh_entropy + qs_adj
                                                    if qh:
                                                        qg_corr = qh_tot - temperature * qs_tot
                                                    else:
                                                        qg_corr = h_tot - temperature * qs_tot
                                            self.g_species_qhgzero[n][i].append(zero_conf)  # Raw data for graphing
                                    except KeyError:
                                        base_name = os.path.relpath(file)
                                        warning("Structure {} has not been defined correctly "
                                                "in {}\n".format(structure, base_name))
                                        raise InvalidDataError("   Please edit {} and try again\n".format(base_name))
                                    self.species[n].append(point)
                                    self.e_abs[n].append(e_abs)
                                    self.spc_abs[n].append(spc_abs)
                                    self.zpe_abs[n].append(zpe_abs)
                                    conformers, single_structure, mix = False, False, False
                                    self.g_rel_val[n].append(rel_val)
                                    for structure in point_structures:
                                        if not isinstance(species[structure], list):
                                            single_structure = True
                                        else:
                                            conformers = True
                                    if conformers and single_structure:
                                        mix = True
                                    if gconf and min_conf:
                                        if mix:
                                            h_mix = h_tot + h_abs
                                            s_mix = s_tot + s_abs
                                            g_mix = g_corr + g_abs
                                            qh_mix = qh_tot + qh_abs
                                            qs_mix = qs_tot + qs_abs
                                            qg_mix = qg_corr + qhg_abs
                                            cosmo_qhg_mix = qg_corr + cosmo_qhg_zero
                                            self.h_abs[n].append(h_mix)
                                            self.s_abs[n].append(s_mix)
                                            self.g_abs[n].append(g_mix)
                                            self.qh_abs[n].append(qh_mix)
                                            self.qs_abs[n].append(qs_mix)
                                            self.qhg_abs[n].append(qg_mix)
                                            self.cosmo_qhg_abs[n].append(cosmo_qhg_mix)
                                        elif conformers:
                                            self.h_abs[n].append(h_tot)
                                            self.s_abs[n].append(s_tot)
                                            self.g_abs[n].append(g_corr)
                                            self.qh_abs[n].append(qh_tot)
                                            self.qs_abs[n].append(qs_tot)
                                            self.qhg_abs[n].append(qg_corr)
                                            self.cosmo_qhg_abs[n].append(qg_corr)
                                    else:
                                        self.h_abs[n].append(h_abs)
                                        self.s_abs[n].append(s_abs)
                                        self.g_abs[n].append(g_abs)

                                        self.qh_abs[n].append(qh_abs)
                                        self.qs_abs[n].append(qs_abs)
                                        self.qhg_abs[n].append(qhg_abs)
                                        self.cosmo_qhg_abs[n].append(cosmo_qhg_abs)
                                else:
                                    self.species[n].append('none')
                                    self.e_abs[n].append(np.nan)

                            n = n + 1
                        except IndexError:
                            pass


def pes_options(e_sum, g_sum, h_sum, i, options, pes, qhg_sum, sels, delim_row, zero_vals):
    for j, e_abs in enumerate(pes.e_abs[i]):
        species = get_species(i, j, options.qh, options.temperature, pes)
        if options.cosmo:
            species.append(pes.cosmo_qhg_abs[i][j])
        relative = [species[x] - zero_vals[x] for x in range(len(zero_vals))]
        if pes.units == 'kJ/mol':
            formatted_list = [AU_TO_J / 1000.0 * x for x in relative]
        else:
            formatted_list = [EHPART_TO_KCAL_MOL * x for x in relative]  # Defaults to kcal/mol
        print("\no  ")
        if options.spc:
            if options.qh and options.cosmo:
                if pes.dec == 1:
                    print('{:<39} {:13.2f} {:13.1f} {:10.1f} {:13.1f} {:13.1f} {:10.1f} {:10.1f} '
                          '{:13.1f} {:13.1f} {:13.1f}'.format(pes.species[i][j], *formatted_list))
                if pes.dec == 2:
                    print('{:<39} {:13.2f} {:13.2f} {:10.2f} {:13.2f} {:13.2f} {:10.2f} {:10.2f} '
                          '{:13.2f} {:13.2f} {:13.2f}'.format(pes.species[i][j], *formatted_list))
            elif options.qh or options.cosmo:
                if pes.dec == 1:
                    print('{:<39} {:13.2f} {:13.1f} {:10.1f} {:13.1f} {:13.1f} {:10.1f} {:10.1f} '
                          '{:13.1f} {:13.1f}'.format(pes.species[i][j], *formatted_list))
                if pes.dec == 2:
                    print('{:<39} {:13.2f} {:13.2f} {:10.2f} {:13.2f} {:13.2f} {:10.2f} {:10.2f} '
                          '{:13.2f} {:13.2f}'.format(pes.species[i][j], *formatted_list))
            else:
                if pes.dec == 1:
                    print('{:<39} {:13.2f} {:13.1f} {:10.1f} {:13.1f} {:10.1f} {:10.1f} {:13.1f} '
                          '{:13.1f}'.format(pes.species[i][j], *formatted_list))
                if pes.dec == 2:
                    print('{:<39} {:13.2f} {:13.2f} {:10.2f} {:13.2f} {:10.2f} {:10.2f} {:13.2f} '
                          '{:13.2f}'.format(pes.species[i][j], *formatted_list))

        else:
            formatted_list = formatted_list[1:]
            if options.qh and options.cosmo:
                if pes.dec == 1:
                    print('{:<39} {:13.2f} {:10.1f} {:13.1f} {:13.1f} {:10.1f} {:10.1f} {:13.1f} '
                          '{:13.1f} {:13.1f}'.format(pes.species[i][j], *formatted_list))
                if pes.dec == 2:
                    print('{:<39} {:13.2f} {:10.2f} {:13.2f} {:13.2f} {:10.2f} {:10.2f} {:13.2f} '
                          '{:13.2f} {:13.2f}'.format(pes.species[i][j], *formatted_list))
            elif options.qh or options.cosmo:
                if pes.dec == 1:
                    print('{:<39} {:13.2f} {:10.1f} {:13.1f} {:13.1f} {:10.1f} {:10.1f} {:13.1f} '
                          '{:13.1f}'.format(pes.species[i][j], *formatted_list))
                if pes.dec == 2:
                    print('{:<39} {:13.2f} {:10.2f} {:13.2f} {:13.2f} {:10.2f} {:10.2f} {:13.2f} '
                          '{:13.2f}'.format(pes.species[i][j], *formatted_list))
            else:
                if pes.dec == 1:
                    print('{:<39} {:13.2f} {:10.1f} {:13.1f} {:10.1f} {:10.1f} {:13.1f} '
                          '{:13.1f}'.format(pes.species[i][j], *formatted_list))
                if pes.dec == 2:
                    print('{:<39} {:13.2f} {:10.2f} {:13.2f} {:10.2f} {:10.2f} {:13.2f} '
                          '{:13.2f}'.format(pes.species[i][j], *formatted_list))
        if pes.boltz:
            boltz = [np.exp(-relative[1] * AU_TO_J / GAS_CONSTANT / options.temperature) / e_sum,
                     np.exp(-relative[3] * AU_TO_J / GAS_CONSTANT / options.temperature) / h_sum,
                     np.exp(-relative[6] * AU_TO_J / GAS_CONSTANT / options.temperature) / g_sum,
                     np.exp(-relative[7] * AU_TO_J / GAS_CONSTANT / options.temperature) / qhg_sum]
            selectivity = [boltz[x] * 100.0 for x in range(len(boltz))]
            print("\n  " + '{:<39} {:13.2f}%{:24.2f}%{:35.2f}%{:13.2f}%'.format('', *selectivity))
            sels.append(selectivity)
    if pes.boltz == 'ee' and len(sels) == 2:
        ee = [sels[0][x] - sels[1][x] for x in range(len(sels[0]))]
        if options.spc:
            print("\n" + delim_row + "\n   " + '{:<39} {:27.2f} {:24.1f} {:35.1f} {:13.1f} '.
                  format('ee (%)', *ee))
        else:
            print("\n" + delim_row + "\n   " + '{:<39} {:13.2f}%{:24.1f}%{:35.1f}%{:13.1f}%'
                  .format('ee (%)', *ee))

    print("\n" + delim_row + "\n")


def output_rel_e_data(options, delim_row, thermo_data):
    # Making IDe happy
    e_sum, g_sum, h_sum, qhg_sum, sels = None, None, None, None, None

    if options.cosmo:
        pes = GetPES(options.pes, thermo_data, options.temperature, options.gconf, options.qh,
                     cosmo=True)
    else:
        pes = GetPES(options.pes, thermo_data, options.temperature, options.gconf, options.qh)
    for i, path in enumerate(pes.path):
        if options.qh:
            zero_vals = [pes.spc_zero[i][0], pes.e_zero[i][0], pes.zpe_zero[i][0], pes.h_zero[i][0],
                         pes.qh_zero[i][0], options.temperature * pes.ts_zero[i][0],
                         options.temperature * pes.qhts_zero[i][0], pes.g_zero[i][0], pes.qhg_zero[i][0]]
        else:
            zero_vals = [pes.spc_zero[i][0], pes.e_zero[i][0], pes.zpe_zero[i][0], pes.h_zero[i][0],
                         options.temperature * pes.ts_zero[i][0], options.temperature * pes.qhts_zero[i][0],
                         pes.g_zero[i][0], pes.qhg_zero[i][0]]
        if options.cosmo:
            zero_vals.append(pes.cosmo_qhg_zero[i][0])
        if pes.boltz:
            e_sum, h_sum, g_sum, qhg_sum, cosmo_qhg_sum = 0.0, 0.0, 0.0, 0.0, 0.0
            sels = []
            for j, e_abs in enumerate(pes.e_abs[i]):
                species = get_species(i, j, options.qh, options.temperature, pes)
                if options.cosmo:
                    species.append(pes.cosmo_qhg_abs[i][j])
                relative = [species[x] - zero_vals[x] for x in range(len(zero_vals))]
                e_sum += np.exp(-relative[1] * AU_TO_J / GAS_CONSTANT / options.temperature)
                h_sum += np.exp(-relative[3] * AU_TO_J / GAS_CONSTANT / options.temperature)
                g_sum += np.exp(-relative[7] * AU_TO_J / GAS_CONSTANT / options.temperature)
                qhg_sum += np.exp(-relative[8] * AU_TO_J / GAS_CONSTANT / options.temperature)
                cosmo_qhg_sum += np.exp(-relative[9] * AU_TO_J / GAS_CONSTANT / options.temperature)

        print_select_pes_output(options, path, pes, delim_row)

        pes_options(e_sum, g_sum, h_sum, i, options, pes, qhg_sum, sels, delim_row, zero_vals)


def get_species(i, j, qh_flag, temperature, pes):
    if qh_flag:
        species = [pes.spc_abs[i][j], pes.e_abs[i][j], pes.zpe_abs[i][j], pes.h_abs[i][j],
                   pes.qh_abs[i][j], temperature * pes.s_abs[i][j],
                   temperature * pes.qs_abs[i][j], pes.g_abs[i][j], pes.qhg_abs[i][j]]
    else:
        species = [pes.spc_abs[i][j], pes.e_abs[i][j], pes.zpe_abs[i][j], pes.h_abs[i][j],
                   temperature * pes.s_abs[i][j],
                   temperature * pes.qs_abs[i][j],
                   pes.g_abs[i][j], pes.qhg_abs[i][j]]
    return species


def print_select_pes_output(options, path, pes, delim_row):
    if options.spc:
        print("\n   " + '{:<40}'.format("RXN: " + path + " (" + pes.units + ") ", ))
        if options.qh and options.cosmo:
            print('{:>13} {:>13} {:>10} {:>13} {:>13} {:>10} {:>10} {:>14} {:>14} {:>14}'.
                  format(" DE_SPC", "DE", "DZPE", "DH_SPC", "qh-DH_SPC", "T.DS", "T.qh-DS",
                         "DG(T)_SPC", "qh-DG(T)_SPC", 'COSMO-qh-G(T)_SPC'))
        elif options.qh:
            print('{:>13} {:>13} {:>10} {:>13} {:>13} {:>10} {:>10} {:>14} {:>14}'.
                  format(" DE_SPC", "DE", "DZPE", "DH_SPC", "qh-DH_SPC", "T.DS", "T.qh-DS",
                         "DG(T)_SPC", "qh-DG(T)_SPC"))
        elif options.cosmo:
            print('{:>13} {:>13} {:>10} {:>13} {:>13} {:>10} {:>10} {:>14} {:>14}'.
                  format(" DE_SPC", "DE", "DZPE", "DH_SPC", "T.DS", "T.qh-DS", "DG(T)_SPC",
                         "qh-DG(T)_SPC", 'COSMO-qh-G(T)_SPC'))
        else:
            print('{:>13} {:>13} {:>10} {:>13} {:>10} {:>10} {:>14} {:>14}'.
                  format(" DE_SPC", "DE", "DZPE", "DH_SPC", "T.DS", "T.qh-DS", "DG(T)_SPC",
                         "qh-DG(T)_SPC"))
    else:
        print("\n   " + '{:<40}'.format("RXN: " + path + " (" + pes.units + ") ", ))
        if options.qh and options.cosmo:
            print('{:>13} {:>10} {:>13} {:>13} {:>10} {:>10} {:>13} {:>13} {:>13}'.
                  format(" DE", "DZPE", "DH", "qh-DH", "T.DS", "T.qh-DS", "DG(T)", "qh-DG(T)",
                         'COSMO-qh-G(T)'))
        elif options.qh:
            print('{:>13} {:>10} {:>13} {:>13} {:>10} {:>10} {:>13} {:>13}'.
                  format(" DE", "DZPE", "DH", "qh-DH", "T.DS", "T.qh-DS", "DG(T)", "qh-DG(T)"))
        elif options.cosmo:
            print('{:>13} {:>10} {:>13} {:>13} {:>10} {:>10} {:>13} {:>13}'.
                  format(" DE", "DZPE", "DH", "T.DS", "T.qh-DS", "DG(T)", "qh-DG(T)", 'COSMO-qh-G(T)'))
        else:
            print('{:>13} {:>10} {:>13} {:>10} {:>10} {:>13} {:>13}'.
                  format(" DE", "DZPE", "DH", "T.DS", "T.qh-DS", "DG(T)", "qh-DG(T)"))
    print("\n" + delim_row)


def calc_enantio_excess(clustering, clusters, dup_list, files, options, thermo_data):
    delim_row = "   " + '-' * 109
    boltz_facts, weighted_free_energy, boltz_sum = get_boltz(files, thermo_data, clustering, clusters,
                                                             options.temperature, dup_list)
    ee, er, ratio, dd_free_energy, failed, preference = get_selectivity(options.ee, files, boltz_facts,
                                                                        boltz_sum,
                                                                        options.temperature, dup_list)
    if not failed:
        print("\n   " + '{:<39} {:>13} {:>13} {:>13} {:>13} {:>13}'.format("Selectivity", "Excess (%)",
                                                                           "Ratio (%)", "Ratio", "Major Iso",
                                                                           "ddG"))
        print("\n" + delim_row)
        print('\no {:<40} {:13.2f} {:>13} {:>13} {:>13} {:13.2f}'.format('', ee, er, ratio, preference,
                                                                         dd_free_energy))
        print("\n" + delim_row + "\n")


def get_boltz(files, thermo_data, clustering, clusters, temperature, dup_list):
    # Obtain Boltzmann factors, Boltzmann sums, and weighted free energy values
    # used for --ee and --boltz options
    boltz_facs, weighted_free_energy, e_rel, e_min, boltz_sum = {}, {}, {}, sys.float_info.max, 0.0

    for file in files:  # Need the most stable structure
        bbe = thermo_data[file]
        if hasattr(bbe, "qh_gibbs_free_energy"):
            if bbe.qh_gibbs_free_energy is not None:
                if bbe.qh_gibbs_free_energy < e_min:
                    e_min = bbe.qh_gibbs_free_energy

    if clustering:
        for n, cluster in enumerate(clusters):
            boltz_facs['cluster-' + ALPHABET[n].upper()] = 0.0
            weighted_free_energy['cluster-' + ALPHABET[n].upper()] = 0.0
    # Calculate E_rel and Boltzmann factors
    for file in files:
        duplicate = False
        if len(dup_list) != 0:
            for dup in dup_list:
                if dup[0] == file:
                    duplicate = True
        if not duplicate:
            bbe = thermo_data[file]
            if hasattr(bbe, "qh_gibbs_free_energy"):
                if bbe.qh_gibbs_free_energy is not None:
                    # noinspection PyTypeChecker
                    e_rel[file] = bbe.qh_gibbs_free_energy - e_min
                    boltz_facs[file] = np.exp(-e_rel[file] * AU_TO_J / GAS_CONSTANT / temperature)
                    if clustering:
                        for n, cluster in enumerate(clusters):
                            for structure in cluster:
                                if structure == file:
                                    boltz_facs['cluster-' + ALPHABET[n].upper()] += np.exp(
                                        -e_rel[file] * AU_TO_J / GAS_CONSTANT / temperature)
                                    weighted_free_energy['cluster-' + ALPHABET[n].upper()] += np.exp(
                                        -e_rel[file] * AU_TO_J / GAS_CONSTANT / temperature) * bbe.qh_gibbs_free_energy
                    boltz_sum += np.exp(-e_rel[file] * AU_TO_J / GAS_CONSTANT / temperature)

    return boltz_facs, weighted_free_energy, boltz_sum


def get_selectivity(pattern, files, boltz_facs, boltz_sum, temperature, dup_list):
    # Calculate selectivity - enantioselectivity/diastereomeric ratio
    # based on boltzmann factors of given stereoisomers
    # Grab files for selectivity calcs
    a_files, b_files, a_sum, b_sum, failed, pref = [], [], 0.0, 0.0, False, ''
    pattern = pattern.split(',')
    a_item = ''.join(a for a in pattern[0] if a.isalnum())
    b_item = ''.join(b for b in pattern[1] if b.isalnum())
    a_files.extend(glob(pattern[0]))
    b_files.extend(glob(pattern[1]))

    if len(a_files) == 0 or len(b_files) == 0:
        warning("\n   Filenames have not been formatted correctly for determining selectivity,\n   Make sure the "
                "filename contains either {} or {}".format(a_item, b_item))
        raise InvalidDataError("   Please edit either your filenames or selectivity pattern argument and try again\n")
    # Grab Boltzmann sums
    for file in files:
        duplicate = False
        if len(dup_list) != 0:
            for dup in dup_list:
                if dup[0] == file:
                    duplicate = True
        if not duplicate:
            if file in a_files:
                a_sum += boltz_facs[file] / boltz_sum
            elif file in b_files:
                b_sum += boltz_facs[file] / boltz_sum
    # Get ratios
    a_round = round(a_sum * 100)
    b_round = round(b_sum * 100)
    r = str(a_round) + ':' + str(b_round)
    if a_sum > b_sum:
        pref = a_item
        try:
            ratio = a_sum / b_sum
            if ratio < 3:
                ratio = str(round(ratio, 1)) + ':1'
            else:
                ratio = str(round(ratio)) + ':1'
        except ZeroDivisionError:
            ratio = '1:0'
    else:
        pref = b_item
        try:
            ratio = b_sum / a_sum
            if ratio < 3:
                ratio = '1:' + str(round(ratio, 1))
            else:
                ratio = '1:' + str(round(ratio))
        except ZeroDivisionError:
            ratio = '0:1'
    ee = (a_sum - b_sum) * 100.
    if ee == 0:
        warning("\n   No files found for an enantioselectivity analysis, adjust the stereodetermining step "
                "name and try again.\n")
        failed = True
    ee = abs(ee)
    if ee > 99.99:
        ee = 99.99
    try:
        dd_free_energy = GAS_CONSTANT / AU_TO_J * temperature * np.log((50 + abs(ee) / 2.0) / (50 - abs(ee) / 2.0)) * \
                         EHPART_TO_KCAL_MOL
    except ZeroDivisionError:
        dd_free_energy = 0.0
    return ee, r, ratio, dd_free_energy, failed, pref


def output_pes_temp_interval(options, delim_row, interval, interval_bbe_data, interval_thermo_data, file_list):
    e_sum, h_sum, g_sum, qhg_sum, sels = None, None, None, None, None  # make IDE happy
    # Interval applied to PES
    delim_row = delim_row + '-' * 22
    for i in range(len(interval)):
        bbe_vals = []
        for j in range(len(interval_bbe_data)):
            bbe_vals.append(interval_bbe_data[j][i])
        interval_thermo_data.append(dict(zip(file_list, bbe_vals)))
    j = 0
    for i in interval:
        temp = float(i)
        if options.cosmo_int:
            pes = GetPES(options.pes, interval_thermo_data[j], temp, options.gconf, options.qh,
                         cosmo=True)
        else:
            pes = GetPES(options.pes, interval_thermo_data[j], temp, options.gconf, options.qh)

        for k, path in enumerate(pes.path):
            if options.qh:
                zero_vals = [pes.spc_zero[k][0], pes.e_zero[k][0], pes.zpe_zero[k][0], pes.h_zero[k][0],
                             pes.qh_zero[k][0], temp * pes.ts_zero[k][0], temp * pes.qhts_zero[k][0],
                             pes.g_zero[k][0], pes.qhg_zero[k][0]]
            else:
                zero_vals = [pes.spc_zero[k][0], pes.e_zero[k][0], pes.zpe_zero[k][0], pes.h_zero[k][0],
                             temp * pes.ts_zero[k][0], temp * pes.qhts_zero[k][0], pes.g_zero[k][0],
                             pes.qhg_zero[k][0]]
            if options.cosmo_int:
                zero_vals.append(pes.cosmo_qhg_abs[k][0])
            if pes.boltz:
                e_sum, h_sum, g_sum, qhg_sum = 0.0, 0.0, 0.0, 0.0
                sels = []
                for n, e_abs in enumerate(pes.e_abs[k]):
                    species = get_species(k, n, options.qh, options.temperature, pes)
                    relative = [species[x] - zero_vals[x] for x in range(len(zero_vals))]
                    e_sum += np.exp(-relative[1] * AU_TO_J / GAS_CONSTANT / temp)
                    h_sum += np.exp(-relative[3] * AU_TO_J / GAS_CONSTANT / temp)
                    g_sum += np.exp(-relative[7] * AU_TO_J / GAS_CONSTANT / temp)
                    qhg_sum += np.exp(-relative[8] * AU_TO_J / GAS_CONSTANT / temp)
            if options.spc:
                print("\n   " + '{:<40}'.format("RXN: " + path + " (" + pes.units + ")  at T: " +
                                                str(temp)))
                if options.qh and options.cosmo_int:
                    print('{:>13} {:>13} {:>10} {:>13} {:>13} {:>10} {:>10} {:>14} {:>14} {:>14}'.
                          format(" DE_SPC", "DE", "DZPE", "DH_SPC", "qh-DH_SPC", "T.DS", "T.qh-DS",
                                 "DG(T)_SPC", "qh-DG(T)_SPC", 'COSMO-qh-G(T)_SPC'))
                elif options.qh:
                    print('{:>13} {:>13} {:>10} {:>13} {:>13} {:>10} {:>10} {:>14} '
                          '{:>14}'.format(" DE_SPC", "DE", "DZPE", "DH_SPC", "qh-DH_SPC", "T.DS",
                                          "T.qh-DS", "DG(T)_SPC", "qh-DG(T)_SPC"))
                elif options.cosmo_int:
                    print('{:>13} {:>13} {:>10} {:>13} {:>13} {:>10} {:>10} {:>14} '
                          '{:>14}'.format(" DE_SPC", "DE", "DZPE", "DH_SPC", "T.DS", "T.qh-DS",
                                          "DG(T)_SPC", "qh-DG(T)_SPC", 'COSMO-qh-G(T)_SPC'))
                else:
                    print('{:>13} {:>13} {:>10} {:>13} {:>10} {:>10} {:>14} '
                          '{:>14}'.format(" DE_SPC", "DE", "DZPE", "DH_SPC", "T.DS", "T.qh-DS",
                                          "DG(T)_SPC", "qh-DG(T)_SPC"))
            else:
                print("\n   " + '{:<40}'.format("RXN: " + path + " (" + pes.units + ")  at T: " + str(temp)))
                if options.qh and options.cosmo_int:
                    print('{:>13} {:>10} {:>13} {:>13} {:>10} {:>10} {:>13} {:>13} '
                          '{:>13}'.format(" DE", "DZPE", "DH", "qh-DH", "T.DS", "T.qh-DS", "DG(T)",
                                          "qh-DG(T)", 'COSMO-qh-G(T)'))
                elif options.qh:
                    print('{:>13} {:>10} {:>13} {:>13} {:>10} {:>10} {:>13} '
                          '{:>13}'.format(" DE", "DZPE", "DH", "qh-DH", "T.DS", "T.qh-DS", "DG(T)",
                                          "qh-DG(T)"))
                elif options.cosmo_int:
                    print('{:>13} {:>10} {:>13} {:>13} {:>10} {:>10} {:>13} '
                          '{:>13}'.format(" DE", "DZPE", "DH", "T.DS", "T.qh-DS", "DG(T)", "qh-DG(T)",
                                          'COSMO-qh-G(T)'))
                else:
                    print('{:>13} {:>10} {:>13} {:>10} {:>10} {:>13} '
                          '{:>13}'.format(" DE", "DZPE", "DH", "T.DS", "T.qh-DS", "DG(T)", "qh-DG(T)"))
            print(delim_row)

            for n, e_abs in enumerate(pes.e_abs[k]):
                if options.qh:
                    species = [pes.spc_abs[k][n], pes.e_abs[k][n], pes.zpe_abs[k][n], pes.h_abs[k][n],
                               pes.qh_abs[k][n], temp * pes.s_abs[k][n], temp * pes.qs_abs[k][n],
                               pes.g_abs[k][n], pes.qhg_abs[k][n]]
                else:
                    species = [pes.spc_abs[k][n], pes.e_abs[k][n], pes.zpe_abs[k][n], pes.h_abs[k][n],
                               temp * pes.s_abs[k][n], temp * pes.qs_abs[k][n], pes.g_abs[k][n],
                               pes.qhg_abs[k][n]]
                if options.cosmo_int:
                    species.append(pes.cosmo_qhg_abs[k][n])
                relative = [species[x] - zero_vals[x] for x in range(len(zero_vals))]
                if pes.units == 'kJ/mol':
                    formatted_list = [AU_TO_J / 1000.0 * x for x in relative]
                else:
                    formatted_list = [EHPART_TO_KCAL_MOL * x for x in relative]  # Defaults to kcal/mol

                if options.spc:
                    if options.qh and options.cosmo_int:
                        if pes.dec == 1:
                            print('o  {:<39} {:13.2f} {:13.1f} {:10.1f} {:13.1f} {:13.1f} {:10.1f} '
                                  '{:10.1f} {:13.1f} {:13.1f} {:13.1f}'.format(pes.species[k][n],
                                                                               *formatted_list))
                        if pes.dec == 2:
                            print('o  {:<39} {:13.2f} {:13.2f} {:10.2f} {:13.2f} {:13.2f} {:10.2f} '
                                  '{:10.2f} {:13.2f} {:13.2f} {:13.2f}'.format(pes.species[k][n],
                                                                               *formatted_list))
                    elif options.qh or options.cosmo_int:
                        if pes.dec == 1:
                            print('o  {:<39} {:13.2f} {:13.1f} {:10.1f} {:13.1f} {:13.1f} {:10.1f} '
                                  '{:10.1f} {:13.1f} {:13.1f}'.format(pes.species[k][n],
                                                                      *formatted_list))
                        if pes.dec == 2:
                            print('o  {:<39} {:13.2f} {:13.2f} {:10.2f} {:13.2f} {:13.2f} {:10.2f} '
                                  '{:10.2f} {:13.2f} {:13.2f}'.format(pes.species[k][n],
                                                                      *formatted_list))
                    else:
                        if pes.dec == 1:
                            print('o  {:<39} {:13.2f} {:13.1f} {:10.1f} {:13.1f} {:10.1f} {:10.1f} '
                                  '{:13.1f} {:13.1f}'.format(pes.species[k][n], *formatted_list))
                        if pes.dec == 2:
                            print('o  {:<39} {:13.2f} {:13.2f} {:10.2f} {:13.2f} {:10.2f} {:10.2f} '
                                  '{:13.2f} {:13.2f}'.format(pes.species[k][n], *formatted_list))
                else:
                    formatted_list = formatted_list[1:]
                    format_1 = 'o  {:<39} {:13.2f} {:10.1f} {:13.1f} {:13.1f} {:10.1f} {:10.1f} {:13.1f} ' \
                               '{:13.1f} {:13.1f}'
                    format_2 = 'o  {:<39} {:13.2f} {:10.2f} {:13.2f} {:13.2f} {:10.2f} {:10.2f} {:13.2f} ' \
                               '{:13.2f} {:13.2f}'
                    if options.qh and options.cosmo_int:
                        if pes.dec == 1:
                            print(format_1.format(pes.species[k][n], *formatted_list))
                        if pes.dec == 2:
                            print(format_2.format(pes.species[k][n], *formatted_list))
                    elif options.qh or options.cosmo_int:
                        if pes.dec == 1:
                            print(format_1.format(pes.species[k][n], *formatted_list))
                        if pes.dec == 2:
                            print(format_2.format(pes.species[k][n], *formatted_list))
                    else:
                        if pes.dec == 1:
                            print(format_1.format(pes.species[k][n], *formatted_list))
                        if pes.dec == 2:
                            print(format_2.format(pes.species[k][n], *formatted_list))

                if pes.boltz:
                    boltz = [np.exp(-relative[1] * AU_TO_J / GAS_CONSTANT / options.temperature) / e_sum,
                             np.exp(-relative[3] * AU_TO_J / GAS_CONSTANT / options.temperature) / h_sum,
                             np.exp(-relative[6] * AU_TO_J / GAS_CONSTANT / options.temperature) / g_sum,
                             np.exp(-relative[7] * AU_TO_J / GAS_CONSTANT / options.temperature) / qhg_sum]
                    selectivity = [boltz[x] * 100.0 for x in range(len(boltz))]
                    print("\n  " + 'o  {:<39} {:13.2f}%{:24.2f}%{:35.2f}%{:13.2f}%'.
                          format('', *selectivity))
                    sels.append(selectivity)
            if pes.boltz == 'ee' and len(sels) == 2:
                ee = [sels[0][x] - sels[1][x] for x in range(len(sels[0]))]
                if options.spc:
                    print(delim_row + '\n   {:<39} {:27.2f} {:24.1f} {:35.1f} {:13.1f}'.format('ee (%)', *ee))
                else:
                    print(delim_row + '\n   {:<39} {:13.2f} {:24.1f} {:35.1f} {:13.1f}'.format('ee (%)', *ee))

            print(delim_row)
        j += 1


def output_cosmos_rs_interval(files, options, s_m, l_o_t):
    gsolv_dicts = None  # make IDE happy
    scaling_data_tuple = namedtuple("ScalingData",
                                    ['level_basis', 'zpe_fac', 'zpe_ref', 'zpe_meth', 'harm_fac', 'harm_ref',
                                     'harm_meth', 'fund_fac', 'fund_ref', 'fund_meth'])
    scaling_data_dict, scaling_data_dict_mod = {}, {}
    for row in SCALING_DATA:
        level_basis, zpe_fac, zpe_ref, zpe_meth, harm_fac, harm_ref, harm_meth, fund_fac, fund_ref, fund_meth = row
        data = scaling_data_tuple(*row)
        scaling_data_dict[level_basis.upper()] = scaling_data_dict_mod[level_basis.replace("-", "").upper()] = data

    # Attempt to automatically obtain frequency scale factor,
    # Application of freq scale factors requires all outputs to be same level of theory
    if options.freq_scale_factor:
        if 'ONIOM' not in l_o_t[0]:
            print("\nUser-defined vibrational scale factor of {} for {} level of theory".
                  format(options.freq_scale_factor, l_o_t[0]))
        else:
            print("\nUser-defined vibrational scale factor for QM region is {}".
                  format(options.freq_scale_factor, l_o_t[0]))
    else:
        # Look for vibrational scaling factor automatically
        if all_same(l_o_t):
            level = l_o_t[0].upper()
            for data in (scaling_data_dict_mod, scaling_data_dict_mod):
                if level in data:
                    options.freq_scale_factor = data[level].zpe_fac
                    ref = SCALING_REFS[data[level].zpe_ref]
                    print("\nFound vibrational scaling factor of {:.3f} for {} level of theory\n"
                          "    REF: {}".format(options.freq_scale_factor, l_o_t[0], ref))
                    break
        else:  # Print files and different levels of theory found
            files_l_o_t, levels_l_o_t, filtered_calcs_l_o_t = [], [], []
            for file in files:
                files_l_o_t.append(file)
            for i in l_o_t:
                levels_l_o_t.append(i)
            filtered_calcs_l_o_t.append(files_l_o_t)
            filtered_calcs_l_o_t.append(levels_l_o_t)
            print_check_fails(filtered_calcs_l_o_t[1], filtered_calcs_l_o_t[0], "levels of theory")

    if options.zpe_scale_factor:
        print("    For calculating the ZPE, a scale factor of {} will be used.\n".format(options.zpe_scale_factor))
    else:
        options.zpe_scale_factor = options.freq_scale_factor
        print("    The same scaling factor will be used to calculate the ZPE.\n")

    # Exit program if a comparison of Boltzmann factors is requested and level of theory is not uniform across all files
    if not all_same(l_o_t) and (not options.boltz or not options.ee):
        raise InvalidDataError("When comparing files with Boltzmann factors (with bolts, ee, dr "
                               "options), the level of theory used should be the same for all files.\n ")
    # Exit program if molecular mechanics scaling factor is given and all files are not ONIOM calculations
    if options.mm_freq_scale_factor:
        if all_same(l_o_t) and 'ONIOM' in l_o_t[0]:
            print("\n   User-defined vibrational scale factor " +
                  str(options.mm_freq_scale_factor) + " for MM region of " + l_o_t[0])
            print("    REF: {}".format(ONIOM_SCALE_REF))
        else:
            raise InvalidDataError("\n   Option --vmm is only for use in ONIOM calculation output files.\n    "
                                   "For help use option '-h'\n")

    if not options.freq_scale_factor:
        options.freq_scale_factor = 1.0  # If no scaling factor is found use 1.0
        if all_same(l_o_t):
            print("\n   Using vibrational scale factor {} for {} level of "
                  "theory".format(options.freq_scale_factor, l_o_t[0]))
        else:
            print("\n   Using vibrational scale factor {}: differing levels of theory "
                  "detected.".format(options.freq_scale_factor))

    # Checks to see whether the available free space of a requested solvent is defined
    freespace = get_free_space(options.freespace)
    if freespace != 1000.0:
        print("\n   Specified solvent " + options.freespace + ": free volume " + str(
            "%.3f" % (freespace / 10.0)) + " (mol/l) corrects the translational entropy")

    # Check for implicit solvation
    printed_solv_warn = False
    for i in s_m:
        if ('smd' in i.lower() or 'cpcm' in i.lower()) and not printed_solv_warn:
            print("\n   Caution! Implicit solvation (SMD/CPCM) detected. Enthalpic and entropic terms "
                  "cannot be safely separated. Use them at your own risk!")
            printed_solv_warn = True

    # COSMO-RS temperature interval
    cosmo_solv = None  # make IDE happy
    t_interval = None  # make IDE happy
    if options.cosmo_int:
        args = options.cosmo_int.split(',')
        c_file = args[0]
        c_interval = args[1:]
        print('\n   Reading COSMO-RS file: ' + c_file + ' over a T range of ' + c_interval[0] + '-' +
              c_interval[1] + ' K.')

        t_interval, gsolv_dicts = cosmo_rs_out(c_file, files, interval=c_interval)
        options.temperature_interval = True

    elif options.cosmo:  # Read from COSMO-RS output
        try:
            cosmo_solv = cosmo_rs_out(options.cosmo, files)
            print('\n\n   Reading COSMO-RS file: ' + options.cosmo)
        except ValueError:
            print('\n\n   Warning! COSMO-RS file ' + options.cosmo + ' requested but not found')

    if options.freq_cutoff != 100.0:
        options.S_freq_cutoff = options.freq_cutoff
        options.h_freq_cutoff = options.freq_cutoff

    # Summary of the quasi-harmonic treatment; print out the relevant reference
    print("Entropic quasi-harmonic treatment: frequency cut-off value of {} wavenumbers will be "
          "applied.".format(options.S_freq_cutoff))
    # Only print reference if actually cutting off wavenumbers
    if options.S_freq_cutoff > 0.0:
        if options.qs == "grimme":
            print("    qs = Grimme: Using a mixture of RRHO and Free-rotor vibrational entropies.")
            qs_ref = GRIMME_REF
        elif options.qs == "truhlar":
            print("    qs = Truhlar: Using an RRHO treatment where low frequencies are adjusted to the cut-off value.")
            qs_ref = TRUHLAR_REF
        else:
            raise InvalidDataError("\n   FATAL ERROR: Unknown quasi-harmonic model " + options.qs +
                                   " specified (qs must = grimme or truhlar).")
        print("    REF: " + qs_ref + '\n')
    else:
        # adding one blank line...
        print("")

    # Check if qh-H correction should be applied
    if options.qh:
        print("Enthalpy quasi-harmonic treatment: frequency cut-off value of " + str(
            options.h_freq_cutoff) + " wavenumbers will be applied.")
        print("    qh = Head-Gordon: Using an RRHO treatment with an approximation term for vibrational energy.")
        print("    REF: {}\n".format(HEAD_GORDON_REF))

    # Check if D3 corrections should be applied
    if options.D3:
        print("D3-Dispersion energy with zero-damping will be calculated and included in the energy and "
              "enthalpy terms.\n    REF: {}\n".format(D3_REF))
    if options.D3BJ:
        print("D3-Dispersion energy with Becke-Johnson damping will be calculated and added to the "
              "energy terms.\n    REF: {}\n".format(D3BJ_REF))
    if options.ATM:
        print("The repulsive Axilrod-Teller-Muto 3-body term will be included in the dispersion "
              "correction.\n    REF: {}\n".format(ATM_REF))

    # Check if entropy symmetry correction should be applied
    if options.ssymm:
        print("Ssymm requested. Symmetry contribution to entropy to be calculated using S. Patchkovskii's \n"
              "    open source software 'Brute Force Symmetry Analyzer' available under GNU General Public "
              'License.\n   REF: (C) 1996, 2003 S. Patchkovskii, Serguei.Patchkovskii@sympatico.ca\n'
              '    Atomic radii used to calculate internal symmetry based on Cambridge Structural Database '
              'covalent radii.\n   REF: {}\n'.format(CSD_REF))

    # Whether linked single-point energies are to be used
    if options.spc:
        print("    Link job: combining final single point energy with thermal corrections.\n")

    return cosmo_solv, gsolv_dicts, t_interval


def cosmo_rs_out(dat_file, names, interval=False):
    # Read solvation free energies from a COSMO-RS dat file
    gsolv = {}
    if os.path.exists(dat_file):
        with open(dat_file) as f:
            data = f.readlines()
    else:
        raise ValueError("File {} does not exist".format(dat_file))

    temp = 0
    t_interval = []
    gsolv_dicts = []
    found = False
    old_temp = 0
    gsolv_temp = {}
    if interval:
        for i, line in enumerate(data):
            for name in names:
                if line.find('(' + name.split('.')[0] + ')') > -1 and line.find('Compound') > -1:
                    if data[i - 5].find('Temperature') > -1:
                        temp = data[i - 5].split()[2]
                    temp = float(temp)
                    # noinspection PyUnresolvedReferences
                    if float(interval[0]) < temp < float(interval[1]):
                        if float(temp) not in t_interval:
                            t_interval.append(float(temp))
                        if data[i + 10].find('Gibbs') > -1:
                            gsolv = float(data[i + 10].split()[6].strip()) / EHPART_TO_KCAL_MOL
                            gsolv_temp[name] = gsolv

                            found = True
            if found:
                if old_temp == 0:
                    old_temp = temp
                if temp != old_temp:
                    gsolv_dicts.append(gsolv)  # Store dict at one temp
                    gsolv = {}  # Clear gsolv
                    gsolv.update(gsolv_temp)  # Grab the first one for the new temp
                    old_temp = temp
                gsolv.update(gsolv_temp)
                gsolv_temp = {}
                found = False
        gsolv_dicts.append(gsolv)  # Grab last dict
    else:
        for i, line in enumerate(data):
            for name in names:
                if line.find('(' + name.split('.')[0] + ')') > -1 and line.find('Compound') > -1:
                    if data[i + 11].find('Gibbs') > -1:
                        gsolv = float(data[i + 11].split()[6].strip()) / EHPART_TO_KCAL_MOL
                        # noinspection PyUnresolvedReferences
                        gsolv[name] = gsolv

    if interval:
        return t_interval, gsolv_dicts
    else:
        return gsolv
