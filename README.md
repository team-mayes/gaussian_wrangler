gaussian_wrangler
==============================
[//]: # (Badges)
![PythonCI](https://github.com/team-mayes/gaussian_wrangler/workflows/PythonCI/badge.svg)

A suite of scripts that have been helpful primarily with work flows involving Gaussian.

## Installation

In <= 3 simple steps!

+ If you do not already have python \geq 3.0 installed, start with that installation. 
[Anaconda](https://www.anaconda.com/distribution/) or [miniconda](https://docs.conda.io/en/latest/miniconda.html) are 
recommended.
+ Install [RDKit](https://www.rdkit.org/docs/Install.html) using conda. 
  -  To install it in a new environment, run: `conda create -c rdkit -n gaussian-wrangler rdkit`, followed by 
     `conda activate lignin-wrangler` 
  -  If you already have a conda environment created and want to add rdkit to that environment, instead run: 
     `conda install -c conda-forge rdkit`  
+ Install this package with pip (`pip install gaussian-wrangler`) or by building it yourself (see below)

## Build

If you prefer to build the package yourself:

+ Clone this package
+ If desired, run the test on your platform with `pytest` (within the `gaussian_wrangler` main directory)
+ Build the tarball via `python setup.py sdist`
+ Install the resulting tarball with `pip install dist/gaussian*tar.gz`

## Included Scripts

Below is a brief description of the scripts included in this package. More detail can be obtained by running the 
script with the `-h` option, which also provides an overview of the scripts, as well as the options available.

**check_gauss**: There are two main functions:
1) Checks for normal termination of Gaussian output files in a specified directory, and moves them to a new location.
You can specify the directory where to look, where to move them two, and the extension name of the output files.
2) Checks for convergence of Gaussian output files: either only the final convergence (`-z` option) or for each step 
(`-s` option).

**gauss_fragment**: Given either a Gaussian input or output file, and a list of pairs of atoms, it will produce files 
to run a counterpoise correction calculation and optimize each fragment. Currently, the script assumes 
the initial molecule (or molecules; see below) is neutral with singlet multiplicity. If you would like to use the script 
for other cases, please contact the developer.
 
By default, the script will assume that the pair of atoms are (close enough to be) bonded, are not both part of a ring, 
and there is one molecule in the file. Then, the fragments will be radicals (with determined by the type of bond being 
broken). 

If the `two_molecules` option is used, the script will assume that there are two molecules, and the pair of atoms 
provided includes
one atom from each of the two molecules. 

**gausscom2com**: This script combines the atomic coordinates from one file, with all 
Gaussian input specifications from another file.

**gausscom2pdb**: As you might expect, this script takes the atoms and coordinates from a Gaussian input file and 
creates a PDB from them. If provided a template PDB file, it will replace the coordinates in that PDB with those
from the Gaussian input file. Otherwise, it will create a generic one. 

**gausslog2com**: As you might expect, this script takes the atoms and coordinates from a Gaussian output file and 
generates a Gaussian input file using the Gaussian input specifications from a template file. By default, the last 
set of coordinates from the output file is used. The "-e" option instead takes the coordinates from the lowest-energy 
step in the log file.

**gausslog_unique**: This script compares results from output files to determine if conformations are identical (as 
often happens when optimizing a large number of postulated conformers) by checking if the differences in dihedral 
angles are smaller than a user-specified threshold (the default is 5 degrees). Currently, the program assumes that 
the atom order is the same in each file; if they are not, they will be considered different conformations. It also does
not check if symmetry makes two conformers identical--it only checks dihedral angles based on the atom numbers. 
It outputs a list of output files with only one file name per unique conformation. When multiple output files reflect 
the same conformation, the output file with the best convergence is shown. The program will also print a warning if 
the final calculations is not fully converged (as can happen due to optimization calculations using estimated 
convergence; a subsequent frequency calculation can then show that the structure is not fully optimized).

**goodvibes_helper**: This script uses goodvibes to calculate entropy and enthalpy at multiple temperatures. It is 
set up to help calculate the enthalpy of reaction at specified temperatures, and/or to determine the kinetic parameters 
(A, E_A) for reactions. Note: this program requires requires [hartree](https://github.com/team-mayes/hartree), and
will look for [this version of goodvibes](https://github.com/team-mayes/GoodVibes), which has a correction for 
calculations in the condensed phase at multiple temperatures.

**pdbs2gausscoms**: This script combines the coordinates from a PDB file (which may have multiple PDB entries) with 
the Gaussian input specifications from a template file to generate Gaussian input files. 

**plot_steps**: This script makes enthalpy and/or free energy diagrams, given a list of values.

**run_gauss**: This script prepares inputs and runs Gaussian jobs, or (with `-s` and `-l` options) prepares and submits
slurm jobs to run Gaussian. It can be used to run (and/or submit) a series of Gaussian jobs.

The script includes functions to determine inputs to Gaussian based on the specifications of the machine running the 
job. For example, if the '{default_route}' parameter is included in the 'job_run_tpl', the script will 
determine a MaxDisk and Cachesize based on the running systems's specifications: specifically, 90% of free diskspace 
and Gaussian's cache size algorithm of ('cache size' * 1024)/num_procs. The calculations assume the whole node will be 
used for the Gaussian job.  These calculations are used to create a Default.Route file for Gaussian to read. This file 
will be created in the current directory, or in a directory specified with the 'scratch_dir' parameter in the 
configuration file.

The program will also determine an appropriate amount of memory to allocate and max number of cores to use (also 
assuming the whole node is used for the Gaussian job) if a `{mem}` and/or `{proc_list}` parameter is included in the 
`job_run_tpl` and not specified in the configuration file.

### Copyright

Copyright (c) 2021, Heather B Mayes


#### Acknowledgements
 
Project based on the 
[Computational Chemistry Python Cookiecutter](https://github.com/choderalab/cookiecutter-python-comp-chem)
