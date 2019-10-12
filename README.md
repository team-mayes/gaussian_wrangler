nrel_tools
==============================
[//]: # (Badges)
[![Travis Build Status](https://travis-ci.org/REPLACE_WITH_OWNER_ACCOUNT/nrel_tools.png)](https://travis-ci.org/REPLACE_WITH_OWNER_ACCOUNT/nrel_tools)
[![AppVeyor Build status](https://ci.appveyor.com/api/projects/status/REPLACE_WITH_APPVEYOR_LINK/branch/master?svg=true)](https://ci.appveyor.com/project/REPLACE_WITH_OWNER_ACCOUNT/nrel_tools/branch/master)
[![codecov](https://codecov.io/gh/REPLACE_WITH_OWNER_ACCOUNT/nrel_tools/branch/master/graph/badge.svg)](https://codecov.io/gh/REPLACE_WITH_OWNER_ACCOUNT/nrel_tools/branch/master)

A suite of scripts that have been helpful primarily with work flows involving Gaussian.

To install, obtain the tarball (from https://www.dropbox.com/sh/spiu46a0mrgtean/AADfeFWQsJNUpD2UkYNe6gCoa?dl=0)
or by creating it from the project using `python setup.py sdist`). Then install it on your machine, e.g.:

`pip install --upgrade nrel_tools-0.0.0.tar.gz --user`
    
You will also need to install common_wrangler. You can copy the tarball from 
https://www.dropbox.com/sh/spiu46a0mrgtean/AADfeFWQsJNUpD2UkYNe6gCoa?dl=0 or download the project 
and built it yourself (available at https://github.com/team-mayes/common_wrangler). You can use a similar command to install it:

`pip install --upgrade nrel_tools-0.0.0.tar.gz --user`

**check_gauss**: Checks for normal termination of Gaussian output files in a specified directory, and moves them to a new location.
You can specify the directory where to look, where to move them two, and the extension name of the output files.

**gauss_fragment**: Given either a Gaussian input or output file, and a list of pairs of atoms, it will produce files 
to run a counterpoise correction calculation and optimize each fragment. Currently, the script assumes 
the initial molecule (or molecules; see below) is neutral with singlet multiplicity. If you would like to use the script 
for other cases, please contact the developer.
 
By default, the script will assume that the pair of atoms are (close enough to be) bonded, are not both part of a ring, and there is one 
molecule in the file. Then, the fragments will be radicals (with determined by the type of bond being broken). 

If the `` option is used, the script will assume that there are two molecules, and the pair of atoms provided includes
one atom from each of the two molecules. 


### Copyright

Copyright (c) 2019, Heather B Mayes


#### Acknowledgements
 
Project based on the 
[Computational Chemistry Python Cookiecutter](https://github.com/choderalab/cookiecutter-python-comp-chem)
