# -*- coding: utf-8 -*-
"""
nrel_tools
Find low energy conformers using RDKit
"""
from setuptools import setup
# import versioneer

DOCLINES = __doc__.split("\n")

requirements = ['progressbar2']

setup(
    # Self-descriptive entries which should always be present
    name='nrel_tools',
    author='Heather B Mayes',
    description=DOCLINES[0],
    long_description="\n".join(DOCLINES[2:]),
    # version=versioneer.get_version(),
    # cmdclass=versioneer.get_cmdclass(),
    license='MIT',

    # Which Python importable modules should be included when your package is installed
    packages=['nrel_tools'],

    # Optional include package data to ship with your package
    # Comment out this line to prevent the files from being packaged with your software
    # Extend/modify the list to include/exclude other items as need be
    package_data={'nrel_tools': ["data/*.dat"]
                  },

    entry_points={'console_scripts': ['gausscom2pdb = nrel_tools.gausscom2pdb:main',
                                      'pdbs2gausscoms = nrel_tools.pdbs2gausscoms:main',
                                      'gausslog2pdb = nrel_tools.gausslog2pdb:main',
                                      'gausslog2com = nrel_tools.gausslog2com:main',
                                      'gauss_fragment = nrel_tools.gauss_fragment:main',
                                      'run_gauss = nrel_tools.run_gauss:main',
                                      'check_gauss = nrel_tools.check_gauss:main',
                                      'goodvibes_helper = nrel_tools.goodvibes_helper:main',
                                      'plot_delta_g = nrel_tools.plot_delta_g:main'
                                      ],
                  },     package_dir={'nrel_tools': 'nrel_tools'},

    test_suite='tests', install_requires=['numpy', 'six', 'matplotlib']
    # Additional entries you may want simply uncomment the lines you want and fill in the data
    # author_email='me@place.org',      # Author email
    # url='http://www.my_package.com',  # Website
    # install_requires=[],              # Required packages, pulls from pip if needed; do not use for Conda deployment
    # platforms=['Linux',
    #            'Mac OS-X',
    #            'Unix',
    #            'Windows'],            # Valid platforms your code works on, adjust to your flavor
    # python_requires=">=3.5" # Python version restrictions

    # Manual control if final package is compressible or not, set False to prevent the .egg from being made
    # zip_safe=False,

)
