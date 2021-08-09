# -*- coding: utf-8 -*-
"""
Gaussian input/output tools
Scripts to set up Gaussian input files and analyze Gaussian output files
"""
from setuptools import setup
import versioneer

DOCLINES = __doc__.split("\n")

requirements = ['progressbar2']

setup(
    # Self-descriptive entries which should always be present
    name='gaussian_wrangler',
    author='Heather B Mayes',
    author_email='hmayes@hmayes.com',
    description=DOCLINES[0],
    long_description="\n".join(DOCLINES[2:]),
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    license='MIT',

    # Which Python importable modules should be included when your package is installed
    packages=['gaussian_wrangler'],

    # Optional include package data to ship with your package
    # Comment out this line to prevent the files from being packaged with your software
    # Extend/modify the list to include/exclude other items as need be
    package_data={'gaussian_wrangler': ["data/*.dat", "hartree/*.*", "good_vibes/*.*"]
                  },

    entry_points={'console_scripts': ['gausscom2pdb = gaussian_wrangler.gausscom2pdb:main',
                                      'gausscom2com = gaussian_wrangler.gausscom2com:main',
                                      'pdbs2gausscoms = gaussian_wrangler.pdbs2gausscoms:main',
                                      'gausslog2com = gaussian_wrangler.gausslog2com:main',
                                      'gausslog2pdb = gaussian_wrangler.gausslog2pdb:main',
                                      'gausslog_unique = gaussian_wrangler.gausslog_unique:main',
                                      'gauss_fragment = gaussian_wrangler.gauss_fragment:main',
                                      'run_gauss = gaussian_wrangler.run_gauss:main',
                                      'check_gauss = gaussian_wrangler.check_gauss:main',
                                      'goodvibes_helper = gaussian_wrangler.goodvibes_helper:main',
                                      'goodvibes_hm = gaussian_wrangler.goodvibes_hm:main',
                                      'plot_steps = gaussian_wrangler.plot_steps:main',
                                      'smi2gausscom = gaussian_wrangler.smi2gausscom:main'
                                      ],
                  },     package_dir={'gaussian_wrangler': 'gaussian_wrangler'},

    test_suite='tests',
    # to install jpype1, use conda: install -c conda-forge jpype1
    # install_requires=requirements,
    # install_requires=[],   # Required packages, pulls from pip if needed; do not use for Conda deployment
    install_requires=['numpy', 'matplotlib', 'scipy', 'common-wrangler>=0.3.7.2', 'jpype1', 'pubchempy'],
    tests_require=['pytest'],
    # Additional entries you may want simply uncomment the lines you want and fill in the data
    # url='http://www.my_package.com',  # Website
    platforms=['Linux',
               'Mac OS-X',
               'Unix',
               'Windows'],            # Valid platforms your code works on, adjust to your flavor
    python_requires=">=3.0",  # Python version restrictions
    # Manual control if final package is compressible or not, set False to prevent the .egg from being made
    # zip_safe=False,

)
