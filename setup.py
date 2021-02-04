#!/bin/env python

#######################################################################
#  Copyright (C) 2020 Vinh Tran
#
#  fdog is the python package for feature-aware directed ortholog
#  search. fdog is a free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  hamstr1s is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with fcat.  If not, see <http://www.gnu.org/licenses/>.
#
#######################################################################

from setuptools import setup, find_packages

with open("README.md", "r") as input:
    long_description = input.read()

setup(
    name="fcat",
    version="0.0.26",
    python_requires='>=3.7.0',
    description="Python package for fCAT, a feature-aware completeness assessment tool",
    long_description=long_description,
    long_description_content_type="text/markdown",
    author="Vinh Tran",
    author_email="tran@bio.uni-frankfurt.de",
    url="https://github.com/BIONF/fCAT",
    packages=find_packages(),
    package_data={'': ['*']},
    install_requires=[
        'biopython',
        'tqdm',
        'ete3',
        'six',
        'greedyFAS>=1.4.0',
        'fdog>=0.0.8',
        'rpy2',
        'tzlocal',
        'scipy'
    ],
    entry_points={
        'console_scripts': ["fcat = fcat.fcat:main",
                            "fcat.cutoff = fcat.calcCutoff:main",
                            "fcat.ortho = fcat.searchOrtho:main",
                            "fcat.report = fcat.assessCompleteness:main",
                            "fcat.createProfile = fcat.createPhyloprofile:main",
                            "fcat.mergeOutput = fcat.mergePhyloprofile:main",
                            "fcat.showTaxa = fcat.showTaxa:main"],
    },
    license="GPL-3.0",
    classifiers=[
        "Environment :: Console",
        "Intended Audience :: End Users/Desktop",
        "License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)",
        "Natural Language :: English",
        "Programming Language :: Python :: 3",
    ],
)
