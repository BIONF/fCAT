# fCAT
[![PyPI version](https://badge.fury.io/py/fcat.svg)](https://pypi.org/project/fcat/)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Build Status](https://travis-ci.com/BIONF/fCAT.svg?branch=main)](https://travis-ci.com/BIONF/fCAT)

Python package for fCAT, a **f**eature-aware **C**ompleteness **A**ssessment **T**ool

# Table of Contents
* [How to install](#how-to-install)
* [Usage](#usage)
* [Output](#output)
* [Bugs](#bugs)
* [Contributors](#contributors)
* [Contact](#contact)

# How to install

*fCAT* tool is distributed as a python package called *fcat*. It is compatible with [Python ≥ v3.7](https://www.python.org/downloads/).

You can install *fcat* using `pip`:
```
python3 -m pip install fcat
```

or, in case you do not have admin rights, and don't use package systems like Anaconda to manage environments you need to use the `--user` option:
```
python3 -m pip install --user fcat
```

and then add the following line to the end of your **~/.bashrc** or **~/.bash_profile** file, restart the current terminal to apply the change (or type `source ~/.bashrc`):

```
export PATH=$HOME/.local/bin:$PATH
```

# Usage

The complete process of *fCAT* can be done using one function `fcat`
```
fcat --coreDir /path/to/fcat_data --coreSet eukaryota --refspecList "HOMSA@9606@2" --querySpecies /path/to/query.fa [--annoQuery /path/to/query.json] [--outDir /path/to/fcat/output]
```

where **eukaryota** is name of the fCAT core set (equivalent to BUSCO set); **HOMSA@9606@2** is the reference species from that core set that will be used for the ortholog search; **query** is the name of species of interest. If `--annoQuery` not specified, *fCAT* fill do the feature annotation for the query proteins using [FAS tool](https://github.com/BIONF/FAS).

# Output
You will find the output in the */path/to/fcat/output/fcatOutput/eukaryota/query/* folder, where */path/to/fcat/output/* could be your current directory if you not specified `--outDir` when running `fcat`. The following important output files/folders can be found:

- *all_summary.txt*: summary of the completeness assessment using all 4 score modes
- *fdogOutput.tar.gz*: a zipped file of the ortholog search result
- *mode_1*, *mode_2*, *mode_3* and *mode_4*: detailed output for each score mode
- *phyloprofileOutput*: folder contains output phylogenetic profile data that can be used with [PhyloProfile tool](https://github.com/BIONF/PhyloProfile)

# Bugs
Any bug reports or comments, suggestions are highly appreciated. Please [open an issue on GitHub](https://github.com/BIONF/fCAT/issues/new) or be in touch via email.

# Contributors
- [Vinh Tran](https://github.com/trvinh)
- [Giang Nguyen](https://github.com/giangnguyen0709)

# Contact
For further support or bug reports please contact: tran@bio.uni-frankfurt.de
