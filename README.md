# fCATpy
[![PyPI version](https://badge.fury.io/py/fcatpy.svg)](https://pypi.org/project/fcatpy/)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Build Status](https://travis-ci.com/BIONF/fCATpy.svg?branch=master)](https://travis-ci.com/BIONF/fCATpy)

Python package for fCAT, a feature-aware completeness assessment tool

# Table of Contents
* [How to install](#how-to-install)
* [Usage](#usage)
* [Bugs](#bugs)
* [Contributors](#contributors)
* [Contact](#contact)

# How to install

*fCAT* tool is distributed as a python package called *fcatpy*. It is compatible with [Python ≥ v3.7](https://www.python.org/downloads/).

You can install *fcatpy* using `pip`:
```
# python3 -m pip install fcatpy
python3 -m pip install git+https://github.com/BIONF/fCATpy
```

or, in case you do not have admin rights, and don't use package systems like Anaconda to manage environments you need to use the `--user` option:
```
# python3 -m pip install --user fcatpy
python3 -m pip install --user git+https://github.com/BIONF/fCATpy
```

and then add the following line to the end of your **~/.bashrc** or **~/.bash_profile** file, restart the current terminal to apply the change (or type `source ~/.bashrc`):

```
export PATH=$HOME/.local/bin:$PATH
```

# Usage

*fCAT* algorithm consists of 3 main steps:

1) Calculate group-specific cutoffs for a core set
```
fcat.cutoff --coreDir /path/to/core/sets --coreSet eukaryota --annoDir /path/to/core/weight_dir --blastDir /path/to/core/blast_dir --cpus 4
```

2) Search for orthologs in a gene set of interst and create phylogenetic profiles
```
fcat.ortho --coreDir /path/to/core/sets --coreSet eukaryota --annoDir /path/to/core/weight_dir --blastDir /path/to/core/blast_dir --refspecList "HOMSA@9606@2" --querySpecies /path/to/query.fa --annoQuery /path/to/query.json --cpus 4 --cleanup
```

3) Create report for completeness assessment
```
fcat.report --coreDir /path/to/core/sets --coreSet eukaryota --outDir /path/to/fcat/output --queryID queryID --mode 1
```

The complete process can be done using one function `fcat`
```
fcat --coreDir /path/to/core/sets --coreSet eukaryota --annoDir /path/to/core/weight_dir --blastDir /path/to/core/blast_dir --refspecList "HOMSA@9606@2" --querySpecies /path/to/query.fa --annoQuery /path/to/query.json --cpus 4 --cleanup
```

*NOTE: currently there is an [issue with rpy2 library](https://github.com/rpy2/rpy2/issues/739), in which step 1 cannot be run using `fcat.cutoff` function. You must instead run the script directly using the python command e.g. `python3 /path/to/fcatpy/calcCutoff.py -h`. The function `fcat` is currently cannot be used for the same reason!*

# Bugs
Any bug reports or comments, suggestions are highly appreciated. Please [open an issue on GitHub](https://github.com/BIONF/fCATpy/issues/new) or be in touch via email.

# Contributors
- [Vinh Tran](https://github.com/trvinh)

# Contact
For further support or bug reports please contact: tran@bio.uni-frankfurt.de
