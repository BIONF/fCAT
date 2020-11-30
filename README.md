# fCAT
[![PyPI version](https://badge.fury.io/py/fcat.svg)](https://pypi.org/project/fcat/)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Build Status](https://travis-ci.com/BIONF/fCAT.svg?branch=main)](https://travis-ci.com/BIONF/fCAT)

Python package for fCAT, a **f**eature-aware **C**ompleteness **A**ssessment **T**ool

# Table of Contents
* [How to install](#how-to-install)
* [Usage](#usage)
* [Output](#output)
* [fCAT score modes](#fcat-score-modes)
* [Bugs](#bugs)
* [Contributors](#contributors)
* [Contact](#contact)

# How to install

*fCAT* tool is distributed as a python package called *fcat*. It is compatible with [Python ≥ v3.7](https://www.python.org/downloads/).

You can install *fcat* using `pip`:
```
python3 -m pip install fcat
```

or, in case you do not have admin rights, and don't use package systems like [Anaconda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/) to manage environments you need to use the `--user` option:
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

where **eukaryota** is name of the fCAT core set (equivalent to [BUSCO](https://busco.ezlab.org/) set); **HOMSA@9606@2** is the reference species from that core set that will be used for the ortholog search; **query** is the name of species of interest. If `--annoQuery` not specified, *fCAT* fill do the feature annotation for the query proteins using [FAS tool](https://github.com/BIONF/FAS).

# Output
You will find the output in the */path/to/fcat/output/fcatOutput/eukaryota/query/* folder, where */path/to/fcat/output/* could be your current directory if you not specified `--outDir` when running `fcat`. The following important output files/folders can be found:

- *all_summary.txt*: summary of the completeness assessment using all 4 score modes
- *all_full.txt*: the complete assessment of 4 score modes in tab delimited file
- *fdogOutput.tar.gz*: a zipped file of the ortholog search result
- *mode_1*, *mode_2*, *mode_3* and *mode_4*: detailed output for each score mode
- *phyloprofileOutput*: folder contains output phylogenetic profile data that can be used with [PhyloProfile tool](https://github.com/BIONF/PhyloProfile)

Besides, if you have already run *fCAT* for several query taxa with the same fCAT core set, you can find the merged phylogentic profiles for all of those taxa within the corresponding core set output (e.g. _/path/to/fcat/output/fcatOutput/eukaryota/*.phyloprofile_).

# fCAT score modes

The table below explains how the *specific ortholog group cutoffs* for each fCAT core set were calculated, and which *value of the query ortholog* is used to assess its completeness, or more precisely, its functional equivalence to the ortholog group it belongs to. If the value of a query ortholog is *not less than* its ortholog group cutoff, that group will be evaluated as **similar** or **complete**. In case co-orthologs have been predicted, the assessment for the core group will be **duplicated**. Depending on the value of each single ortholog, a *duplicated* group can be seen as **duplicated (similar)** or **duplicated (dissimilar)** in the full report (e.g. *all_full.txt*).

| Score mode | Cutoff | Value |
|---|---|---|
| Mode 1 | Mean of FAS scores between all core orthologs | Mean of FAS scores between query ortholog and all core proteins |
| Mode 2 | Mean of FAS scores between refspec and all other core orthologs | Mean of FAS scores between query ortholog and refspec protein |
| Mode 3 | The lower bound of the confidence interval calculated by the distribution of all-vs-all FAS score in a core group | Mean of FAS scores between query ortholog and refspec protein |
| Mode 4 | Mean and standard deviation of all core protein lengths | Length of query ortholog |

*Note: __FAS scores__ are bidirectional FAS scors; __core protein__ or __core ortholog__ is protein in the core ortholog groups; __query protein__ or __query ortholog__ is ortholog protein of query species; __refspec__ is the specified reference species*

# Bugs
Any bug reports or comments, suggestions are highly appreciated. Please [open an issue on GitHub](https://github.com/BIONF/fCAT/issues/new) or be in touch via email.

# Contributors
- [Vinh Tran](https://github.com/trvinh)
- [Giang Nguyen](https://github.com/giangnguyen0709)

# Contact
For further support or bug reports please contact: tran@bio.uni-frankfurt.de
