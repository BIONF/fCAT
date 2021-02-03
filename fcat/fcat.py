# -*- coding: utf-8 -*-

#######################################################################
#  Copyright (C) 2020 Vinh Tran
#
#  This script will do the complete completeness check for a given
#  protein set.
#
#  This script is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License <http://www.gnu.org/licenses/> for
#  more details
#
#  Contact: tran@bio.uni-frankfurt.de
#
#######################################################################

import sys
import os
import argparse
import time
import multiprocessing as mp
import shutil
# import fcat.calcCutoff as fcatC
import fcat.searchOrtho as fcatO
import fcat.assessCompleteness as fcatR
import fcat.createPhyloprofile as fcatP
import fcat.mergePhyloprofile as fcatM
import fcat.functions as fcatFn

def checkRefspec(coreDir, coreSet, refspecList):
    taxaFile = '%s/core_orthologs/%s/done.txt' % (coreDir, coreSet)
    fcatFn.checkFileExist(taxaFile, '')
    allTaxa = fcatFn.readFile(taxaFile)
    invalid = []
    for refspec in refspecList:
        if not refspec in allTaxa:
            invalid.append(refspec)
    if len(invalid) > 0:
        print('ERROR: Invalid refspec found: %s' % '; '.join(invalid))
        print('Please use fcat.showTaxa for a list of valid reference taxa!')
        sys.exit()

def fcat(args):
    # calculate group specific cutoffs
    # print('##### Calculating group specific cutoffs...')
    # fcatC.calcGroupCutoff(args)

    # check for valid refspec
    checkRefspec(args.coreDir, args.coreSet, str(args.refspecList).split(","))

    # search for orthologs and create phylognetic profile files
    print('##### Searching for orthologs...')
    fcatO.searchOrtho(args)

    # get queryID
    outDir = args.outDir
    if outDir == '':
        outDir = os.getcwd()
    else:
        Path(outDir).mkdir(parents=True, exist_ok=True)
    annoDir = args.annoDir
    if annoDir == '':
        annoDir = '%s/weight_dir' % args.coreDir
    annoDir = os.path.abspath(annoDir)
    cpus = args.cpus
    if cpus >= mp.cpu_count():
        cpus = mp.cpu_count()-1
    (doAnno, queryTaxId) = fcatO.checkQueryAnno(args.annoQuery, annoDir, args.taxid, args.querySpecies)
    args.queryID = fcatO.parseQueryFa(args.coreSet, os.path.abspath(args.querySpecies), args.annoQuery, str(args.taxid), outDir, doAnno, annoDir, cpus)
    if doAnno == False:
        if os.path.exists( '%s/query_%s.json' % (annoDir, queryTaxId)):
            os.remove('%s/query_%s.json' % (annoDir, queryTaxId))

    # create phyloprofile files
    print('##### Creating phylogenetic profiles....')
    fcatP.createPhyloProfile(args)

    # do completeness assessment
    print('##### Generating reports...')
    flag = fcatR.assessCompteness(args) # flag=0: no new ppFile needed

    # merge phyloprofile
    if flag == 1:
        print('##### Merging phylogenetic profiles of all query taxa...')
        fcatM.mergePP(args)

def main():
    version = '0.0.25'
    parser = argparse.ArgumentParser(description='You are running fcat version ' + str(version) + '.')
    required = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('optional arguments')
    required.add_argument('-d', '--coreDir', help='Path to core set directory, where folder core_orthologs can be found', action='store', default='', required=True)
    required.add_argument('-c', '--coreSet', help='Name of core set, which is subfolder within coreDir/core_orthologs/ directory', action='store', default='', required=True)
    required.add_argument('-r', '--refspecList', help='List of reference species', action='store', default='', required=True)
    required.add_argument('-q', '--querySpecies', help='Path to gene set for species of interest', action='store', default='', required=True)
    optional.add_argument('-o', '--outDir', help='Path to output directory', action='store', default='')
    optional.add_argument('-b', '--blastDir', help='Path to BLAST directory of all core species', action='store', default='')
    optional.add_argument('-a', '--annoDir', help='Path to FAS annotation directory', action='store', default='')
    optional.add_argument('--annoQuery', help='Path to FAS annotation for species of interest', action='store', default='')
    optional.add_argument('-i', '--taxid', help='Taxonomy ID of gene set for species of interest', action='store', default=0, type=int)
    optional.add_argument('-m', '--mode', help='Score cutoff mode. (0) all modes, (1) all-vs-all FAS scores, (2) mean FAS of refspec seed, (3) confidence interval of all group FAS scores, (4) mean and stdev of sequence length',
                            action='store', default=0, choices=[0,1,2,3,4], type=int)
    optional.add_argument('--cpus', help='Number of CPUs used for annotation. Default = 4', action='store', default=4, type=int)
    optional.add_argument('--force', help='Force overwrite existing ortholog search and assessment output', action='store_true', default=False)
    optional.add_argument('--forceProfile', help='Force overwrite phylogenetic profiles output', action='store_true', default=False)
    optional.add_argument('--noDomain', help='Not output feature domain files', action='store_true', default=False)
    optional.add_argument('--keep', help='Keep temporary phyloprofile data', action='store_true', default=False)
    optional.add_argument('--bidirectional', help=argparse.SUPPRESS, action='store_true', default=False)

    args = parser.parse_args()

    start = time.time()
    fcat(args)
    ende = time.time()
    print('##### Finished in ' + '{:5.3f}s'.format(ende-start))

if __name__ == '__main__':
    main()
