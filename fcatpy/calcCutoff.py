# -*- coding: utf-8 -*-

#######################################################################
#  Copyright (C) 2020 Vinh Tran
#
#  Calculate FAS cutoff for each core ortholog group of the core set
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
from pathlib import Path
from Bio import SeqIO
import subprocess
import multiprocessing as mp
import shutil
from tqdm import tqdm
# from ete3 import NCBITaxa
# import re
# from datetime import datetime

def checkFileExist(file):
    if not os.path.exists(os.path.abspath(file)):
        sys.exit('%s not found' % file)

def prepareJob(coreDir, coreSet, annoDir, blastDir):
    groups = os.listdir(coreDir + '/core_orthologs/' + coreSet)
    fasJobs = []
    if len(groups) > 0:
        for groupID in groups:
            group = '%s/core_orthologs/%s/%s' % (coreDir, coreSet, groupID)
            if os.path.isdir(group):
                print(group)
                query = '%s/%s.fa' % (group, groupID)
                annoDirTmp = '%s/fas_dir/annotation_dir/' % (group)
                Path(annoDirTmp).mkdir(parents=True, exist_ok=True)
                outDir = '%s/fas_dir/outTmp/' % (group)
                Path(outDir).mkdir(parents=True, exist_ok=True)
                for s in SeqIO.parse(query, 'fasta'):
                    ref = s.id.split('|')[1]
                    if not os.path.exists('%s/%s.json' % (annoDirTmp, ref)):
                        if os.path.exists('%s/%s.json' % (annoDir, ref)):
                            src = '%s/%s.json' % (annoDir, ref)
                            dst = '%s/%s.json' % (annoDirTmp, ref)
                            os.symlink(src, dst)
                    refGenome = '%s/%s/%s.fa' % (blastDir, ref, ref)
                    checkFileExist(refGenome)
                    fasJobs.append([s.id, s.seq, ref, query, annoDirTmp, outDir, refGenome])
    return(fasJobs)

# def parseFasOut(fasOut, refSpec):
#     score = 0
#     with open(fasOut, 'r') as file:
#         for l in file.readlines():
#             tmp = l.split('\t')
#             if refSpec in tmp[0]:
#                 if not refSpec in tmp[1]:



def calcFAS(args):
    (seqId, seq, refSpec, query, annoDir, outputDir, ref) = args
    # write seed fasta file
    seed = '%s/seed_%s.fa' % (outputDir, refSpec)
    tmpFile = open(seed, 'w')
    tmpFile.write(str('>' + seqId + '\n' + seq))
    tmpFile.close()
    # calculate fas scores for seed vs all
    fasCmd = 'calcFAS -s %s -q %s -a %s -o %s -n %s --bidirectional --domain -r %s' % (seed, query, annoDir, outputDir, refSpec, ref)
    try:
        subprocess.run([fasCmd], shell=True, check=True)
    except:
        print('Problem occurred while running calcFAS\n%s' % fasCmd)

def main():
    version = '0.0.1'
    parser = argparse.ArgumentParser(description='You are running calcCutoff version ' + str(version) + '.')
    required = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('optional arguments')
    required.add_argument('-d', '--coreDir', help='Path to core set directory, where folder core_orthologs can be found', action='store', default='', required=True)
    required.add_argument('-s', '--coreSet', help='Name of core set, which is subfolder within coreDir/core_orthologs/ directory', action='store', default='', required=True)
    optional.add_argument('-a', '--annoDir', help='Path to FAS annotation directory', action='store', default='')
    optional.add_argument('-b', '--blastDir', help='path to BLAST directory of all core species', action='store', default='')
    # optional.add_argument('-v', '--verProt', help='Proteome version', action='store', default=1, type=int)
    # optional.add_argument('-c', '--coreTaxa', help='Include this taxon to core taxa (i.e. taxa in blast_dir folder)', action='store_true', default=False)
    # optional.add_argument('-a', '--noAnno', help='Do NOT annotate this taxon using annoFAS', action='store_true', default=False)
    # optional.add_argument('--oldFAS', help='Use old verion of FAS (annoFAS ≤ 1.2.0)', action='store_true', default=False)
    optional.add_argument('--cpus', help='Number of CPUs used for annotation. Default = 4', action='store', default=4, type=int)
    # optional.add_argument('--replace', help='Replace special characters in sequences by "X"', action='store_true', default=False)
    # optional.add_argument('--delete', help='Delete special characters in sequences', action='store_true', default=False)
    # optional.add_argument('--force', help='Force overwrite existing data', action='store_true', default=False)

    args = parser.parse_args()

    coreDir = os.path.abspath(args.coreDir)
    coreSet = args.coreSet
    checkFileExist(coreDir + '/core_orthologs/' + coreSet)
    annoDir = os.path.abspath(args.annoDir)
    if annoDir == '':
        annoDir = '%s/weight_dir' % coreDir
    blastDir = os.path.abspath(args.blastDir)
    if blastDir == '':
        blastDir = '%s/blast_dir' % coreDir
    cpus = args.cpus
    # force = args.force

    fasJobs = prepareJob(coreDir, coreSet, annoDir, blastDir)
    pool = mp.Pool(cpus)
    fasOut = []
    for _ in tqdm(pool.imap_unordered(calcFAS, fasJobs), total=len(fasJobs)):
        fasOut.append(_)
    # for j in fasJobs:
    #     calcFAS(j)


if __name__ == '__main__':
    main()
