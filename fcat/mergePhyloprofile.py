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
# from Bio import SeqIO
# import subprocess
# import multiprocessing as mp
# import shutil
# from tqdm import tqdm
# import time
# import statistics
# import collections

def checkFileExist(file, msg):
    if not os.path.exists(os.path.abspath(file)):
        sys.exit('%s not found! %s' % (file, msg))

def readFile(file):
    if os.path.exists(file):
        with open(file, 'r') as f:
            lines = f.read().splitlines()
            f.close()
            return(lines)
    else:
        sys.exit('%s not found' % file)

def mergePP(args):
    coreDir = os.path.abspath(args.coreDir)
    coreSet = args.coreSet
    outDir = os.path.abspath(args.outDir)
    if not 'fcatOutput' in outDir:
        outDir = outDir + '/fcatOutput/' + coreSet
    checkFileExist(outDir,'')
    coreTaxaId = []
    for coreSpec in os.listdir('%s/blast_dir' % coreDir):
        coreTaxaId.append(coreSpec.split('@')[1])
    mode1out = []
    mode2out = []
    mode3out = []
    domainOut = []
    for query in os.listdir(outDir):
        if os.path.isdir(outDir + '/' + query):
            if not query.split('@')[1] in coreTaxaId:
                if os.path.exists('%s/%s/phyloprofileOutput/%s_mode1.phyloprofile' % (outDir, query, query)):
                    for line in readFile('%s/%s/phyloprofileOutput/%s_mode1.phyloprofile' % (outDir, query, query)):
                        if not line in mode1out:
                            mode1out.append(line)
                if os.path.exists('%s/%s/phyloprofileOutput/%s_mode2.phyloprofile' % (outDir, query, query)):
                    for line in readFile('%s/%s/phyloprofileOutput/%s_mode2.phyloprofile' % (outDir, query, query)):
                        if not line in mode2out:
                            mode2out.append(line)
                if os.path.exists('%s/%s/phyloprofileOutput/%s_mode3.phyloprofile' % (outDir, query, query)):
                    for line in readFile('%s/%s/phyloprofileOutput/%s_mode3.phyloprofile' % (outDir, query, query)):
                        if not line in mode3out:
                            mode3out.append(line)
                if os.path.exists('%s/%s/phyloprofileOutput/%s.domains' % (outDir, query, query)):
                    for line in readFile('%s/%s/phyloprofileOutput/%s.domains' % (outDir, query, query)):
                        if not line in domainOut:
                            domainOut.append(line)
    if len(mode1out) > 0:
        mode1File = open('%s/%s_mode1.phyloprofile' % (outDir, coreSet), 'w')
        mode1File.write('%s\n' % '\n'.join(mode1out))
        mode1File.close()
    if len(mode2out) > 0:
        mode2File = open('%s/%s_mode2.phyloprofile' % (outDir, coreSet), 'w')
        mode2File.write('%s\n' % '\n'.join(mode2out))
        mode2File.close()
    if len(mode3out) > 0:
        mode3File = open('%s/%s_mode3.phyloprofile' % (outDir, coreSet), 'w')
        mode3File.write('%s\n' % '\n'.join(mode3out))
        mode3File.close()
    if len(domainOut) > 0:
        domainFile = open('%s/%s.domains' % (outDir, coreSet), 'w')
        domainFile.write('%s\n' % '\n'.join(domainOut))
        domainFile.close()

def main():
    version = '0.0.19'
    parser = argparse.ArgumentParser(description='You are running fcat version ' + str(version) + '.')
    parser.add_argument('-d', '--coreDir', help='Path to core set directory, where folder core_orthologs can be found', action='store', default='', required=True)
    parser.add_argument('-c', '--coreSet', help='Name of core set, which is subfolder within coreDir/core_orthologs/ directory', action='store', default='', required=True)
    parser.add_argument('-o', '--outDir', help='Path to output directory', action='store', default='', required=True)
    args = parser.parse_args()
    mergePP(args)

if __name__ == '__main__':
    main()
