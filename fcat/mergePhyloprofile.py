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
import fcat.functions as fcatFn

def mergePP(args):
    coreDir = os.path.abspath(args.coreDir)
    coreSet = args.coreSet
    outDir = os.path.abspath(args.outDir)
    if not 'fcatOutput' in outDir:
        outDir = outDir + '/fcatOutput/' + coreSet
    else:
        if not coreSet in outDir:
            outDir = outDir + '/' +coreSet
    fcatFn.checkFileExist(outDir,'')
    coreTaxaId = []
    for coreSpec in os.listdir('%s/blast_dir' % coreDir):
        coreTaxaId.append(coreSpec.split('@')[1])

    mode1out = []
    mode3out = []
    mode2out = []
    domainOut = []

    for query in os.listdir(outDir):
        if os.path.isdir(outDir + '/' + query):
            if '@' in query:
                if not query.split('@')[1] in coreTaxaId:
                    if os.path.exists('%s/%s/phyloprofileOutput/%s_mode1.phyloprofile' % (outDir, query, query)):
                        with open('%s/%s/phyloprofileOutput/%s_mode1.phyloprofile' % (outDir, query, query),'r') as fi:
                            mode1out += fi.readlines()[1:]
                    if os.path.exists('%s/%s/phyloprofileOutput/%s_mode2.phyloprofile' % (outDir, query, query)):
                        with open('%s/%s/phyloprofileOutput/%s_mode2.phyloprofile' % (outDir, query, query),'r') as fi:
                            mode2out += fi.readlines()[1:]
                    if os.path.exists('%s/%s/phyloprofileOutput/%s_mode3.phyloprofile' % (outDir, query, query)):
                        with open('%s/%s/phyloprofileOutput/%s_mode3.phyloprofile' % (outDir, query, query),'r') as fi:
                            mode3out += fi.readlines()[1:]
                    if os.path.exists('%s/%s/phyloprofileOutput/%s.domains' % (outDir, query, query)):
                        with open('%s/%s/phyloprofileOutput/%s.domains' % (outDir, query, query),'r') as fi:
                            domainOut += fi.readlines()

    if len(mode1out) > 0:
        mode1out = set(mode1out)
        with open('%s/%s_mode1.phyloprofile' % (outDir, coreSet),'w') as fo:
            fo.write('geneID\tncbiID\torthoID\tFAS\tAssessment\n')
            fo.write(''.join(mode1out))
    if len(mode2out) > 0:
        mode2out = set(mode2out)
        with open('%s/%s_mode2.phyloprofile' % (outDir, coreSet),'w') as fo:
            fo.write('geneID\tncbiID\torthoID\tFAS\tAssessment\n')
            fo.write("".join(mode2out))
    if len(mode3out) > 0:
        mode3out = set(mode3out)
        with open('%s/%s_mode3.phyloprofile' % (outDir, coreSet),'w') as fo:
            fo.write('geneID\tncbiID\torthoID\tFAS\tAssessment\n')
            fo.write("".join(mode3out))
    if len(domainOut) > 0:
        domainOut = set(domainOut)
        with open('%s/%s.domains' % (outDir, coreSet),'w') as fo:
            fo.write("".join(domainOut))

def main():
    version = '0.0.25'
    parser = argparse.ArgumentParser(description='You are running fcat version ' + str(version) + '.')
    parser.add_argument('-d', '--coreDir', help='Path to core set directory, where folder core_orthologs can be found', action='store', default='', required=True)
    parser.add_argument('-c', '--coreSet', help='Name of core set, which is subfolder within coreDir/core_orthologs/ directory', action='store', default='', required=True)
    parser.add_argument('-o', '--outDir', help='Path to output directory', action='store', default='', required=True)
    args = parser.parse_args()
    mergePP(args)

if __name__ == '__main__':
    main()
