# -*- coding: utf-8 -*-

#######################################################################
#  Copyright (C) 2020 Vinh Tran
#
#  This script show possible reference taxa for each core set
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
from ete3 import NCBITaxa
import fcat.functions as fcatFn

def getNcbiName(taxonName):
    ncbi = NCBITaxa()
    taxId = taxonName.split('@')[1]
    try:
        name = ncbi.get_taxid_translator([taxId])[int(taxId)]
    except:
        name = taxonName
    return(name)

def printRefTaxa(args):
    coreDir = os.path.abspath(args.coreDir)
    coreSet = args.coreSet
    fcatFn.checkFileExist(coreDir + '/core_orthologs/' + coreSet, '')
    taxaFile = '%s/core_orthologs/%s/done.txt' % (coreDir, coreSet)
    fcatFn.checkFileExist(taxaFile, '')
    print('##### Taxa in the core sets, which can be used as reference species #####')
    for taxon in fcatFn.readFile(taxaFile):
        print('%s\t%s' % (taxon, getNcbiName(taxon)))
    sys.exit()

def main():
    version = '0.0.28'
    parser = argparse.ArgumentParser(description='You are running fcat version ' + str(version) + '.')
    required = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('optional arguments')
    required.add_argument('-d', '--coreDir', help='Path to core set directory, where folder core_orthologs can be found', action='store', default='', required=True)
    required.add_argument('-c', '--coreSet', help='Name of core set, which is subfolder within coreDir/core_orthologs/ directory', action='store', default='', required=True)

    args = parser.parse_args()
    printRefTaxa(args)

if __name__ == '__main__':
    main()
