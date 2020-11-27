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
import shutil
import time
import statistics
import glob

def checkFileExist(file, msg):
    if not os.path.exists(os.path.abspath(file)):
        sys.exit('%s not found! %s' % (file, msg))

def roundTo4(number):
    return("%.4f" % round(float(number), 4))

def readFile(file):
    if os.path.exists(file):
        with open(file, 'r') as f:
            lines = f.readlines()
            f.close()
            return(lines)
    else:
        sys.exit('%s not found' % file)

def removeDup(file):
    if os.path.exists(file):
        lines_seen = set() # holds lines already seen
        outfile = open(file+'.temp', 'w')
        for line in open(file, 'r'):
            if line not in lines_seen: # not a duplicate
                outfile.write(line)
                lines_seen.add(line)
        outfile.close()
        os.replace(file+'.temp', file)

def readRefspecFile(refspecFile):
    groupRefspec = {}
    for line in readFile(refspecFile):
        groupRefspec[line.split('\t')[0]] = line.split('\t')[1].strip()
    return(groupRefspec)

def outputMode(outDir, coreSet, queryID, force, approach):
    phyloprofileDir = '%s/fcatOutput/%s/%s/phyloprofileOutput' % (outDir, coreSet, queryID)
    Path(phyloprofileDir).mkdir(parents=True, exist_ok=True)
    if not os.path.exists('%s/%s_%s.phyloprofile' % (phyloprofileDir, coreSet, approach)):
        mode = 3
    else:
        if force:
            mode = 1
        else:
            mode = 0
    return(mode, phyloprofileDir)

def createProfile23(coreDir, outDir, coreSet, queryID, force):
    # output files
    (mode, phyloprofileDir) = outputMode(outDir, coreSet, queryID, force, 'other')
    if mode == 1 or mode == 3:
        finalPhyloprofile = open('%s/mode2.phyloprofile' % (phyloprofileDir), 'w')
        finalPhyloprofile.write('geneID\tncbiID\torthoID\tFAS\n')
        finalDomain = open('%s/mode23.domains' % (phyloprofileDir), 'w')
    # parse into phyloprofile file
    fdogOutDir = '%s/fcatOutput/%s/%s/fdogOutput' % (outDir, coreSet, queryID)
    out = os.listdir(fdogOutDir)
    for refSpec in out:
        if os.path.isdir(fdogOutDir + '/' + refSpec):
            refDir = fdogOutDir + '/' + refSpec
            groups = os.listdir(refDir)
            # move to phyloprofile output dir
            if not mode == 0:
                if os.path.exists('%s/%s.phyloprofile' % (refDir, refSpec)):
                    for line in readFile('%s/%s.phyloprofile' % (refDir, refSpec)):
                        if queryID in line:
                            tmpQuery = line.split('\t')
                            # statistics.mean((float(line.split('\t')[3], float(line.split('\t')[4]))
                            finalPhyloprofile.write('%s\t%s\t%s\t%s\n' % (tmpQuery[0], tmpQuery[1], tmpQuery[2], roundTo4(statistics.mean((float(line.split('\t')[3]), float(line.split('\t')[4]))))))
                # append profile of core sequences
                for groupID in groups:
                    coreFasDir = '%s/core_orthologs/%s/%s/fas_dir/fasscore_dir' % (coreDir, coreSet, groupID)
                    for fasFile in glob.glob('%s/*.tsv' % coreFasDir):
                        if queryID.split('@')[1] == refSpec.split('@')[1]:
                            if not refSpec in fasFile:
                                for fLine in readFile(fasFile):
                                    if refSpec in fLine.split('\t')[0]:
                                        tmp = fLine.split('\t')
                                        revFAS = 0
                                        revFile = '%s/%s.tsv' % (coreFasDir, tmp[0].split('|')[1])
                                        for revLine in readFile(revFile):
                                            if tmp[1] == revLine.split('\t')[0]:
                                                revFAS = revLine.split('\t')[2].split('/')[0]
                                        coreLine = '%s\t%s\t%s\t%s\n' % (groupID, 'ncbi' + str(tmp[1].split('|')[1].split('@')[1]), tmp[1], roundTo4(statistics.mean((float(tmp[2].split('/')[0]), float(revFAS)))))
                                        finalPhyloprofile.write(coreLine)
                        else:
                            for fLine in readFile(fasFile):
                                if refSpec in fLine.split('\t')[0]:
                                    tmp = fLine.split('\t')
                                    revFAS = 0
                                    revFile = '%s/%s.tsv' % (coreFasDir, tmp[0].split('|')[1])
                                    for revLine in readFile(revFile):
                                        if tmp[1] == revLine.split('\t')[0]:
                                            revFAS = revLine.split('\t')[2].split('/')[0]
                                    coreLine = '%s\t%s\t%s\t%s\n' % (groupID, 'ncbi' + str(tmp[1].split('|')[1].split('@')[1]), tmp[1], roundTo4(statistics.mean((float(tmp[2].split('/')[0]), float(revFAS)))))
                                    finalPhyloprofile.write(coreLine)
    # parse domain file
    # if os.path.exists('%s/%s_forward.domains' % (refDir, refSpec)):
    #     for line in readFile('%s/%s_forward.domains' % (refDir, refSpec)):
    #         finalDomain.write(line)
    # finalize
    if not mode == 0:
        finalPhyloprofile.close()
        finalDomain.close()
    # delete duplicate lines
    removeDup('%s/mode2.phyloprofile' % (phyloprofileDir))
    removeDup('%s/mode23.domains' % (phyloprofileDir))
    shutil.copy('%s/mode2.phyloprofile' % (phyloprofileDir), '%s/mode3.phyloprofile' % (phyloprofileDir))

def createProfile1(coreDir, outDir, coreSet, queryID, force, groupRefspec):
    # output files
    phyloprofileDir = '%s/fcatOutput/%s/%s/phyloprofileOutput' % (outDir, coreSet, queryID)
    (mode, phyloprofileDir) = outputMode(outDir, coreSet, queryID, force, 'mode1')
    if mode == 1 or mode == 3:
        finalFa = open('%s/%s.mod.fa' % (phyloprofileDir, coreSet), 'w')
        finalDomain = open('%s/mode1.domains' % (phyloprofileDir), 'w')
        finalPhyloprofile = open('%s/mode1.phyloprofile' % (phyloprofileDir), 'w')
        finalPhyloprofile.write('geneID\tncbiID\torthoID\tFAS\n')
        finalLen = open('%s/length.phyloprofile' % (phyloprofileDir), 'w')
        finalLen.write('geneID\tncbiID\torthoID\tLength\n')
    # parse into phyloprofile files
    fdogOutDir = '%s/fcatOutput/%s/%s/fdogOutput' % (outDir, coreSet, queryID)
    mergedFa = '%s/%s_all.extended.fa' % (fdogOutDir, queryID)
    if not mode == 0:
        # get fas scores for each group
        groupScoreFwd = {} # all fwd fas scores of query ortholog vs core proteins
        groupScoreRev = {} # all rev fas scores of query ortholog vs core proteins
        groupOrtho = {} # query ortholog ID of each group, used for orthoID column of phyloprofile output
        for line in readFile('%s/%s_all.phyloprofile' % (fdogOutDir, queryID)):
            if not line.split('\t')[0] == 'geneID':
                groupID = line.split('\t')[0]
                if not groupID in groupScoreFwd:
                    groupScoreFwd[groupID] = []
                    groupScoreRev[groupID] = []
                if queryID in line.split('\t')[2]:
                    groupOrtho[groupID] = line.split('\t')[2]
                else:
                    groupScoreFwd[groupID].append(float(line.split('\t')[3]))
                    groupScoreRev[groupID].append(float(line.split('\t')[4]))

        # calculate mean fas score for query ortholog against all core proteins and add to phyloprofile output
        for groupID in groupOrtho:
            groupIDmod = '_'.join(groupID.split('_')[1:])
            groupOrthoMod = '_'.join(groupOrtho[groupID].split('_')[1:])
            newline = '%s\t%s\t%s\t%s\n' % (groupIDmod, 'ncbi' + str(queryID.split('@')[1]), groupOrthoMod, roundTo4(statistics.mean((statistics.mean(groupScoreFwd[groupID]), statistics.mean(groupScoreRev[groupID])))))
            finalPhyloprofile.write(newline)
            # append profile of core sequences for this group
            meanCoreFile = '%s/core_orthologs/%s/%s/fas_dir/cutoff_dir/2.cutoff' % (coreDir, coreSet, groupIDmod)
            for tax in readFile(meanCoreFile):
                if not tax.split('\t')[0] == 'taxa':
                    # not include core taxon that have the same taxonomy ID as query
                    if queryID.split('@')[1] == groupRefspec[groupIDmod].split('@')[1]:
                        if not tax.split('\t')[0] == groupRefspec[groupIDmod]:
                            ppCore = '%s\t%s\t%s|1\t%s\n' % (groupIDmod, 'ncbi' + str(tax.split('\t')[0].split('@')[1]), tax.split('\t')[2].strip(), roundTo4(tax.split('\t')[1]))
                            finalPhyloprofile.write(ppCore)
                    else:
                        ppCore = '%s\t%s\t%s|1\t%s\n' % (groupIDmod, 'ncbi' + str(tax.split('\t')[0].split('@')[1]), tax.split('\t')[2].strip(), roundTo4(tax.split('\t')[1]))
                        finalPhyloprofile.write(ppCore)
        # add missing groups
        if os.path.exists('%s/fcatOutput/%s/%s/missing.txt' % (outDir, coreSet, queryID)):
            for missingGr in readFile('%s/fcatOutput/%s/%s/missing.txt' % (outDir, coreSet, queryID)):
                meanCoreFile = '%s/core_orthologs/%s/%s/fas_dir/cutoff_dir/2.cutoff' % (coreDir, coreSet, missingGr.strip())
                for tax in readFile(meanCoreFile):
                    if not tax.split('\t')[0] == 'taxa':
                        # not include core taxon that have the same taxonomy ID as query
                        if queryID.split('@')[1] == groupRefspec[missingGr.strip()].split('@')[1]:
                            if not tax.split('\t')[0] == groupRefspec[missingGr.strip()]:
                                ppCore = '%s\t%s\t%s|1\t%s\n' % (missingGr.strip(), 'ncbi' + str(tax.split('\t')[0].split('@')[1]), tax.split('\t')[2].strip(), roundTo4(tax.split('\t')[1]))
                                finalPhyloprofile.write(ppCore)
                        else:
                            ppCore = '%s\t%s\t%s|1\t%s\n' % (missingGr.strip(), 'ncbi' + str(tax.split('\t')[0].split('@')[1]), tax.split('\t')[2].strip(), roundTo4(tax.split('\t')[1]))
                            finalPhyloprofile.write(ppCore)
        finalPhyloprofile.close()

        # length phyloprofile file and final fasta file
        for s in SeqIO.parse(mergedFa, 'fasta'):
            idMod = '_'.join(s.id.split('_')[1:])
            if queryID.split('@')[1] == groupRefspec[groupIDmod].split('@')[1]:
                if not idMod.split('|')[1] == groupRefspec[idMod.split('|')[0]]:
                    finalFa.write('>%s\n%s\n' % (idMod, s.seq))
                    ppLen = '%s\t%s\t%s\t%s\n' % (idMod.split('|')[0], 'ncbi' + str(idMod.split('|')[1].split('@')[1]), idMod, len(s.seq))
                    finalLen.write(ppLen)
            else:
                finalFa.write('>%s\n%s\n' % (idMod, s.seq))
                ppLen = '%s\t%s\t%s\t%s\n' % (idMod.split('|')[0], 'ncbi' + str(idMod.split('|')[1].split('@')[1]), idMod, len(s.seq))
                finalLen.write(ppLen)
        finalFa.close()
        finalLen.close()

        # parse domain file
        # shutil.copyfileobj(open('%s/%s_all_forward.domains' % (fdogOutDir, queryID), 'rb'), finalFwdDomain)
        # finalFwdDomain.close()
        # finalDomain = open('%s/mode1.domains' % (phyloprofileDir), 'w')
        # for domains in readFile('%s/FAS_forward.domains' % (phyloprofileDir)):
        #     tmp = domains.split('\t')
        #     mGroup = '_'.join(tmp[0].split('#')[0].split('_')[1:])
        #     mQuery = '_'.join(tmp[0].split('#')[1].split('_')[1:])
        #     mSeed = '_'.join(tmp[1].split('_')[1:])
        #     # if queryID in mSeed:
        #     domainLine = '%s\t%s\t%s\t%s\t%s\t%s\tNA\tN\n' % (mGroup+'#'+mQuery, mSeed, tmp[2], tmp[3], tmp[4], tmp[5])
        #     finalDomain.write(domainLine)
        # finalDomain.close()
        # os.remove('%s/FAS_forward.domains' % (phyloprofileDir))
    # delete duplicate lines
    removeDup('%s/mode1.phyloprofile' % (phyloprofileDir))
    removeDup('%s/mode1.domains' % (phyloprofileDir))
    removeDup('%s/length.phyloprofile' % (phyloprofileDir))

def deleteFolder(folder):
    if os.path.exists(folder):
        if os.path.isfile(folder):
            os.remove(folder)
        else:
            shutil.rmtree(folder)

def createPhyloProfile(args):
    coreDir = os.path.abspath(args.coreDir)
    coreSet = args.coreSet
    checkFileExist(coreDir + '/core_orthologs/' + coreSet, '')
    queryID = args.queryID
    outDir = args.outDir
    if outDir == '':
        outDir = os.getcwd()
    else:
        Path(outDir).mkdir(parents=True, exist_ok=True)
    force = args.forceProfile
    keep = args.keep

    # check old output files
    fcatOut = '%s/fcatOutput/%s/%s' % (outDir, coreSet, queryID)

    if force:
        deleteFolder('%s/phyloprofileOutput' % fcatOut)
    if not os.path.exists('%s/phyloprofileOutput/mode1.phyloprofile' % fcatOut):
        if not os.path.exists('%s/fdogOutput'):
            if os.path.exists('%s/fdogOutput.tar.gz' % fcatOut):
                shutil.unpack_archive('%s/fdogOutput.tar.gz' % fcatOut, fcatOut + '/', 'gztar')
            else:
                sys.exit('No ortholog output found!')
        if os.path.exists('%s/last_refspec.txt' % fcatOut):
            groupRefspec = readRefspecFile('%s/last_refspec.txt' % fcatOut)
            createProfile1(coreDir, outDir, coreSet, queryID, force, groupRefspec)
            createProfile23(coreDir, outDir, coreSet, queryID, force)
        else:
            sys.exit('No last_refspec.txt file found!')

        if keep == False:
            print('Cleaning up...')
            if os.path.exists('%s/fdogOutput/' % (fcatOut)):
                shutil.rmtree('%s/fdogOutput/' % (fcatOut))

def main():
    version = '0.0.10'
    parser = argparse.ArgumentParser(description='You are running fcat version ' + str(version) + '.')
    required = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('optional arguments')
    required.add_argument('-d', '--coreDir', help='Path to core set directory, where folder core_orthologs can be found', action='store', default='', required=True)
    required.add_argument('-c', '--coreSet', help='Name of core set, which is subfolder within coreDir/core_orthologs/ directory', action='store', default='', required=True)
    required.add_argument('--queryID', help='ID of taxon of interest (e.g. HUMAN@9606@3)', action='store', default='', type=str)
    optional.add_argument('-o', '--outDir', help='Path to output directory', action='store', default='')
    optional.add_argument('--force', help='Force overwrite existing data', action='store_true', default=False)
    optional.add_argument('--keep', help='Keep temporary phyloprofile data', action='store_true', default=False)
    args = parser.parse_args()

    start = time.time()
    searchOrtho(args)
    ende = time.time()
    print('Finished in ' + '{:5.3f}s'.format(ende-start))

if __name__ == '__main__':
    main()
