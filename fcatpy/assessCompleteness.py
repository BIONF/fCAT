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
import time
import statistics
import collections

def checkFileExist(file, msg):
    if not os.path.exists(os.path.abspath(file)):
        sys.exit('%s not found! %s' % (file, msg))

def readFile(file):
    if os.path.exists(file):
        with open(file, 'r') as f:
            lines = f.readlines()
            f.close()
            return(lines)
    else:
        return('NA')

def addToDict(dict, groupID, seqID, type):
    if not groupID in dict:
        dict[groupID] = '%s\t%s\t%s' % (groupID, type, seqID)
    else:
        tmp = dict[groupID].split('\t')
        dict[groupID] = '%s\tduplicated (%s)\t%s' % (tmp[0], tmp[1], tmp[2])
        dict[groupID] = '%s\n%s\tduplicated (%s)\t%s' % (dict[groupID], groupID, type, seqID)
    return(dict)

def mode1(ppFile, coreDir, coreSet, queryID):
    noCutoff = []
    assessment = {}
    for line in readFile(ppFile):
        groupID = line.split('\t')[0]
        if queryID in line.split('\t')[2]:
            # meanFas = statistics.mean((float(line.split('\t')[3]), float(line.split('\t')[4].strip())))
            meanFas = float(line.split('\t')[3].strip())
            scoreFile = '%s/core_orthologs/%s/%s/fas_dir/cutoff_dir/1.cutoff' % (coreDir, coreSet, groupID)
            if os.path.exists(scoreFile):
                meanGroup = 0
                for l in readFile(scoreFile):
                    if l.split('\t')[0] == 'mean':
                        meanGroup = float(l.split('\t')[1].strip())
                if meanFas >= meanGroup:
                    assessment = addToDict(assessment, groupID, line.split('\t')[2], 'similar')
                else:
                    assessment = addToDict(assessment, groupID, line.split('\t')[2], 'dissimilar')
            else:
                noCutoff.append(groupID)
    return(assessment, noCutoff)

def mode2(ppFile, coreDir, coreSet, queryID, outDir):
    noCutoff = []
    assessment = {}
    # get refspec for each group
    groupRefspec = {}
    refspecFile = '%s/%s/%s/last_refspec.txt' % (outDir, coreSet, queryID)
    for g in readFile(refspecFile):
        groupRefspec[g.split('\t')[0]] = g.split('\t')[1]
    # do assessment
    for line in readFile(ppFile):
        groupID = line.split('\t')[0]
        if queryID in line.split('\t')[2]:
            meanFas = statistics.mean((float(line.split('\t')[3]), float(line.split('\t')[4].strip())))
            scoreFile = '%s/core_orthologs/%s/%s/fas_dir/cutoff_dir/2.cutoff' % (coreDir, coreSet, groupID)
            if os.path.exists(scoreFile):
                meanRefspec = 0
                for l in readFile(scoreFile):
                    if l.split('\t')[0] == groupRefspec[groupID].strip():
                        meanRefspec = float(l.split('\t')[1].strip())
                if meanFas >= meanRefspec:
                    assessment = addToDict(assessment, groupID, line.split('\t')[2], 'similar')
                else:
                    assessment = addToDict(assessment, groupID, line.split('\t')[2], 'dissimilar')
            else:
                noCutoff.append(groupID)
    return(assessment, noCutoff)

def mode3(ppFile, coreDir, coreSet, queryID):
    noCutoff = []
    assessment = {}
    for line in readFile(ppFile):
        groupID = line.split('\t')[0]
        if queryID in line.split('\t')[2]:
            meanFas = statistics.mean((float(line.split('\t')[3]), float(line.split('\t')[4].strip())))
            scoreFile = '%s/core_orthologs/%s/%s/fas_dir/cutoff_dir/1.cutoff' % (coreDir, coreSet, groupID)
            if os.path.exists(scoreFile):
                LCL = 0
                UCL = 0
                for l in readFile(scoreFile):
                    if l.split('\t')[0] == 'LCL':
                        LCL = float(l.split('\t')[1].strip())
                    if l.split('\t')[0] == 'UCL':
                        UCL = float(l.split('\t')[1].strip())
                # if LCL <= meanFas <= UCL:
                if LCL <= meanFas:
                    assessment = addToDict(assessment, groupID, line.split('\t')[2], 'similar')
                else:
                    assessment = addToDict(assessment, groupID, line.split('\t')[2], 'dissimilar')
            else:
                noCutoff.append(groupID)
    return(assessment, noCutoff)

def mode4(ppFile, coreDir, coreSet, queryID):
    noCutoff = []
    assessment = {}
    for line in readFile(ppFile):
        groupID = line.split('\t')[0]
        if queryID in line.split('\t')[2]:
            length = float(line.split('\t')[3].strip())
            scoreFile = '%s/core_orthologs/%s/%s/fas_dir/cutoff_dir/1.cutoff' % (coreDir, coreSet, groupID)
            if os.path.exists(scoreFile):
                meanLen = 0
                stdevLen = 0
                for l in readFile(scoreFile):
                    if l.split('\t')[0] == 'meanLen':
                        meanLen = float(l.split('\t')[1].strip())
                    if l.split('\t')[0] == 'stdevLen':
                        stdevLen = float(l.split('\t')[1].strip())
                if not stdevLen == 0:
                    check = abs((length - meanLen) / stdevLen)
                    if check <= 1:
                        assessment = addToDict(assessment, groupID, line.split('\t')[2], 'complete')
                    else:
                        assessment = addToDict(assessment, groupID, line.split('\t')[2], 'fragmented')
                else:
                    noCutoff.append(groupID)
            else:
                noCutoff.append(groupID)
    return(assessment, noCutoff)

def mode5(ppFile, coreDir, coreSet, queryID):
    noCutoff = []
    assessment = {}
    for line in readFile(ppFile):
        groupID = line.split('\t')[0]
        # if queryID in line.split('\t')[2]:
        if line.split('\t')[1] == 'ncbi'+str(queryID.split('@')[1]):
            # meanFas = statistics.mean((float(line.split('\t')[3]), float(line.split('\t')[4].strip())))
            meanFas = float(line.split('\t')[3].strip())
            scoreFile = '%s/core_orthologs/%s/%s/fas_dir/cutoff_dir/3.cutoff' % (coreDir, coreSet, groupID)
            if os.path.exists(scoreFile):
                meanGroup = 0
                for l in readFile(scoreFile):
                    if l.split('\t')[0] == 'meanCons':
                        meanGroup = float(l.split('\t')[1].strip())
                if meanFas >= meanGroup:
                    assessment = addToDict(assessment, groupID, line.split('\t')[2], 'similar')
                else:
                    assessment = addToDict(assessment, groupID, line.split('\t')[2], 'dissimilar')
            else:
                noCutoff.append(groupID)
    return(assessment, noCutoff)

def mode6(ppFile, coreDir, coreSet, queryID, outDir):
    noCutoff = []
    assessment = {}
    # get refspec for each group
    groupRefspec = {}
    refspecFile = '%s/%s/%s/last_refspec.txt' % (outDir, coreSet, queryID)
    for g in readFile(refspecFile):
        groupRefspec[g.split('\t')[0]] = g.split('\t')[1]
    # do assessment
    for line in readFile(ppFile):
        groupID = line.split('\t')[0]
        # if queryID in line.split('\t')[2]:
        if line.split('\t')[1] == 'ncbi'+str(queryID.split('@')[1]):
            # meanFas = statistics.mean((float(line.split('\t')[3]), float(line.split('\t')[4].strip())))
            meanFas = float(line.split('\t')[3].strip())
            scoreFile = '%s/core_orthologs/%s/%s/fas_dir/cutoff_dir/4.cutoff' % (coreDir, coreSet, groupID)
            if os.path.exists(scoreFile):
                meanRefspec = 0
                for l in readFile(scoreFile):
                    if l.split('\t')[0] == groupRefspec[groupID].strip():
                        meanRefspec = float(l.split('\t')[1].strip())
                if meanFas >= meanRefspec:
                    assessment = addToDict(assessment, groupID, line.split('\t')[2], 'similar')
                else:
                    assessment = addToDict(assessment, groupID, line.split('\t')[2], 'dissimilar')
            else:
                noCutoff.append(groupID)
    return(assessment, noCutoff)

def doAssessment(ppDir, coreDir, coreSet, queryID, outDir, mode):
    assessment = {}
    noCutoff = []
    # assess completeness
    if mode == 1:
        ppFile = '%s/mode1.phyloprofile' % (ppDir)
        (assessment, noCutoff) = mode1(ppFile, coreDir, coreSet, queryID)
    elif mode == 2:
        ppFile = '%s/mode23.phyloprofile' % (ppDir)
        (assessment, noCutoff) = mode2(ppFile, coreDir, coreSet, queryID, outDir)
    elif mode == 3:
        ppFile = '%s/mode23.phyloprofile' % (ppDir)
        (assessment, noCutoff) = mode3(ppFile, coreDir, coreSet, queryID)
    elif mode == 4:
        ppFile = '%s/length.phyloprofile' % (ppDir)
        (assessment, noCutoff) = mode4(ppFile, coreDir, coreSet, queryID)
    elif mode == 5:
        ppFile = '%s/mode4.phyloprofile' % (ppDir)
        (assessment, noCutoff) = mode5(ppFile, coreDir, coreSet, queryID)
    elif mode == 6:
        ppFile = '%s/mode4.phyloprofile' % (ppDir)
        (assessment, noCutoff) = mode6(ppFile, coreDir, coreSet, queryID, outDir)
    # print full report
    stat = writeReport(assessment, outDir, coreDir, coreSet, queryID, mode)
    return(noCutoff, stat)

def writeReport(assessment, outDir, coreDir, coreSet, queryID, mode):
    missing = '%s/%s/%s/missing.txt' % (outDir, coreSet, queryID)
    ignored = '%s/%s/%s/ignored.txt' % (outDir, coreSet, queryID)
    Path('%s/%s/%s/mode_%s' % (outDir, coreSet, queryID, mode)).mkdir(parents=True, exist_ok=True)
    # write full report
    fullFile = open('%s/%s/%s/mode_%s/full.txt' % (outDir, coreSet, queryID, mode), 'w')
    for group in assessment:
        fullFile.write('%s\n' % assessment[group])
    for m in readFile(missing):
        fullFile.write(m.strip() + '\tmissing\n')
    for i in readFile(ignored):
        fullFile.write(i.strip() + '\tignored\n')
    fullFile.close()
    # write summary report
    summaryFile = open('%s/%s/%s/mode_%s/summary.txt' % (outDir, coreSet, queryID, mode), 'w')
    type = [x.split('\t')[1] for x in open('%s/%s/%s/mode_%s/full.txt' % (outDir, coreSet, queryID, mode)).readlines()]
    groupID = [x.split('\t')[0] for x in open('%s/%s/%s/mode_%s/full.txt' % (outDir, coreSet, queryID, mode)).readlines()]
    dup = [item for item, count in collections.Counter(groupID).items() if count > 1]
    coreGroups = os.listdir(coreDir + '/core_orthologs/' + coreSet)
    header = 'genomeID\tsimilar\tdissimilar\tduplicated\tmissing\tignored\ttotal'
    stat = '%s\n%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (header, queryID, type.count('similar'), type.count('dissimilar'), len(dup), len(readFile(missing)), len(readFile(ignored)), len(coreGroups))
    if mode == 4:
        header = 'genomeID\tcomplete\tfragmented\tduplicated\tmissing\tignored\ttotal'
        stat = '%s\n%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (header, queryID, type.count('complete'), type.count('fragmented'), len(dup), len(readFile(missing)), len(readFile(ignored)), len(coreGroups))
    summaryFile.write(stat)
    summaryFile.close()
    return(stat)

def assessCompteness(args):
    coreDir = os.path.abspath(args.coreDir)
    coreSet = args.coreSet
    checkFileExist(coreDir + '/core_orthologs/' + coreSet, '')
    mode = args.mode
    queryID = args.queryID
    outDir = os.path.abspath(args.outDir)
    if not 'fcatOutput' in outDir:
        outDir = outDir + '/fcatOutput'
    ppDir = '%s/%s/%s/phyloprofileOutput' % (outDir, coreSet, queryID)
    checkFileExist(os.path.abspath(ppDir), 'No phylogenetic profile folder found!')
    mode = args.mode
    force = args.force
    # do assessment and print report
    if mode == 0:
        for m in range(1,5):
            print('Mode %s:' % m)
            (noCutoff, stat) = doAssessment(ppDir, coreDir, coreSet, queryID, outDir, m)
            if len(noCutoff) > 0:
                print('\033[92mWARNING: No cutoff for %s group(s):\033[0m\n%s\n' % (len(noCutoff), ','.join(noCutoff)))
            if stat:
                print(stat)
    else:
        (noCutoff, stat) = doAssessment(ppDir, coreDir, coreSet, queryID, outDir, mode)
        if len(noCutoff) > 0:
            print('\033[92mWARNING: No cutoff for %s group(s):\033[0m\n%s\n' % (len(noCutoff), ','.join(noCutoff)))
        if stat:
            print(stat)

    # merge all report (if available)
    modes = ('mode_1', 'mode_2', 'mode_3', 'mode_4')
    mergedStat = open('%s/%s/%s/all_summary.txt' % (outDir, coreSet, queryID), 'w')
    mergedStat.write('mode\tgenomeID\tsimilar\tdissimilar\tduplicated\tmissing\tignored\ttotal\n')
    for m in modes:
        if os.path.exists('%s/%s/%s/%s/summary.txt' % (outDir, coreSet, queryID, m)):
            for line in readFile('%s/%s/%s/%s/summary.txt' % (outDir, coreSet, queryID, m)):
                if not line.split('\t')[0] == 'genomeID':
                    mergedStat.write('%s\t%s\n' % (m, line.strip()))
    mergedStat.close()

def main():
    version = '0.0.1'
    parser = argparse.ArgumentParser(description='You are running assessCompleteness version ' + str(version) + '.')
    required = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('optional arguments')
    required.add_argument('-d', '--coreDir', help='Path to core set directory, where folder core_orthologs can be found', action='store', default='', required=True)
    required.add_argument('-c', '--coreSet', help='Name of core set, which is subfolder within coreDir/core_orthologs/ directory', action='store', default='', required=True)
    required.add_argument('-o', '--outDir', help='Path to output directory', action='store', default='')
    required.add_argument('--queryID', help='ID of taxon of interest (e.g. HUMAN@9606@3)', action='store', default='', type=str)
    optional.add_argument('-m', '--mode',
                        help='Score cutoff mode. (0) all modes, (1) mean of all-vs-all FAS scores, (2) mean FAS of refspec seed, (3) lower endpoint of CI of all-vs-all FAS scores, (4) mean and stdev of sequence length, (5) mean consensus, (6) reference consensus',
                        action='store', default=0, choices=[0,1,2,3,4,5,6], type=int)
    optional.add_argument('--force', help='Force overwrite existing data', action='store_true', default=False)
    args = parser.parse_args()

    start = time.time()
    assessCompteness(args)
    ende = time.time()
    print('Finished in ' + '{:5.3f}s'.format(ende-start))

if __name__ == '__main__':
    main()
