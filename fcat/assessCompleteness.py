# -*- coding: utf-8 -*-

#######################################################################
#  Copyright (C) 2020 Vinh Tran
#
#  Assign completeness status for each core gene in the core set
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
import shutil
import time
import statistics
import collections
import fcat.functions as fcatFn

def addToDict(dict, groupID, seqID, type, value, cutoff):
    if not groupID in dict:
        dict[groupID] = '%s\t%s\t%s\t%s\t%s' % (groupID, type, seqID, value, cutoff)
    else:
        if 'duplicated' in dict[groupID]:
            dict[groupID] = '%s\n%s\tduplicated (%s)\t%s\t%s\t%s' % (dict[groupID], groupID, type, seqID, value, cutoff)
        else:
            tmp = dict[groupID].split('\t')
            dict[groupID] = '%s\tduplicated (%s)\t%s\t%s\t%s' % (tmp[0], tmp[1], tmp[2], tmp[3], tmp[4])
            dict[groupID] = '%s\n%s\tduplicated (%s)\t%s\t%s\t%s' % (dict[groupID], groupID, type, seqID, value, cutoff)
    return(dict)

def mode1(ppFile, missingGr, coreDir, coreSet, queryID):
    noCutoff = []
    assessment = {}
    flag = 0 # used for checking final profile file
    # dict used to group genes into similar/dissimilar/missing
    geneCat = {}
    geneCat['similar'] = []
    geneCat['dissimilar'] = []
    geneCat['missing'] = []
    # validate orthologs
    for line in fcatFn.readFile(ppFile):
        if len(line) > 0:
            if len(line.split('\t')) < 5:
                flag = 1
            groupID = line.split('\t')[0]
            if queryID in line.split('\t')[2]:
                meanFas = float(line.split('\t')[3])
                scoreFile = '%s/core_orthologs/%s/%s/fas_dir/cutoff_dir/1.cutoff' % (coreDir, coreSet, groupID)
                if os.path.exists(scoreFile):
                    meanGroup = 0
                    for l in fcatFn.readFile(scoreFile):
                        if l.split('\t')[0] == 'mean':
                            meanGroup = float(l.split('\t')[1])
                    if meanFas >= meanGroup:
                        assessment = addToDict(assessment, groupID, line.split('\t')[2], 'similar', meanFas, meanGroup)
                        geneCat['similar'].append(line+'\t0')
                    else:
                        assessment = addToDict(assessment, groupID, line.split('\t')[2], 'dissimilar', meanFas, meanGroup)
                        geneCat['dissimilar'].append(line+'\t1')
                else:
                    noCutoff.append(groupID)
            else:
                if not groupID == 'geneID':
                    if groupID in missingGr:
                        geneCat['missing'].append(line+'\t0')
                    else:
                        geneCat['similar'].append(line+'\t0')
    if flag == 1:
        # create new pp file
        newPP = '%s.mod' % ppFile
        newPPfile = open(newPP, 'w')
        newPPfile.write('geneID\tncbiID\torthoID\tFAS\tAssessment\n')
        newPPfile.write('%s\n' % '\n'.join(geneCat['dissimilar']))
        newPPfile.write('%s\n' % '\n'.join(geneCat['similar']))
        newPPfile.write('%s\n' % '\n'.join(geneCat['missing']))
        newPPfile.close()
        shutil.move('%s.mod' % ppFile, ppFile)
    return(assessment, noCutoff, flag)

def mode2(ppFile, missingGr, coreDir, coreSet, queryID, outDir):
    noCutoff = []
    assessment = {}
    flag = 0 # used for checking final profile file
    # dict used to group genes into similar/dissimilar/missing
    geneCat = {}
    geneCat['similar'] = []
    geneCat['dissimilar'] = []
    geneCat['missing'] = []
    # get refspec for each group
    groupRefspec = {}
    refspecFile = '%s/%s/%s/last_refspec.txt' % (outDir, coreSet, queryID)
    for g in fcatFn.readFile(refspecFile):
        groupRefspec[g.split('\t')[0]] = g.split('\t')[1]
    # validate orthologs
    for line in fcatFn.readFile(ppFile):
        if len(line) > 0:
            if len(line.split('\t')) < 5:
                flag = 1
            groupID = line.split('\t')[0]
            if queryID in line.split('\t')[2]:
                meanFas = float(line.split('\t')[3]) #statistics.mean((float(line.split('\t')[3]), float(line.split('\t')[4].strip())))
                scoreFile = '%s/core_orthologs/%s/%s/fas_dir/cutoff_dir/2.cutoff' % (coreDir, coreSet, groupID)
                if os.path.exists(scoreFile):
                    meanRefspec = 0
                    for l in fcatFn.readFile(scoreFile):
                        if l.split('\t')[0] == groupRefspec[groupID].strip():
                            meanRefspec = float(l.split('\t')[1].strip())
                    if meanFas >= meanRefspec:
                        assessment = addToDict(assessment, groupID, line.split('\t')[2], 'similar', meanFas, meanRefspec)
                        geneCat['similar'].append(line+'\t0')
                    else:
                        assessment = addToDict(assessment, groupID, line.split('\t')[2], 'dissimilar', meanFas, meanRefspec)
                        geneCat['dissimilar'].append(line+'\t1')
                else:
                    noCutoff.append(groupID)
            else:
                if not groupID == 'geneID':
                    if groupID in missingGr:
                        geneCat['missing'].append(line+'\t0')
                    else:
                        geneCat['similar'].append(line+'\t0')
    if flag == 1:
        # create new pp file
        newPP = '%s.mod' % ppFile
        newPPfile = open(newPP, 'w')
        newPPfile.write('geneID\tncbiID\torthoID\tFAS\tAssessment\n')
        newPPfile.write('%s\n' % '\n'.join(geneCat['dissimilar']))
        newPPfile.write('%s\n' % '\n'.join(geneCat['similar']))
        newPPfile.write('%s\n' % '\n'.join(geneCat['missing']))
        newPPfile.close()
        shutil.move('%s.mod' % ppFile, ppFile)
    return(assessment, noCutoff, flag)

def mode3(ppFile, missingGr, coreDir, coreSet, queryID):
    noCutoff = []
    assessment = {}
    flag = 0 # used for checking final profile file
    # dict used to group genes into similar/dissimilar/missing
    geneCat = {}
    geneCat['similar'] = []
    geneCat['dissimilar'] = []
    geneCat['missing'] = []
    # validate orthologs
    for line in fcatFn.readFile(ppFile):
        if len(line) > 0:
            if len(line.split('\t')) < 5:
                flag = 1
            groupID = line.split('\t')[0]
            if queryID in line.split('\t')[2]:
                meanFas = float(line.split('\t')[3]) #statistics.mean((float(line.split('\t')[3]), float(line.split('\t')[4].strip())))
                scoreFile = '%s/core_orthologs/%s/%s/fas_dir/cutoff_dir/1.cutoff' % (coreDir, coreSet, groupID)
                if os.path.exists(scoreFile):
                    LCL = 0
                    UCL = 0
                    for l in fcatFn.readFile(scoreFile):
                        if l.split('\t')[0] == 'LCL':
                            LCL = float(l.split('\t')[1])
                        if l.split('\t')[0] == 'UCL':
                            UCL = float(l.split('\t')[1])
                    # if LCL <= meanFas <= UCL:
                    if LCL <= meanFas:
                        assessment = addToDict(assessment, groupID, line.split('\t')[2], 'similar', meanFas, LCL)
                        geneCat['similar'].append(line+'\t0')
                    else:
                        assessment = addToDict(assessment, groupID, line.split('\t')[2], 'dissimilar', meanFas, LCL)
                        geneCat['dissimilar'].append(line+'\t1')
                else:
                    noCutoff.append(groupID)
            else:
                if not groupID == 'geneID':
                    if groupID in missingGr:
                        geneCat['missing'].append(line+'\t0')
                    else:
                        geneCat['similar'].append(line+'\t0')
    if flag == 1:
        # create new pp file
        newPP = '%s.mod' % ppFile
        newPPfile = open(newPP, 'w')
        newPPfile.write('geneID\tncbiID\torthoID\tFAS\tAssessment\n')
        newPPfile.write('%s\n' % '\n'.join(geneCat['dissimilar']))
        newPPfile.write('%s\n' % '\n'.join(geneCat['similar']))
        newPPfile.write('%s\n' % '\n'.join(geneCat['missing']))
        newPPfile.close()
        shutil.move('%s.mod' % ppFile, ppFile)
    return(assessment, noCutoff, flag)

def mode4(ppFile, missingGr, coreDir, coreSet, queryID):
    noCutoff = []
    assessment = {}
    flag = 0 # used for checking final profile file
    # dict used to group genes into similar/dissimilar/missing
    geneCat = {}
    geneCat['similar'] = []
    geneCat['dissimilar'] = []
    geneCat['missing'] = []
    # validate orthologs
    for line in fcatFn.readFile(ppFile):
        if len(line) > 0:
            if len(line.split('\t')) < 5:
                flag = 1
            if len(line) > 0:
                groupID = line.split('\t')[0]
                if queryID in line.split('\t')[2]:
                    length = float(line.split('\t')[3])
                    scoreFile = '%s/core_orthologs/%s/%s/fas_dir/cutoff_dir/1.cutoff' % (coreDir, coreSet, groupID)
                    if os.path.exists(scoreFile):
                        meanLen = 0
                        stdevLen = 0
                        for l in fcatFn.readFile(scoreFile):
                            if l.split('\t')[0] == 'meanLen':
                                meanLen = float(l.split('\t')[1])
                            if l.split('\t')[0] == 'stdevLen':
                                stdevLen = float(l.split('\t')[1])
                        ppValue = abs(length - meanLen) / max(length, meanLen)
                        if stdevLen > 0:
                            check = abs((length - meanLen) / (2 * stdevLen))
                            if check <= 1:
                                assessment = addToDict(assessment, groupID, line.split('\t')[2], 'complete', length, '%s (sd=%s)' % (fcatFn.roundTo4(meanLen), fcatFn.roundTo4(stdevLen)))
                                # geneCat['similar'].append(line+'\t0')
                                geneCat['similar'].append('\t'.join(line.split('\t')[0:-1])+'\t1')
                            else:
                                assessment = addToDict(assessment, groupID, line.split('\t')[2], 'fragmented', length, '%s (sd=%s)' % (fcatFn.roundTo4(meanLen), fcatFn.roundTo4(stdevLen)))
                                # geneCat['dissimilar'].append(line+'\t1')
                                geneCat['similar'].append('\t'.join(line.split('\t')[0:-1])+'\t0')
                        else:
                            if (length - meanLen) >= 0:
                                assessment = addToDict(assessment, groupID, line.split('\t')[2], 'complete', length, '%s (sd=%s)' % (fcatFn.roundTo4(meanLen), fcatFn.roundTo4(stdevLen)))
                                # geneCat['similar'].append(line+'\t0')
                                geneCat['similar'].append('\t'.join(line.split('\t')[0:-1])+'\t1')
                            else:
                                assessment = addToDict(assessment, groupID, line.split('\t')[2], 'fragmented', length, '%s (sd=%s)' % (fcatFn.roundTo4(meanLen), fcatFn.roundTo4(stdevLen)))
                                # geneCat['dissimilar'].append(line+'\t1')
                                geneCat['similar'].append('\t'.join(line.split('\t')[0:-1])+'\t0')
                    else:
                        noCutoff.append(groupID)
                else:
                    if not groupID == 'geneID':
                        if groupID in missingGr:
                            # geneCat['missing'].append(line+'\t0')
                            geneCat['similar'].append('\t'.join(line.split('\t')[0:-1])+'\t0.5')
                        else:
                            # geneCat['similar'].append(line+'\t0')
                            geneCat['similar'].append('\t'.join(line.split('\t')[0:-1])+'\t1')
    if flag == 1:
        # create new pp file
        newPP = '%s.mod' % ppFile
        newPPfile = open(newPP, 'w')
        newPPfile.write('geneID\tncbiID\torthoID\tAssessment\n')
        newPPfile.write('%s\n' % '\n'.join(geneCat['dissimilar']))
        newPPfile.write('%s\n' % '\n'.join(geneCat['similar']))
        newPPfile.write('%s\n' % '\n'.join(geneCat['missing']))
        newPPfile.close()
        shutil.move('%s.mod' % ppFile, ppFile)
    return(assessment, noCutoff, flag)

# def mode5(ppFile, coreDir, coreSet, queryID):
#     noCutoff = []
#     assessment = {}
#     for line in fcatFn.readFile(ppFile):
#         groupID = line.split('\t')[0]
#         # if queryID in line.split('\t')[2]:
#         if line.split('\t')[1] == 'ncbi'+str(queryID.split('@')[1]):
#             meanFas = float(line.split('\t')[3])
#             scoreFile = '%s/core_orthologs/%s/%s/fas_dir/cutoff_dir/3.cutoff' % (coreDir, coreSet, groupID)
#             if os.path.exists(scoreFile):
#                 meanGroup = 0
#                 for l in fcatFn.readFile(scoreFile):
#                     if l.split('\t')[0] == 'meanCons':
#                         meanGroup = float(l.split('\t')[1])
#                 if meanFas >= meanGroup:
#                     assessment = addToDict(assessment, groupID, line.split('\t')[2], 'similar')
#                 else:
#                     assessment = addToDict(assessment, groupID, line.split('\t')[2], 'dissimilar')
#             else:
#                 noCutoff.append(groupID)
#     return(assessment, noCutoff)
#
# def mode6(ppFile, coreDir, coreSet, queryID, outDir):
#     noCutoff = []
#     assessment = {}
#     # get refspec for each group
#     groupRefspec = {}
#     refspecFile = '%s/%s/%s/last_refspec.txt' % (outDir, coreSet, queryID)
#     for g in fcatFn.readFile(refspecFile):
#         groupRefspec[g.split('\t')[0]] = g.split('\t')[1]
#     # do assessment
#     for line in fcatFn.readFile(ppFile):
#         groupID = line.split('\t')[0]
#         # if queryID in line.split('\t')[2]:
#         if line.split('\t')[1] == 'ncbi'+str(queryID.split('@')[1]):
#             # meanFas = statistics.mean((float(line.split('\t')[3]), float(line.split('\t')[4].strip())))
#             meanFas = float(line.split('\t')[3])
#             scoreFile = '%s/core_orthologs/%s/%s/fas_dir/cutoff_dir/4.cutoff' % (coreDir, coreSet, groupID)
#             if os.path.exists(scoreFile):
#                 meanRefspec = 0
#                 for l in fcatFn.readFile(scoreFile):
#                     if l.split('\t')[0] == groupRefspec[groupID].strip():
#                         meanRefspec = float(l.split('\t')[1].strip())
#                 if meanFas >= meanRefspec:
#                     assessment = addToDict(assessment, groupID, line.split('\t')[2], 'similar')
#                 else:
#                     assessment = addToDict(assessment, groupID, line.split('\t')[2], 'dissimilar')
#             else:
#                 noCutoff.append(groupID)
#     return(assessment, noCutoff)

def writeReport(assessment, outDir, coreDir, coreSet, queryID, mode):
    missing = '%s/%s/%s/missing.txt' % (outDir, coreSet, queryID)
    ignored = '%s/%s/%s/ignored.txt' % (outDir, coreSet, queryID)
    Path('%s/%s/%s/mode_%s' % (outDir, coreSet, queryID, mode)).mkdir(parents=True, exist_ok=True)
    # write full report
    fullFile = open('%s/%s/%s/mode_%s/full.txt' % (outDir, coreSet, queryID, mode), 'w')
    fullFile.write('GroupID\tAssessment\tOrthoID\tValue\tCutoff\n')
    for group in assessment:
        fullFile.write('%s\n' % assessment[group])
    if os.path.exists(missing):
        for m in fcatFn.readFile(missing):
            fullFile.write(m + '\tmissing\n')
    if os.path.exists(ignored):
        for i in fcatFn.readFile(ignored):
            fullFile.write(i + '\tignored\n')
    fullFile.close()
    # write summary report
    summaryFile = open('%s/%s/%s/mode_%s/summary.txt' % (outDir, coreSet, queryID, mode), 'w')
    type = [x.split('\t')[1] for x in open('%s/%s/%s/mode_%s/full.txt' % (outDir, coreSet, queryID, mode)).readlines()]
    groupID = [x.split('\t')[0] for x in open('%s/%s/%s/mode_%s/full.txt' % (outDir, coreSet, queryID, mode)).readlines()]
    dup = [item for item, count in collections.Counter(groupID).items() if count > 1]
    coreGroups = os.listdir(coreDir + '/core_orthologs/' + coreSet)
    header = 'genomeID\tsimilar\tdissimilar\tduplicated\tmissing\tignored\ttotal'
    stat = '%s\n%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (header, queryID, type.count('similar'), type.count('dissimilar'), len(dup), len(fcatFn.readFile(missing)), len(fcatFn.readFile(ignored)), len(coreGroups)-1)
    if mode == 4:
        header = 'genomeID\tcomplete\tfragmented\tduplicated\tmissing\tignored\ttotal'
        stat = '%s\n%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (header, queryID, type.count('complete'), type.count('fragmented'), len(dup), len(fcatFn.readFile(missing)), len(fcatFn.readFile(ignored)), len(coreGroups)-1)
    summaryFile.write(stat)
    summaryFile.close()
    return(stat)

def doAssessment(ppDir, coreDir, coreSet, queryID, outDir, mode):
    assessment = {}
    noCutoff = []
    flag = 0 # flag=0: no new ppFile needed
    # assess completeness
    missingFile = '%s/%s/%s/missing.txt' % (outDir, coreSet, queryID)
    missingGr = []
    if os.path.exists(missingFile):
        with open(missingFile) as f:
            missingGr = [line.rstrip('\n') for line in f]
    if mode == 1:
        ppFile = '%s/%s_mode1.phyloprofile' % (ppDir, queryID)
        (assessment, noCutoff, flag) = mode1(ppFile, missingGr, coreDir, coreSet, queryID)
    elif mode == 2:
        ppFile = '%s/%s_mode2.phyloprofile' % (ppDir, queryID)
        (assessment, noCutoff, flag) = mode2(ppFile, missingGr, coreDir, coreSet, queryID, outDir)
    elif mode == 3:
        ppFile = '%s/%s_mode3.phyloprofile' % (ppDir, queryID)
        (assessment, noCutoff, flag) = mode3(ppFile, missingGr, coreDir, coreSet, queryID)
    elif mode == 4:
        ppFile = '%s/%s_length.phyloprofile' % (ppDir, queryID)
        (assessment, noCutoff, flag) = mode4(ppFile, missingGr, coreDir, coreSet, queryID)
    # elif mode == 5:
    #     ppFile = '%s/mode4.phyloprofile' % (ppDir)
    #     (assessment, noCutoff) = mode5(ppFile, coreDir, coreSet, queryID)
    # elif mode == 6:
    #     ppFile = '%s/mode4.phyloprofile' % (ppDir)
    #     (assessment, noCutoff) = mode6(ppFile, coreDir, coreSet, queryID, outDir)
    # print full report
    stat = writeReport(assessment, outDir, coreDir, coreSet, queryID, mode)
    # remove empty lines in pp file
    fcatFn.removeDupEmpty(ppFile)
    return(noCutoff, stat, flag)

def mergeReports(args):
    (outDir, coreSet, queryID) = args
    modes = ('mode_1', 'mode_2', 'mode_3', 'mode_4')
    mergedStat = open('%s/%s/%s/report_summary.txt' % (outDir, coreSet, queryID), 'w')
    mergedStat.write('mode\tgenomeID\tsimilar\tdissimilar\tduplicated\tmissing\tignored\ttotal\n')
    mergedFull = open('%s/%s/%s/report_full.txt' % (outDir, coreSet, queryID), 'w')
    mergedDismiss = open('%s/%s/%s/report_dismiss.txt' % (outDir, coreSet, queryID), 'w')
    fullTypeDict = {}
    dismissTypeDict = {}
    availableMode = []
    for m in modes:
        availableMode.append('%s_cat\t%s_value\t%s_cutoff\t%s_delta' % (m,m,m,m))
        if os.path.exists('%s/%s/%s/%s/summary.txt' % (outDir, coreSet, queryID, m)):
            for line in fcatFn.readFile('%s/%s/%s/%s/summary.txt' % (outDir, coreSet, queryID, m)):
                if not line.split('\t')[0] == 'genomeID':
                    mergedStat.write('%s\t%s\n' % (m, line))
        if os.path.exists('%s/%s/%s/%s/full.txt' % (outDir, coreSet, queryID, m)):
            for line in fcatFn.readFile('%s/%s/%s/%s/full.txt' % (outDir, coreSet, queryID, m)):
                tmp = line.split('\t') # 1488235at2759	dissimilar	1488235at2759|TRINE@6336@1|T07_14825.1|1	0.5566	0.7034
                if not tmp[0] == 'GroupID':
                    orthoID = 'NA'
                    value = 'NA'
                    cutoff = 'NA'
                    delta = 'NA'
                    if len(tmp) == 5:
                        orthoID = tmp[2]
                        value = tmp[-2]
                        cutoff = tmp[-1]
                        delta = abs(float(value) - float(cutoff.split(' ')[0]))
                    key = '%s\t%s' % (tmp[0], orthoID)
                    if not key in fullTypeDict:
                        fullTypeDict[key] = []
                    fullTypeDict[key].append('%s\t%s\t%s\t%s' % (tmp[1], value, cutoff, delta))
                    if tmp[1] == 'missing' or tmp[1] == 'dissimilar' or tmp[1] == 'fragmented':
                        dismissTypeDict[key] = 1

    mergedFull.write('groupID\torthoID\t%s\n' % '\t'.join(availableMode))
    mergedDismiss.write('groupID\torthoID\t%s\n' % '\t'.join(availableMode))
    for gid in fullTypeDict:
        mergedFull.write('%s\t%s\n' % (gid, '\t'.join(fullTypeDict[gid])))
        if gid in dismissTypeDict:
            mergedDismiss.write('%s\t%s\n' % (gid, '\t'.join(fullTypeDict[gid])))
    mergedStat.close()
    mergedFull.close()
    # delete report for singe mode
    for m in modes:
        if os.path.exists('%s/%s/%s/%s/summary.txt' % (outDir, coreSet, queryID, m)):
            shutil.rmtree('%s/%s/%s/%s' % (outDir, coreSet, queryID, m))

def assessCompteness(args):
    coreDir = os.path.abspath(args.coreDir)
    coreSet = args.coreSet
    fcatFn.checkFileExist(coreDir + '/core_orthologs/' + coreSet, '')
    mode = args.mode
    queryID = args.queryID
    outDir = os.path.abspath(args.outDir)
    if not 'fcatOutput' in outDir:
        outDir = outDir + '/fcatOutput'
    ppDir = '%s/%s/%s/phyloprofileOutput' % (outDir, coreSet, queryID)
    fcatFn.checkFileExist(os.path.abspath(ppDir), 'No phylogenetic profile folder found!')
    mode = args.mode
    force = args.force
    # do assessment and print report
    flag = 0 # flag=0: no new ppFile needed
    if mode == 0:
        for m in range(1,5):
            print('Mode %s:' % m)
            (noCutoff, stat, flag) = doAssessment(ppDir, coreDir, coreSet, queryID, outDir, m)
            if len(noCutoff) > 0:
                print('\033[92mWARNING: No cutoff for %s group(s):\033[0m\n%s\n' % (len(noCutoff), ','.join(noCutoff)))
            if stat:
                print(stat)
    else:
        (noCutoff, stat, flag) = doAssessment(ppDir, coreDir, coreSet, queryID, outDir, mode)
        if len(noCutoff) > 0:
            print('\033[92mWARNING: No cutoff for %s group(s):\033[0m\n%s\n' % (len(noCutoff), ','.join(noCutoff)))
        if stat:
            print(stat)
    # merge all report (if available)
    mergeReports((outDir, coreSet, queryID))
    return(flag)

def main():
    version = '0.0.28'
    parser = argparse.ArgumentParser(description='You are running fcat version ' + str(version) + '.')
    required = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('optional arguments')
    required.add_argument('-d', '--coreDir', help='Path to core set directory, where folder core_orthologs can be found', action='store', default='', required=True)
    required.add_argument('-c', '--coreSet', help='Name of core set, which is subfolder within coreDir/core_orthologs/ directory', action='store', default='', required=True)
    required.add_argument('-o', '--outDir', help='Path to output directory', action='store', default='', required=True)
    required.add_argument('--queryID', help='ID of taxon of interest (e.g. HUMAN@9606@3)', action='store', default='', type=str, required=True)
    optional.add_argument('-m', '--mode',
                        help='Score cutoff mode. (0) all modes, (1) mean of all-vs-all FAS scores, (2) mean FAS of refspec seed, (3) lower endpoint of CI of all-vs-all FAS scores, (4) mean and stdev of sequence length',
                        action='store', default=0, choices=[0,1,2,3,4], type=int)
    optional.add_argument('--force', help='Force overwrite existing data', action='store_true', default=False)
    args = parser.parse_args()

    start = time.time()
    flag = assessCompteness(args)
    ende = time.time()
    print('Finished in ' + '{:5.3f}s'.format(ende-start))

if __name__ == '__main__':
    main()
