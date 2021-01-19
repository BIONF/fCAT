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
from Bio import AlignIO
from Bio.Align import AlignInfo
import subprocess
import multiprocessing as mp
import shutil
from tqdm import tqdm
import time
from datetime import datetime
import statistics
from scipy import stats
from rpy2.robjects import FloatVector
from rpy2.robjects.packages import importr

def checkFileExist(file):
    if not os.path.exists(os.path.abspath(file)):
        sys.exit('%s not found' % file)

def roundTo4(number):
    return("%.4f" % round(number, 4))

def annoFAS(groupFa, annoDir, cpus, force):
    ### CAN BE IMPROVED!!!
    ### by modify seq IDs in groupFa, then use extract option of annoFAS
    ### and replace mod IDs by original IDs again in the annotaion json file
    annoFAS = 'annoFAS -i %s -o %s --cpus %s > /dev/null 2>&1' % (groupFa, annoDir, cpus)
    if force:
        annoFAS = annoFAS + ' --force'
    try:
        subprocess.run([annoFAS], shell=True, check=True)
    except:
        print('\033[91mProblem occurred while running annoFAS\033[0m\n%s' % annoFAS)

# def getConsensus(alignmentFile, cov):
#     alignment = AlignIO.read(alignmentFile, 'fasta')
#     summary_align = AlignInfo.SummaryInfo(alignment)
#     consensus = summary_align.dumb_consensus(cov) # cov = 0.5 means that only characters that are present in at least 50% sequences will be considered
#     return(consensus)

def prepareJob(coreDir, coreSet, annoDir, blastDir, bidirectional, force, forceCutoff, cpus):
    groups = os.listdir(coreDir + '/core_orthologs/' + coreSet)
    fasJobs = []
    fasJobsCons = []
    groupRefSpec = {}
    gc = 1
    if force or forceCutoff:
        if os.path.exists('%s/core_orthologs/%s/done.txt' % (coreDir, coreSet)):
            os.remove('%s/core_orthologs/%s/done.txt' % (coreDir, coreSet))
            groups = os.listdir(coreDir + '/core_orthologs/' + coreSet)
    if not os.path.exists('%s/core_orthologs/%s/done.txt' % (coreDir, coreSet)):
        if len(groups) > 0:
            for groupID in groups:
                groupRefSpec[groupID] = []
                group = '%s/core_orthologs/%s/%s' % (coreDir, coreSet, groupID)
                if os.path.isdir(group):
                    print('%s/%s \t %s' % (gc, len(groups), groupID))
                    if force:
                        if os.path.exists('%s/fas_dir' % (group)):
                            shutil.rmtree('%s/fas_dir' % (group))
                    groupFa = '%s/%s.fa' % (group, groupID)
                    annoDirTmp = '%s/fas_dir/annotation_dir/' % (group)
                    Path(annoDirTmp).mkdir(parents=True, exist_ok=True)
                    outDir = '%s/fas_dir/fasscore_dir/' % (group)
                    Path(outDir).mkdir(parents=True, exist_ok=True)
                    # check existing cutoff files
                    if forceCutoff:
                        try:
                            shutil.rmtree('%s/fas_dir/cutoff_dir' % (group))
                        except OSError:
                            pass
                    flag = 0
                    if os.path.exists('%s/fas_dir/cutoff_dir/1.cutoff' % (group)):
                        if not os.stat('%s/fas_dir/cutoff_dir/1.cutoff' % (group)).st_size == 0:
                            flag = 1
                    if flag == 0:
                        # do annotation for this group
                        if not os.path.exists('%s/%s.json' % (annoDirTmp, groupID)) or force:
                            annoFAS(groupFa, annoDirTmp, cpus, force)
                        # get annotation for ref genomes and path to ref genomes
                        for s in SeqIO.parse(groupFa, 'fasta'):
                            ref = s.id.split('|')[1]
                            if not os.path.exists('%s/%s.json' % (annoDirTmp, ref)):
                                if os.path.exists('%s/%s.json' % (annoDir, ref)):
                                    src = '%s/%s.json' % (annoDir, ref)
                                    dst = '%s/%s.json' % (annoDirTmp, ref)
                                    if not os.path.exists(dst):
                                        os.symlink(src, dst)
                            refGenome = '%s/%s/%s.fa' % (blastDir, ref, ref)
                            if not os.path.exists(refGenome):
                                if os.path.islink(refGenome):
                                    refGenome = os.path.realpath(refGenome)
                                else:
                                    sys.exit('%s not found!' % refGenome)
                            checkFileExist(refGenome)
                            fasJobs.append([s.id, ref, groupID, groupFa, annoDirTmp, outDir, refGenome, bidirectional, force])
                            groupRefSpec[groupID].append(ref)
                        # ###### consensus approach
                        # # get consensus sequence
                        # groupAln = '%s/%s.aln' % (group, groupID)
                        # consensus = getConsensus(groupAln, 0.5)
                        # consensusFa = '%s/cons.fa' % annoDirTmp
                        # with open(consensusFa, 'w') as cf:
                        #     cf.write('>consensus\n%s\n' % consensus)
                        # # do annotation for consensus sequence
                        # if not os.path.exists('%s/consensus.json' % (annoDirTmp)) or force:
                        #     annoFAS(consensusFa, annoDirTmp, cpus, force)
                        # # add to fasJobsCons
                        # fasJobsCons.append([coreDir, coreSet, groupID, groupFa, consensusFa, annoDirTmp, outDir, force])
                gc = gc + 1
        else:
            sys.exit('No core group found at %s' % (coreDir + '/core_orthologs/' + coreSet))
    return(fasJobs, fasJobsCons, groupRefSpec)

def calcFAS(args):
    (queryID, refSpec, groupID, groupFa, annoDir, outputDir, ref, bidirectional, force) = args
    flag = 0
    if not os.path.exists('%s/%s.tsv' % (outputDir, refSpec)):
        flag = 1
    else:
        if force:
            os.remove('%s/%s.tsv' % (outputDir, refSpec))
            flag = 1
    if flag == 1:
        # calculate fas scores for each sequence vs all
        fasCmd = 'calcFAS -s \"%s\" -q \"%s\" --query_id \"%s\" -a %s -o %s -n %s --domain -r %s -t 10' % (groupFa, groupFa, queryID, annoDir, outputDir, refSpec, ref)
        if bidirectional:
            fasCmd = fasCmd + ' --bidirectional'
        fasCmd = fasCmd + ' > /dev/null 2>&1'
        try:
            subprocess.run([fasCmd], shell=True, check=True)
        except:
            print('\033[91mProblem occurred while running calcFAS\033[0m\n%s' % fasCmd)

def parseFasOut(fasOutDir, refSpecList):
    fasScores = {}
    fasScores['all'] = {}
    for refSpec in refSpecList:
        fasOut = fasOutDir + '/' + refSpec + '.tsv'
        if not os.path.exists(fasOut):
            sys.exit('%s not found! Probably calcFAS could not run correctly. Please check again!' % fasOut)
        else:
            if os.stat(fasOut).st_size == 0:
                sys.exit('%s is empty! Probably calcFAS could not run correctly. Please check again!' % fasOut)
        if not refSpec in fasScores:
            fasScores[refSpec] = {}
            fasScores[refSpec]['score'] = []
        if not refSpec in fasScores['all']:
            fasScores['all'][refSpec] = {}
        with open(fasOut, 'r') as file:
            for l in file.readlines():
                tmp = l.split('\t')
                if refSpec in tmp[1]:
                    if not refSpec in tmp[0]:
                        # get query spec ID
                        querySpec = tmp[0].split('|')[1]
                        if not querySpec in fasScores:
                            fasScores[querySpec] = {}
                            fasScores[querySpec]['score'] = []
                        if not querySpec in fasScores['all'][refSpec]:
                            fasScores['all'][refSpec][querySpec] = []
                        # get scores for refSpec vs others
                        scores = tmp[2].split('/')
                        if scores[1] == 'NA':
                            fasScores[refSpec]['score'].append(float(scores[0]))
                            fasScores[refSpec]['gene'] = tmp[1]
                            fasScores[querySpec]['score'].append(float(scores[0]))
                            fasScores[querySpec]['gene'] = tmp[0]
                            fasScores['all'][refSpec][querySpec].append(float(scores[0]))
                        else:
                            scores = list(map(float, scores))
                            fasScores[refSpec]['score'].append(statistics.mean(scores))
                            fasScores[refSpec]['gene'] = tmp[1]
                            fasScores[querySpec]['score'].append(statistics.mean(scores))
                            fasScores[querySpec]['gene'] = tmp[0]
                            fasScores['all'][refSpec][querySpec].append(statistics.mean(scores))
    return(fasScores)

def getGroupPairs(scoreDict):
    donePair = []
    out = []
    for s in scoreDict:
        for q in scoreDict[s]:
            if not (s+'_'+q in donePair or q+'_'+s in donePair):
                out.append(statistics.mean((scoreDict[s][q][0], scoreDict[q][s][0])))
                donePair.append(s+'_'+q)
    return(out)

# def parseConsFas(args):
#     (coreDir, coreSet, groupID, groupFa, consensusFa, annoDirTmp, outDir, force) = args
#     # calculate fas scores for each sequence (seed) vs consensus (query)
#     fasCmd = 'calcFAS -s \"%s\" -q \"%s\" -a %s -o %s -t 10 --raw --tsv --domain' % (groupFa, consensusFa, annoDirTmp, outDir)
#     try:
#         fasOut = subprocess.run([fasCmd], shell=True, capture_output=True, check=True)
#     except:
#         print('\033[91mProblem occurred while running calcFAS\033[0m\n%s' % fasCmd)
#     cutoffDir = '%s/core_orthologs/%s/%s/fas_dir/cutoff_dir' % (coreDir, coreSet, groupID)
#     Path(cutoffDir).mkdir(parents=True, exist_ok=True)
#     singleOut = open(cutoffDir + '/4.cutoff', 'w')
#     singleOut.write('taxa\tcutoff\tgene\n')
#     groupOut = open(cutoffDir + '/3.cutoff', 'w')
#     groupOut.write('label\tvalue\n')
#     allFas = []
#     # save each score into 3.scores
#     for line in fasOut.stdout.decode().split('\n'):
#         if '#\t' in line:
#             tmp = line.split('\t')
#             singleOut.write('%s\t%s\t%s\n' % (tmp[1].split('|')[1], roundTo4(float(tmp[3])), tmp[1]))
#             allFas.append(float(tmp[3]))
#     # and mean to 1.scores
#     groupOut.write('meanCons\t%s\n' % roundTo4(statistics.mean(allFas)))
#     groupOut.write('medianCons\t%s\n' % roundTo4(statistics.median(allFas)))
#     singleOut.close()
#     groupOut.close()

def calcCutoff(args):
    (coreDir, coreSet, groupRefSpec, groupID) = args
    EnvStats = importr('EnvStats')
    cutoffDir = '%s/core_orthologs/%s/%s/fas_dir/cutoff_dir' % (coreDir, coreSet, groupID)
    Path(cutoffDir).mkdir(parents=True, exist_ok=True)
    singleOut = open(cutoffDir + '/2.cutoff', 'w')
    singleOut.write('taxa\tcutoff\tgene\n')
    groupOut = open(cutoffDir + '/1.cutoff', 'w')
    groupOut.write('label\tvalue\n')

    # parse fas output into cutoffs
    fasOutDir = '%s/core_orthologs/%s/%s/fas_dir/fasscore_dir' % (coreDir, coreSet, groupID)
    fasScores = parseFasOut(fasOutDir, groupRefSpec[groupID])
    # print(groupID)
    # print(groupRefSpec[groupID])
    # print(fasScores)
    for key in fasScores:
        if key == 'all':
            groupPair = getGroupPairs(fasScores[key])
            # print(groupPair)
            if all(v == 0 for v in groupPair):
                groupOut.write('median\t0.0\n')
                groupOut.write('mean\t0.0\n')
                groupOut.write('LCL\t0.0\n')
                groupOut.write('UCL\t0.0\n')
            else:
                tmp = FloatVector(groupPair)
                # print(tmp)
                ci = EnvStats.eexp(tmp, ci = 'TRUE')
                limits = ci.rx2('interval').rx2('limits')
                rateLCL = list(limits.rx2[1])
                rateUCL = list(limits.rx2[2])
                UCL = 1/rateLCL[0]
                LCL = 1/rateUCL[0]
                groupOut.write('median\t%s\n' % roundTo4(statistics.median(groupPair)))
                groupOut.write('mean\t%s\n' % roundTo4(statistics.mean(groupPair)))
                groupOut.write('LCL\t%s\n' % roundTo4(LCL))
                groupOut.write('UCL\t%s\n' % roundTo4(UCL))
        else:
            singleOut.write('%s\t%s\t%s\n' % (key, roundTo4(statistics.mean(fasScores[key]['score'])), fasScores[key]['gene']))
    # get mean and stddev length for each group
    coreTaxa = []
    groupFa = '%s/core_orthologs/%s/%s/%s.fa' % (coreDir, coreSet, groupID, groupID)
    groupLen = []
    for s in SeqIO.parse(groupFa, 'fasta'):
        groupLen.append(len(s.seq))
        taxId = s.id.split('|')[1]
        if not taxId in coreTaxa:
            coreTaxa.append(taxId)
    groupOut.write('meanLen\t%s\n' % round(statistics.mean(groupLen), 2))
    groupOut.write('stdevLen\t%s\n' % round(statistics.stdev(groupLen), 2))
    singleOut.close()
    groupOut.close()
    # return list of taxa for this core group
    return(coreTaxa)

def calcGroupCutoff(args):
    coreDir = os.path.abspath(args.coreDir)
    coreSet = args.coreSet
    checkFileExist(coreDir + '/core_orthologs/' + coreSet)
    annoDir = args.annoDir
    if annoDir == '':
        annoDir = '%s/weight_dir' % coreDir
    annoDir = os.path.abspath(annoDir)
    checkFileExist(annoDir)
    blastDir = args.blastDir
    if blastDir == '':
        blastDir = '%s/blast_dir' % coreDir
    blastDir = os.path.abspath(blastDir)
    checkFileExist(blastDir)
    cpus = args.cpus
    if cpus >= mp.cpu_count():
        cpus = mp.cpu_count()-1
    bidirectional = args.bidirectional
    force = args.forceCutoffFas
    forceCutoff = args.forceCutoff

    print('Preparing...')
    (fasJobs, fasJobsCons, groupRefSpec) = prepareJob(coreDir, coreSet, annoDir, blastDir, bidirectional, force, forceCutoff, cpus)

    print('Calculating fas scores...')
    pool = mp.Pool(cpus)
    if len(fasJobs) > 0:
        fasOut = []
        for _ in tqdm(pool.imap_unordered(calcFAS, fasJobs), total=len(fasJobs)):
            fasOut.append(_)
    # if len(fasJobsCons) > 0:
    #     fasOutCons = []
    #     for _ in tqdm(pool.imap_unordered(parseConsFas, fasJobsCons), total=len(fasJobsCons)):
    #         fasOutCons.append(_)

    if len(groupRefSpec) > 0:
        print('Calculating cutoffs...')
        cutoffJobs = []
        for groupID in groupRefSpec:
            cutoffJobs.append([coreDir, coreSet, groupRefSpec, groupID])
        cutoffOut = []
        if len(cutoffJobs) > 0:
            for _ in tqdm(pool.imap_unordered(calcCutoff, cutoffJobs), total=len(cutoffJobs)):
                cutoffOut = cutoffOut + _ #.append(_)
        if len(cutoffOut) > 0:
            with open('%s/core_orthologs/%s/done.txt' % (coreDir, coreSet), 'w') as f:
                f.write('%s\n' % '\n'.join(list(dict.fromkeys(cutoffOut))))
                f.close()
    pool.close()
    pool.join()

def main():
    version = '0.0.5'
    parser = argparse.ArgumentParser(description='You are running fcat.cutoff version ' + str(version) + '.')
    required = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('optional arguments')
    required.add_argument('-d', '--coreDir', help='Path to core set directory, where folder core_orthologs can be found', action='store', default='', required=True)
    required.add_argument('-c', '--coreSet', help='Name of core set, which is subfolder within coreDir/core_orthologs/ directory', action='store', default='', required=True)
    optional.add_argument('-a', '--annoDir', help='Path to FAS annotation directory', action='store', default='')
    optional.add_argument('-b', '--blastDir', help='Path to BLAST directory of all core species', action='store', default='')
    optional.add_argument('--cpus', help='Number of CPUs used for annotation. Default = 4', action='store', default=4, type=int)
    optional.add_argument('--bidirectional', help=argparse.SUPPRESS, action='store_true', default=False)
    optional.add_argument('--forceCutoffFas', help='Force overwrite existing data (FAS and cutoff)', action='store_true', default=False)
    optional.add_argument('--forceCutoff', help='Force overwrite cutoff data', action='store_true', default=False)
    args = parser.parse_args()

    start = time.time()
    calcGroupCutoff(args)
    ende = time.time()
    print('Finished in ' + '{:5.3f}s'.format(ende-start))

if __name__ == '__main__':
    main()
