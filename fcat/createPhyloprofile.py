# -*- coding: utf-8 -*-

#######################################################################
#  Copyright (C) 2020 Vinh Tran
#
#  Create phylogenetic profile outputs from the ortholog search result
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
import multiprocessing as mp
from tqdm import tqdm
from greedyFAS.mainFAS import fasInput
import fcat.functions as fcatFn

def readRefspecFile(refspecFile):
    groupRefspec = {}
    for line in fcatFn.readFile(refspecFile):
        groupRefspec[line.split('\t')[0]] = line.split('\t')[1]
    return(groupRefspec)

def outputMode(outDir, coreSet, queryID, force, ppFile):
    phyloprofileDir = '%s/fcatOutput/%s/%s/phyloprofileOutput' % (outDir, coreSet, queryID)
    Path(phyloprofileDir).mkdir(parents=True, exist_ok=True)
    if not os.path.exists('%s/%s' % (phyloprofileDir, ppFile)):
        mode = 3
    else:
        if force:
            mode = 1
        else:
            mode = 0
    return(mode, phyloprofileDir)

def createProfile1(coreDir, outDir, coreSet, queryID, force, groupRefspec):
    (mode, phyloprofileDir) = outputMode(outDir, coreSet, queryID, force, queryID+'_mode1.phyloprofile')
    if mode == 1 or mode == 3:
        # output files
        finalFa = open('%s/%s.mod.fa' % (phyloprofileDir, queryID), 'w')
        finalPhyloprofile = open('%s/%s_mode1.phyloprofile' % (phyloprofileDir, queryID), 'w')
        lines_seen = set() # holds lines already seen
        header = 'geneID\tncbiID\torthoID\tFAS\n'
        finalPhyloprofile.write(header)
        lines_seen.add(header)
        finalLen = open('%s/%s_length.phyloprofile' % (phyloprofileDir, queryID), 'w')
        linesLen_seen = set() # holds lines already seen
        finalLen.write('geneID\tncbiID\torthoID\tLength\n')
        linesLen_seen.add(header)
        # parse into phyloprofile files
        fdogOutDir = '%s/fcatOutput/%s/%s/fdogOutput' % (outDir, coreSet, queryID)
        mergedFa = '%s/%s_all.extended.fa' % (fdogOutDir, queryID)
        # get fas scores for each group
        groupScoreFwd = {} # all fwd fas scores of query ortholog vs core proteins
        groupScoreRev = {} # all rev fas scores of query ortholog vs core proteins
        groupOrtho = {} # query ortholog ID of each group, used for orthoID column of phyloprofile output
        for line in fcatFn.readFile('%s/%s_all.phyloprofile' % (fdogOutDir, queryID)):
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
            newline = '%s\t%s\t%s\t%s\n' % (groupIDmod, 'ncbi' + str(queryID.split('@')[1]), groupOrthoMod, fcatFn.roundTo4(statistics.mean((statistics.mean(groupScoreFwd[groupID]), statistics.mean(groupScoreRev[groupID])))))
            if newline not in lines_seen: # not a duplicate
                finalPhyloprofile.write(newline)
                lines_seen.add(newline)
            # append profile of core sequences for this group
            meanCoreFile = '%s/core_orthologs/%s/%s/fas_dir/cutoff_dir/2.cutoff' % (coreDir, coreSet, groupIDmod)
            for tax in fcatFn.readFile(meanCoreFile):
                if not tax.split('\t')[0] == 'taxa':
                    # not include core taxon that have the same taxonomy ID as query
                    if queryID.split('@')[1] == groupRefspec[groupIDmod].split('@')[1]:
                        if not tax.split('\t')[0] == groupRefspec[groupIDmod]:
                            ppCore = '%s\t%s\t%s|1\t%s\n' % (groupIDmod, 'ncbi' + str(tax.split('\t')[0].split('@')[1]), tax.split('\t')[2].strip(), fcatFn.roundTo4(float(tax.split('\t')[1])))
                            if ppCore not in lines_seen: # not a duplicate
                                finalPhyloprofile.write(ppCore)
                                lines_seen.add(ppCore)
                    else:
                        ppCore = '%s\t%s\t%s|1\t%s\n' % (groupIDmod, 'ncbi' + str(tax.split('\t')[0].split('@')[1]), tax.split('\t')[2].strip(), fcatFn.roundTo4(float(tax.split('\t')[1])))
                        if ppCore not in lines_seen: # not a duplicate
                            finalPhyloprofile.write(ppCore)
                            lines_seen.add(ppCore)
        # length phyloprofile file and final fasta file
        checkFa = []
        for s in SeqIO.parse(mergedFa, 'fasta'):
            idMod = '_'.join(s.id.split('_')[1:])
            if queryID.split('@')[1] == groupRefspec[groupIDmod].split('@')[1]:
                if not idMod.split('|')[1] == groupRefspec[idMod.split('|')[0]]:
                    if not idMod in checkFa:
                        finalFa.write('>%s\n%s\n' % (idMod, s.seq))
                        checkFa.append(idMod)
                    ppLen = '%s\t%s\t%s\t%s\n' % (idMod.split('|')[0], 'ncbi' + str(idMod.split('|')[1].split('@')[1]), idMod, len(s.seq))
                    if ppLen not in linesLen_seen: # not a duplicate
                        finalLen.write(ppLen)
                        linesLen_seen.add(ppLen)
            else:
                if not idMod in checkFa:
                    finalFa.write('>%s\n%s\n' % (idMod, s.seq))
                    checkFa.append(idMod)
                ppLen = '%s\t%s\t%s\t%s\n' % (idMod.split('|')[0], 'ncbi' + str(idMod.split('|')[1].split('@')[1]), idMod, len(s.seq))
                if ppLen not in linesLen_seen: # not a duplicate
                    finalLen.write(ppLen)
                    linesLen_seen.add(ppLen)
        # add missing groups
        if os.path.exists('%s/fcatOutput/%s/%s/missing.txt' % (outDir, coreSet, queryID)):
            for missingGr in fcatFn.readFile('%s/fcatOutput/%s/%s/missing.txt' % (outDir, coreSet, queryID)):
                # add to mode2 profile
                meanCoreFile = '%s/core_orthologs/%s/%s/fas_dir/cutoff_dir/2.cutoff' % (coreDir, coreSet, missingGr)
                for tax in fcatFn.readFile(meanCoreFile):
                    if not tax.split('\t')[0] == 'taxa':
                        # not include core taxon that have the same taxonomy ID as query
                        if queryID.split('@')[1] == groupRefspec[missingGr].split('@')[1]:
                            if not tax.split('\t')[0] == groupRefspec[missingGr]:
                                ppCore = '%s\t%s\t%s|1\t%s\n' % (missingGr, 'ncbi' + str(tax.split('\t')[0].split('@')[1]), tax.split('\t')[2], fcatFn.roundTo4(float(tax.split('\t')[1])))
                                if ppCore not in lines_seen: # not a duplicate
                                    finalPhyloprofile.write(ppCore)
                                    lines_seen.add(ppCore)
                        else:
                            ppCore = '%s\t%s\t%s|1\t%s\n' % (missingGr, 'ncbi' + str(tax.split('\t')[0].split('@')[1]), tax.split('\t')[2], fcatFn.roundTo4(float(tax.split('\t')[1])))
                            if ppCore not in lines_seen: # not a duplicate
                                finalPhyloprofile.write(ppCore)
                                lines_seen.add(ppCore)
                # add to length profile and fasta file
                faCoreFile = '%s/core_orthologs/%s/%s/%s.fa' % (coreDir, coreSet, missingGr, missingGr)
                for sm in SeqIO.parse(faCoreFile, 'fasta'):
                    finalFa.write('>%s\n%s\n' % (sm.id+'|1', sm.seq))
                    ppLen = '%s\t%s\t%s\t%s\n' % (sm.id.split('|')[0], 'ncbi' + str(sm.id.split('|')[1].split('@')[1]), sm.id+'|1', len(sm.seq))
                    if ppLen not in linesLen_seen: # not a duplicate
                        finalLen.write(ppLen)
                        linesLen_seen.add(ppLen)
        finalPhyloprofile.close()
        finalFa.close()
        finalLen.close()
        # delete duplicate lines
        # removeDup('%s/%s_mode1.phyloprofile' % (phyloprofileDir, queryID))
        # removeDup('%s/%s_length.phyloprofile' % (phyloprofileDir, queryID))

def createProfile23(coreDir, outDir, coreSet, queryID, force):
    # output files
    (mode, phyloprofileDir) = outputMode(outDir, coreSet, queryID, force, queryID+'_mode2.phyloprofile')
    if not mode == 0:
        finalPhyloprofile = open('%s/%s_mode2.phyloprofile' % (phyloprofileDir, queryID), 'w')
        lines_seen = set() # holds lines already seen
        header = 'geneID\tncbiID\torthoID\tFAS\n'
        finalPhyloprofile.write(header)
        lines_seen.add(header)
        # parse into phyloprofile file
        fdogOutDir = '%s/fcatOutput/%s/%s/fdogOutput' % (outDir, coreSet, queryID)
        out = os.listdir(fdogOutDir)
        for refSpec in out:
            if os.path.isdir(fdogOutDir + '/' + refSpec):
                refDir = fdogOutDir + '/' + refSpec
                groups = os.listdir(refDir)
                # move to phyloprofile output dir
                if os.path.exists('%s/%s.phyloprofile' % (refDir, refSpec)):
                    for line in fcatFn.readFile('%s/%s.phyloprofile' % (refDir, refSpec)):
                        if queryID in line:
                            tmpQuery = line.split('\t')
                            # statistics.mean((float(line.split('\t')[3], float(line.split('\t')[4]))
                            line = '%s\t%s\t%s\t%s\n' % (tmpQuery[0], tmpQuery[1], tmpQuery[2], fcatFn.roundTo4(statistics.mean((float(line.split('\t')[3]), float(line.split('\t')[4])))))
                            if line not in lines_seen: # not a duplicate
                                finalPhyloprofile.write(line)
                                lines_seen.add(line)
                # append profile of core sequences
                for groupID in groups:
                    coreFasDir = '%s/core_orthologs/%s/%s/fas_dir/fasscore_dir' % (coreDir, coreSet, groupID)
                    for fasFile in glob.glob('%s/*.tsv' % coreFasDir):
                        if queryID.split('@')[1] == refSpec.split('@')[1]:
                            if not refSpec in fasFile:
                                for fLine in fcatFn.readFile(fasFile):
                                    if refSpec in fLine.split('\t')[0]:
                                        tmp = fLine.split('\t')
                                        revFAS = 0
                                        revFile = '%s/%s.tsv' % (coreFasDir, tmp[0].split('|')[1])
                                        for revLine in fcatFn.readFile(revFile):
                                            if tmp[1] == revLine.split('\t')[0]:
                                                revFAS = revLine.split('\t')[2].split('/')[0]
                                        coreLine = '%s\t%s\t%s\t%s\n' % (groupID, 'ncbi' + str(tmp[1].split('|')[1].split('@')[1]), tmp[1]+'|1', fcatFn.roundTo4(statistics.mean((float(tmp[2].split('/')[0]), float(revFAS)))))
                                        if coreLine not in lines_seen: # not a duplicate
                                            finalPhyloprofile.write(coreLine)
                                            lines_seen.add(coreLine)
                        else:
                            for fLine in fcatFn.readFile(fasFile):
                                if refSpec in fLine.split('\t')[0]:
                                    tmp = fLine.split('\t')
                                    revFAS = 0
                                    revFile = '%s/%s.tsv' % (coreFasDir, tmp[0].split('|')[1])
                                    for revLine in fcatFn.readFile(revFile):
                                        if tmp[1] == revLine.split('\t')[0]:
                                            revFAS = revLine.split('\t')[2].split('/')[0]
                                    coreLine = '%s\t%s\t%s\t%s\n' % (groupID, 'ncbi' + str(tmp[1].split('|')[1].split('@')[1]), tmp[1]+'|1', fcatFn.roundTo4(statistics.mean((float(tmp[2].split('/')[0]), float(revFAS)))))
                                    if coreLine not in lines_seen: # not a duplicate
                                        finalPhyloprofile.write(coreLine)
                                        lines_seen.add(coreLine)
        # finalize
        finalPhyloprofile.close()
        # delete duplicate lines
        # removeDup('%s/%s_mode2.phyloprofile' % (phyloprofileDir, queryID))
        shutil.copy('%s/%s_mode2.phyloprofile' % (phyloprofileDir, queryID), '%s/%s_mode3.phyloprofile' % (phyloprofileDir, queryID))

def getDomain(args):
    (jsonFile, groupID, seedIDmod, orthoID, orthoIDmod) = args
    proteome = fasInput.read_json(jsonFile)["feature"]
    out = []
    if orthoID in proteome:
        for tool in proteome[orthoID]:
            if not tool == 'length':
                for feature in proteome[orthoID][tool]:
                    for instance in proteome[orthoID][tool][feature]["instance"]:
                        out.append(groupID + "#" + seedIDmod + "\t" + orthoIDmod + "\t" + str(proteome[orthoID]["length"]) +
                                "\t" + feature + "\t" + str(instance[0]) + "\t" + str(instance[1]) + "\tNA\tN")
    else:
        sys.exit('No annotation for %s found in %s' % (orthoID, jsonFile))
    return('\n'.join(out))

def createDomainFile(coreDir, outDir, coreSet, queryID, refspecFile, ppFile, cpus, force):
    # output files
    (mode, phyloprofileDir) = outputMode(outDir, coreSet, queryID, force, queryID+'.domains')
    if not mode == 0:
        # get refspec for each group
        refTaxid = readRefspecFile(refspecFile)
        # get parsing job
        domainJobs = []
        for line in fcatFn.readFile(ppFile):
            if not "geneID" in line:
                groupID = line.split('\t')[2].split('|')[0]
                # get ortho ID
                orthoID = line.split('\t')[2].split('|')[2]
                taxID = line.split('\t')[2].split('|')[1]
                if taxID == queryID:
                    jsonFile = '%s/weight_dir/%s.json' % (coreDir, queryID)
                    domainJobs.append((jsonFile, groupID, line.split('\t')[2], orthoID, line.split('\t')[2]))
                else:
                    jsonFile = '%s/core_orthologs/%s/%s/fas_dir/annotation_dir/%s.json' % (coreDir, coreSet, groupID, groupID)
                    domainJobs.append((jsonFile, groupID, line.split('\t')[2], '%s|%s|%s' % (groupID, taxID, orthoID), line.split('\t')[2]))
                # get prot ID for refspec
                cutoff2File = '%s/core_orthologs/%s/%s/fas_dir/cutoff_dir/2.cutoff' % (coreDir, coreSet, groupID)
                for ll in fcatFn.readFile(cutoff2File):
                    if refTaxid[groupID] in ll:
                        groupJsonFile = '%s/core_orthologs/%s/%s/fas_dir/annotation_dir/%s.json' % (coreDir, coreSet, groupID, groupID)
                        domainJobs.append((groupJsonFile, groupID, line.split('\t')[2], ll.split('\t')[-1], ll.split('\t')[-1]+'|1'))
                        break
        # parse domains
        pool = mp.Pool(cpus)
        domainOut = []
        for _ in tqdm(pool.imap_unordered(getDomain, list(dict.fromkeys(domainJobs))), total=len(list(dict.fromkeys(domainJobs)))):
            domainOut.append(_)
        pool.close()
        pool.join()
        # write domain output
        with open('%s/fcatOutput/%s/%s/phyloprofileOutput/%s.domains' % (outDir, coreSet, queryID, queryID), "w") as finalDomain:
            value = '\n'.join(set(list(filter(None, domainOut))))
            finalDomain.write('%s\n' % value)
        # removeDup('%s/fcatOutput/%s/%s/phyloprofileOutput/%s.domains' % (outDir, coreSet, queryID, queryID))

def deleteFolder(folder):
    if os.path.exists(folder):
        if os.path.isfile(folder):
            os.remove(folder)
        else:
            shutil.rmtree(folder)

def createPhyloProfile(args):
    coreDir = os.path.abspath(args.coreDir)
    coreSet = args.coreSet
    fcatFn.checkFileExist(coreDir + '/core_orthologs/' + coreSet, '')
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
    if not os.path.exists('%s/phyloprofileOutput/%s_mode1.phyloprofile' % (fcatOut, queryID)):
        if not os.path.exists('%s/fdogOutput'):
            if os.path.exists('%s/fdogOutput.tar.gz' % fcatOut):
                shutil.unpack_archive('%s/fdogOutput.tar.gz' % fcatOut, fcatOut + '/', 'gztar')
            else:
                sys.exit('No ortholog output found!')
        if os.path.exists('%s/last_refspec.txt' % fcatOut):
            groupRefspec = readRefspecFile('%s/last_refspec.txt' % fcatOut)
            createProfile1(coreDir, outDir, coreSet, queryID, force, groupRefspec)
            createProfile23(coreDir, outDir, coreSet, queryID, force)
            if args.noDomain == False:
                createDomainFile(coreDir, outDir, coreSet, queryID, '%s/last_refspec.txt' % fcatOut, '%s/phyloprofileOutput/%s_mode1.phyloprofile' % (fcatOut, queryID), args.cpus, force)
        else:
            sys.exit('No last_refspec.txt file found!')

        if keep == False:
            print('Cleaning up...')
            if os.path.exists('%s/fdogOutput/' % (fcatOut)):
                shutil.rmtree('%s/fdogOutput/' % (fcatOut))

def main():
    version = '0.0.28'
    parser = argparse.ArgumentParser(description='You are running fcat version ' + str(version) + '.')
    required = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('optional arguments')
    required.add_argument('-d', '--coreDir', help='Path to core set directory, where folder core_orthologs can be found', action='store', default='', required=True)
    required.add_argument('-c', '--coreSet', help='Name of core set, which is subfolder within coreDir/core_orthologs/ directory', action='store', default='', required=True)
    required.add_argument('--queryID', help='ID of taxon of interest (e.g. HUMAN@9606@3)', action='store', default='', type=str)
    optional.add_argument('-o', '--outDir', help='Path to output directory', action='store', default='')
    optional.add_argument('--forceProfile', help='Force overwrite existing phylogenetic profile data', action='store_true', default=False)
    optional.add_argument('--noDomain', help='Not output feature domain files', action='store_true', default=False)
    optional.add_argument('--keep', help='Keep temporary phyloprofile data', action='store_true', default=False)
    optional.add_argument('--cpus', help='Number of CPUs used for annotation. Default = 4', action='store', default=4, type=int)
    args = parser.parse_args()

    start = time.time()
    createPhyloProfile(args)
    ende = time.time()
    print('Finished in ' + '{:5.3f}s'.format(ende-start))

if __name__ == '__main__':
    main()
