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

def checkFileExist(file, msg):
    if not os.path.exists(os.path.abspath(file)):
        sys.exit('%s not found! %s' % (file, msg))

def parseQueryFa(query, taxid, outDir, doAnno):
    queryID = ''
    addTaxon = 'fdog.addTaxon -f %s -i %s -o %s --replace --force' % (query, taxid, outDir)
    if doAnno == False:
        addTaxon = addTaxon + ' --noAnno'
    try:
        addTaxonOut = subprocess.run([addTaxon], shell=True, capture_output=True, check=True)
    except:
        sys.exit('Problem occurred while parsing query fasta file\n%s' % addTaxon)
    lines = addTaxonOut.stdout.decode().split('\n')
    queryID = lines[1].split('\t')[1]
    return(queryID)

def checkRefspec(refspecList, groupFa):
    coreSpec = []
    for s in SeqIO.parse(groupFa, 'fasta'):
        ref = s.id.split('|')[1]
        coreSpec.append(ref)
    for r in refspecList:
        if r in coreSpec:
            return(r)
    return('')

def prepareJob(coreDir, coreSet, query, taxid, refspecList, outDir, blastDir, annoDir, annoQuery, force, cpus):
    fdogJobs = []
    ignored = []
    groupRefspec = {}
    queryID = ''
    hmmPath = coreDir + '/core_orthologs/' + coreSet
    groups = os.listdir(hmmPath)
    if len(groups) > 0:
        # get query spec ID and searchpath
        doAnno = True
        if not annoQuery == '':
            annoQuery = os.path.abspath(annoQuery)
            checkFileExist(annoQuery, '')
            try:
                os.symlink(annoQuery, annoDir+'/query.json')
            except FileExistsError:
                os.remove(annoDir+'/query.json')
                os.symlink(annoQuery, annoDir+'/query.json')
            doAnno = False
        queryID = parseQueryFa(query, taxid, outDir, doAnno)
        if doAnno == False:
            os.rename(annoDir+'/query.json', annoDir+'/'+queryID+'.json')
        searchPath = '%s/genome_dir' % (outDir)
        # create single fdog job for each core group
        for groupID in groups:
            if os.path.isdir(hmmPath + '/' + groupID):
                groupFa = '%s/core_orthologs/%s/%s/%s.fa' % (coreDir, coreSet, groupID, groupID)
                # check refspec
                refspec = checkRefspec(refspecList, groupFa)
                if refspec == '':
                    ignored.append(groupID)
                else:
                    outPath = '%s/fcatOutput/%s/%s/fdogOutput/%s' % (outDir, coreSet, queryID, refspec)
                    if not os.path.exists('%s/%s/%s.phyloprofile' % (outPath, groupID, groupID)) or force:
                        fdogJobs.append([groupFa, groupID, refspec, outPath, blastDir, hmmPath, searchPath, force])
                    groupRefspec[groupID] = refspec
    else:
        sys.exit('No core group found at %s' % (coreDir + '/core_orthologs/' + coreSet))
    return(fdogJobs, ignored, queryID, groupRefspec)

def runFdog(args):
    (seqFile, seqName, refSpec, outPath, blastPath, hmmPath, searchPath, force) = args
    fdog = 'fdog.run --seqFile %s --seqName %s --refspec %s --outpath %s --blastpath %s --hmmpath %s --searchpath %s --fasoff --reuseCore --cpu 1 > /dev/null 2>&1' % (seqFile, seqName, refSpec, outPath, blastPath, hmmPath, searchPath)
    if force:
        fdog = fdog + ' --force'
    try:
        subprocess.run([fdog], shell=True, check=True)
        os.remove(seqName + '.fa')
    except:
        print('\033[91mProblem occurred while running fDOG for \'%s\' core group\033[0m\n%s' % (seqName, fdog))

def outputMode(outDir, coreSet, queryID, force, append, approach):
    phyloprofileDir = '%s/fcatOutput/%s/%s/phyloprofileOutput' % (outDir, coreSet, queryID)
    Path(phyloprofileDir).mkdir(parents=True, exist_ok=True)
    if not os.path.exists('%s/%s_%s.phyloprofile' % (phyloprofileDir, coreSet, approach)):
        mode = 3
    else:
        if force:
            mode = 1
        elif append:
            mode = 2
        else:
            mode = 0
    return(mode, phyloprofileDir)

def calcFAS(outDir, coreSet, queryID, annoDir, cpus, force, append, cleanup):
    # output files
    missingFile = open('%s/fcatOutput/%s/%s/%s_missing.txt' % (outDir, coreSet, queryID, queryID), 'w')
    (mode, phyloprofileDir) = outputMode(outDir, coreSet, queryID, force, append, 'other')
    if mode == 1 or mode == 3:
        finalPhyloprofile = open('%s/%s_other.phyloprofile' % (phyloprofileDir, queryID), 'w')
        finalPhyloprofile.write('geneID\tncbiID\torthoID\tFAS_F\tFAS_B\n')
        finalFwdDomain = open('%s/%s_other_forward.domains' % (phyloprofileDir, queryID), 'wb')
        finalRevDomain = open('%s/%s_other_reverse.domains' % (phyloprofileDir, queryID), 'wb')
    elif mode == 2:
        finalPhyloprofile = open('%s/%s_other.phyloprofile' % (phyloprofileDir, queryID), 'a')
        finalFwdDomain = open('%s/%s_other_forward.domains' % (phyloprofileDir, queryID), 'ab')
        finalRevDomain = open('%s/%s_other_reverse.domains' % (phyloprofileDir, queryID), 'ab')
    # parse single fdog output
    fdogOutDir = '%s/fcatOutput/%s/%s/fdogOutput' % (outDir, coreSet, queryID)
    out = os.listdir(fdogOutDir)
    for refSpec in out:
        if os.path.isdir(fdogOutDir + '/' + refSpec):
            # merge single extended.fa files for each refspec
            refDir = fdogOutDir + '/' + refSpec
            groups = os.listdir(refDir)
            mergedFa = '%s/%s.extended.fa' % (refDir, refSpec)
            if not os.path.exists(mergedFa) or force:
                mergedFaFile = open(mergedFa, 'wb')
                for groupID in groups:
                    if os.path.isdir(refDir + '/' + groupID):
                        singleFa = '%s/%s/%s.extended.fa' % (refDir, groupID, groupID)
                        if os.path.exists(singleFa):
                            shutil.copyfileobj(open(singleFa, 'rb'), mergedFaFile)
                        else:
                            missingFile.write(groupID + '\n')
                mergedFaFile.close()
                # calculate fas scores for merged extended.fa using fdogFAS
                fdogFAS = 'fdogFAS -i %s -w %s --cores %s' % (mergedFa, annoDir, cpus)
                try:
                    subprocess.run([fdogFAS], shell=True, check=True)
                except:
                    print('\033[91mProblem occurred while running fdogFAS for \'%s\'\033[0m\n%s' % (mergedFa, fdogFAS))
            # move to phyloprofile output dir
            if not mode == 0:
                if os.path.exists('%s/%s.phyloprofile' % (refDir, refSpec)):
                    with open('%s/%s.phyloprofile' % (refDir, refSpec), 'r') as f:
                        lines = f.readlines()
                        f.close()
                        for line in lines:
                            if queryID in line:
                                finalPhyloprofile.write(line)
                    shutil.copyfileobj(open('%s/%s_forward.domains' % (refDir, refSpec), 'rb'), finalFwdDomain)
                    shutil.copyfileobj(open('%s/%s_reverse.domains' % (refDir, refSpec), 'rb'), finalRevDomain)
            # remove files
            if cleanup:
                for r in Path(refDir).glob('%s.*' % refSpec):
                    r.unlink()
                for r in Path(refDir).glob('%s*.domains' % refSpec):
                    r.unlink()
    if not mode == 0:
        finalPhyloprofile.close()
        finalFwdDomain.close()
        finalRevDomain.close()
    missingFile.close()

def calcFASall(coreDir, outDir, coreSet, queryID, annoDir, cpus, force, append, cleanup):
    # output files
    phyloprofileDir = '%s/fcatOutput/%s/%s/phyloprofileOutput' % (outDir, coreSet, queryID)
    (mode, phyloprofileDir) = outputMode(outDir, coreSet, queryID, force, append, 'mode1')
    if mode == 1 or mode == 3:
        finalFa = open('%s/%s.mod.fa' % (phyloprofileDir, queryID), 'w')
        finalLen = open('%s/%s_len.phyloprofile' % (phyloprofileDir, queryID), 'w')
        finalLen.write('geneID\tncbiID\torthoID\tFAS_F\tFAS_B\n')
        finalFwdDomain = open('%s/%s_mode1_forward.domains' % (phyloprofileDir, queryID), 'wb')
        finalRevDomain = open('%s/%s_mode1_reverse.domains' % (phyloprofileDir, queryID), 'wb')
        finalPhyloprofile = open('%s/%s_mode1.phyloprofile' % (phyloprofileDir, queryID), 'w')
        finalPhyloprofile.write('geneID\tncbiID\torthoID\tFAS_F\tFAS_B\n')
    elif mode == 2:
        finalFa = open('%s/%s.mod.fa' % (phyloprofileDir, queryID), 'a')
        finalFa = open('%s/%s_len.phyloprofile' % (phyloprofileDir, coreSet), 'a')
        finalFwdDomain = open('%s/%s_mode1_forward.domains' % (phyloprofileDir, queryID), 'ab')
        finalRevDomain = open('%s/%s_mode1_reverse.domains' % (phyloprofileDir, queryID), 'ab')
        finalPhyloprofile = open('%s/%s_mode1.phyloprofile' % (phyloprofileDir, queryID), 'a')
    # create file for fdogFAS
    fdogOutDir = '%s/fcatOutput/%s/%s/fdogOutput' % (outDir, coreSet, queryID)
    mergedFa = '%s/%s_all.extended.fa' % (fdogOutDir, queryID)
    count = {}
    if not os.path.exists(mergedFa) or force:
        mergedFaFile = open(mergedFa, 'w')
        out = os.listdir(fdogOutDir)
        for refSpec in out:
            if os.path.isdir(fdogOutDir + '/' + refSpec):
                refDir = fdogOutDir + '/' + refSpec
                groups = os.listdir(refDir)
                for groupID in groups:
                    if os.path.isdir(refDir + '/' + groupID):
                        # merge each ortholog seq in single extended.fa file with core group fasta file
                        # and write into mergedFaFile
                        groupFa = '%s/core_orthologs/%s/%s/%s.fa' % (coreDir, coreSet, groupID, groupID)
                        singleFa = '%s/%s/%s.extended.fa' % (refDir, groupID, groupID)
                        if os.path.exists(singleFa):
                            for s in SeqIO.parse(singleFa, 'fasta'):
                                specID = s.id.split('|')[1]
                                if specID == queryID:
                                    if not groupID in count:
                                        count[groupID] = 1
                                    else:
                                        count[groupID] = count[groupID]  + 1
                                    id  = str(count[groupID]) + '_' + s.id
                                    mergedFaFile.write('>%s\n%s\n' % (id, s.seq))
                                    for c in SeqIO.parse(groupFa, 'fasta'):
                                        mergedFaFile.write('>%s_%s|1\n%s\n' % (count[groupID], c.id, c.seq))
                                    # also write to final fasta file
                                    if not mode == 0:
                                        finalFa.write('>%s\n%s\n' % (s.id, s.seq))
                                        finalLen.write('%s\t%s\t%s\t%s\n' % (groupID, 'ncbi' + str(queryID.split('@')[1]), s.id, len(s.seq)))
        mergedFaFile.close()
        # calculate fas scores for merged _all.extended.fa using fdogFAS
        fdogFAS = 'fdogFAS -i %s -w %s --cores %s' % (mergedFa, annoDir, cpus)
        try:
            subprocess.run([fdogFAS], shell=True, check=True)
        except:
            print('\033[91mProblem occurred while running fdogFAS for \'%s\'\033[0m\n%s' % (mergedFa, fdogFAS))
    # move to phyloprofile output dir
    if not mode == 0:
        groupScoreFwd = {}
        groupScoreRev = {}
        groupOrtho = {}
        with open('%s/%s_all.phyloprofile' % (fdogOutDir, queryID), 'r') as f:
            lines = f.readlines()
            f.close()
            for line in lines:
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
        for groupID in groupOrtho:
            groupIDmod = '_'.join(groupID.split('_')[1:])
            groupOrthoMod = '_'.join(groupOrtho[groupID].split('_')[1:])
            newline = '%s\t%s\t%s\t%s\t%s\n' % (groupIDmod, 'ncbi' + str(queryID.split('@')[1]), groupOrthoMod, statistics.mean(groupScoreFwd[groupID]), statistics.mean(groupScoreRev[groupID]))
            finalPhyloprofile.write(newline)
        shutil.copyfileobj(open('%s/%s_all_forward.domains' % (fdogOutDir, queryID), 'rb'), finalFwdDomain)
        shutil.copyfileobj(open('%s/%s_all_reverse.domains' % (fdogOutDir, queryID), 'rb'), finalRevDomain)

        finalFa.close()
        finalLen.close()
        finalPhyloprofile.close()
        finalFwdDomain.close()
        finalRevDomain.close()
    # remove files
    if cleanup:
        for r in Path(fdogOutDir).glob('%s_all*' % queryID):
            r.unlink()

def main():
    version = '0.0.1'
    parser = argparse.ArgumentParser(description='You are running searchOrtho version ' + str(version) + '.')
    required = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('optional arguments')
    required.add_argument('-d', '--coreDir', help='Path to core set directory, where folder core_orthologs can be found', action='store', default='', required=True)
    required.add_argument('-c', '--coreSet', help='Name of core set, which is subfolder within coreDir/core_orthologs/ directory', action='store', default='', required=True)
    required.add_argument('-r', '--refspecList', help='List of reference species', action='store', default='')
    required.add_argument('-q', '--querySpecies', help='Path to gene set for species of interest', action='store', default='')
    required.add_argument('-i', '--taxid', help='Taxonomy ID of gene set for species of interest', action='store', default='', required=True, type=int)
    optional.add_argument('-o', '--outDir', help='Path to output directory', action='store', default='')
    optional.add_argument('-b', '--blastDir', help='Path to BLAST directory of all core species', action='store', default='')
    optional.add_argument('-w', '--annoDir', help='Path to FAS annotation directory', action='store', default='')
    optional.add_argument('-a', '--annoQuery', help='Path to FAS annotation for species of interest', action='store', default='')
    optional.add_argument('--cpus', help='Number of CPUs used for annotation. Default = 4', action='store', default=4, type=int)
    optional.add_argument('--force', help='Force overwrite existing data', action='store_true', default=False)
    optional.add_argument('--append', help='Append to existing phyloprofile data', action='store_true', default=False)
    optional.add_argument('--cleanup', help='Delete temporary phyloprofile data', action='store_true', default=False)

    args = parser.parse_args()

    coreDir = os.path.abspath(args.coreDir)
    coreSet = args.coreSet
    checkFileExist(coreDir + '/core_orthologs/' + coreSet, '')
    refspecList = str(args.refspecList).split(",")
    query = args.querySpecies
    checkFileExist(os.path.abspath(query), '')
    query = os.path.abspath(query)
    taxid = str(args.taxid)
    outDir = args.outDir
    if outDir == '':
        outDir = os.getcwd()
    else:
        Path(outDir).mkdir(parents=True, exist_ok=True)
    blastDir = args.blastDir
    if blastDir == '':
        blastDir = '%s/blast_dir' % coreDir
    blastDir = os.path.abspath(blastDir)
    checkFileExist(blastDir, 'Please set path to blastDB using --blastDir option.')
    annoDir = args.annoDir
    if annoDir == '':
        annoDir = '%s/weight_dir' % coreDir
    annoDir = os.path.abspath(annoDir)
    checkFileExist(annoDir, 'Please set path to annotation directory using --annoDir option.')
    annoQuery = args.annoQuery

    cpus = args.cpus
    if cpus >= mp.cpu_count():
        cpus = mp.cpu_count()-1
    force = args.force
    append = args.append
    cleanup = args.cleanup

    start = time.time()
    print('Preparing...')
    (fdogJobs, ignored, queryID, groupRefspec) = prepareJob(coreDir, coreSet, query, taxid, refspecList, outDir, blastDir, annoDir, annoQuery, force, cpus)

    print('Searching orthologs...')
    pool = mp.Pool(cpus)
    fdogOut = []
    for _ in tqdm(pool.imap_unordered(runFdog, fdogJobs), total=len(fdogJobs)):
        fdogOut.append(_)
    if cleanup:
        shutil.rmtree('%s/genome_dir' % (outDir))

    print('Calculating pairwise FAS scores between query orthologs and sequences of refspec...')
    calcFAS(outDir, coreSet, queryID, annoDir, cpus, force, append, cleanup)
    print('Calculating FAS scores between query orthologs and all sequences in each core group...')
    calcFASall(coreDir, outDir, coreSet, queryID, annoDir, cpus, force, append, cleanup)

    if len(ignored) > 0:
        print('\033[92mNo species in %s found in core set(s): %s\033[0m' % (refspecList, ','.join(ignored)))
        ignoredFile = open('%s/fcatOutput/%s/%s/%s_ignored.txt' % (outDir, coreSet, queryID, queryID), 'w')
        ignoredFile.write('\n'.join(ignored))
        ignoredFile.close()

    if len(groupRefspec) > 0:
        refspecFile = open('%s/fcatOutput/%s/%s/%s_refspec.txt' % (outDir, coreSet, queryID, queryID), 'w')
        for g in groupRefspec:
            refspecFile.write('%s\t%s\n' % (g, groupRefspec[g]))
        refspecFile.close()

    ende = time.time()
    print('Finished in ' + '{:5.3f}s'.format(ende-start))

if __name__ == '__main__':
    main()
