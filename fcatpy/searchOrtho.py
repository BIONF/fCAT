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
import datetime
import statistics
import glob
import tarfile

def checkFileExist(file, msg):
    if not os.path.exists(os.path.abspath(file)):
        sys.exit('%s not found! %s' % (file, msg))

def roundTo4(number):
    return("%.4f" % round(number, 4))

def readFile(file):
    with open(file, 'r') as f:
        lines = f.readlines()
        f.close()
        return(lines)

def make_archive(source, destination, format):
        base = os.path.basename(destination)
        name = base.split('.')[0]
        ext = '.'.join(base.split('.')[1:3] )
        archive_from = os.path.dirname(source)
        archive_to = os.path.basename(source.strip(os.sep))
        shutil.make_archive(name, format, archive_from, archive_to)
        shutil.move('%s.%s' % (name, ext), destination)

def isInt(s):
    try:
        int(s)
        return True
    except ValueError:
        return False

def checkQueryAnno(annoQuery, annoDir):
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
    return(doAnno)

def parseQueryFa(query, taxid, outDir, doAnno, annoDir, cpus):
    queryID = query.split('/')[-1].split('.')[0]
    queryIDtmp = queryID.split('@')
    if not (len(queryIDtmp) == 3 and isInt(queryIDtmp[1])):
        if taxid == '0':
            sys.exit('Query taxon does not have suitable ID format (e.g. HUMAN@9606@3). Please provide its taxonomy ID additionaly using --taxid option!')
        else:
            addTaxon = 'fdog.addTaxon -f %s -i %s -o %s --replace --force' % (query, taxid, outDir)
            if doAnno == False:
                addTaxon = addTaxon + ' --noAnno'
            else:
                print('Annotation for %s not given!' % queryID)
            try:
                addTaxonOut = subprocess.run([addTaxon], shell=True, capture_output=True, check=True)
            except:
                sys.exit('Problem occurred while parsing query fasta file\n%s' % addTaxon)
            lines = addTaxonOut.stdout.decode().split('\n')
            queryID = lines[1].split('\t')[1]
    else:
        Path('%s/genome_dir/%s' % (outDir, queryID)).mkdir(parents=True, exist_ok=True)
        shutil.copy(query, '%s/genome_dir/%s/%s.fa' % (outDir, queryID, queryID))
        checkedFile = open('%s/genome_dir/%s/%s.fa.checked' % (outDir, queryID, queryID), 'w')
        now = datetime.datetime.now()
        checkedFile.write(now.strftime("%Y-%m-%d %H:%M:%S"))
        checkedFile.close()
        if doAnno:
            annoFAS = 'annoFAS -i %s -o %s --cpus %s > /dev/null 2>&1' % (query, annoDir, cpus)
            try:
                subprocess.run([annoFAS], shell=True, check=True)
            except:
                print('\033[91mProblem occurred while running annoFAS for query protein set\033[0m\n%s' % annoFAS)
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

def readRefspecFile(refspecFile):
    groupRefspec = {}
    for line in readFile(refspecFile):
        groupRefspec[line.split('\t')[0]] = line.split('\t')[1].strip()
    return(groupRefspec)

def prepareJob(coreDir, coreSet, queryID, refspecList, outDir, blastDir, annoDir, annoQuery, force, cpus):
    fdogJobs = []
    ignored = []
    groupRefspec = {}
    hmmPath = coreDir + '/core_orthologs/' + coreSet
    groups = os.listdir(hmmPath)
    if len(groups) > 0:
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
    return(fdogJobs, ignored, groupRefspec)

def runFdog(args):
    (seqFile, seqName, refSpec, outPath, blastPath, hmmPath, searchPath, force) = args
    fdog = 'fdog.run --seqFile %s --seqName %s --refspec %s --outpath %s --blastpath %s --hmmpath %s --searchpath %s --fasoff --reuseCore --checkCoorthologsRef --cpu 1 > /dev/null 2>&1' % (seqFile, seqName, refSpec, outPath, blastPath, hmmPath, searchPath)
    if force:
        fdog = fdog + ' --force'
    try:
        subprocess.run([fdog], shell=True, check=True)
        os.remove(seqName + '.fa')
    except:
        print('\033[91mProblem occurred while running fDOG for \'%s\' core group\033[0m\n%s' % (seqName, fdog))

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

def calcFAS(coreDir, outDir, coreSet, queryID, annoDir, cpus, force):
    # output files
    (mode, phyloprofileDir) = outputMode(outDir, coreSet, queryID, force, 'other')
    if mode == 1 or mode == 3:
        finalPhyloprofile = open('%s/mode23.phyloprofile' % (phyloprofileDir), 'w')
        finalPhyloprofile.write('geneID\tncbiID\torthoID\tFAS_F\tFAS_B\n')
    elif mode == 2:
        finalPhyloprofile = open('%s/mode23.phyloprofile' % (phyloprofileDir), 'a')
    # parse single fdog output
    missing = []
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
                            missing.append(groupID)
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
                    for line in readFile('%s/%s.phyloprofile' % (refDir, refSpec)):
                        if queryID in line:
                            finalPhyloprofile.write(line)
                # append profile of core sequences
                for groupID in groups:
                    coreFasDir = '%s/core_orthologs/%s/%s/fas_dir/fasscore_dir' % (coreDir, coreSet, groupID)
                    for fasFile in glob.glob('%s/*.tsv' % coreFasDir):
                        if not refSpec in fasFile:
                            for fLine in readFile(fasFile):
                                if refSpec in fLine.split('\t')[0]:
                                    tmp = fLine.split('\t')
                                    revFAS = 0
                                    revFile = '%s/%s.tsv' % (coreFasDir, tmp[0].split('|')[1])
                                    for revLine in readFile(revFile):
                                        if tmp[1] == revLine.split('\t')[0]:
                                            revFAS = revLine.split('\t')[2].split('/')[0]
                                    coreLine = '%s\t%s\t%s\t%s\t%s\n' % (groupID, 'ncbi' + str(tmp[1].split('|')[1].split('@')[1]), tmp[1], tmp[2].split('/')[0], revFAS)
                                    finalPhyloprofile.write(coreLine)
    if not mode == 0:
        finalPhyloprofile.close()
    return(missing)

def calcFASall(coreDir, outDir, coreSet, queryID, annoDir, cpus, force, groupRefspec):
    # output files
    phyloprofileDir = '%s/fcatOutput/%s/%s/phyloprofileOutput' % (outDir, coreSet, queryID)
    (mode, phyloprofileDir) = outputMode(outDir, coreSet, queryID, force, 'mode1')
    if mode == 1 or mode == 3:
        finalFa = open('%s/%s.mod.fa' % (phyloprofileDir, coreSet), 'w')
        finalFwdDomain = open('%s/FAS_forward.domains' % (phyloprofileDir), 'wb')
        finalPhyloprofile = open('%s/mode1.phyloprofile' % (phyloprofileDir), 'w')
        finalPhyloprofile.write('geneID\tncbiID\torthoID\tFAS_MEAN\n')
        finalLen = open('%s/length.phyloprofile' % (phyloprofileDir), 'w')
        finalLen.write('geneID\tncbiID\torthoID\tLength\n')
    elif mode == 2:
        finalFa = open('%s/%s.mod.fa' % (phyloprofileDir, coreSet), 'a')
        finalFwdDomain = open('%s/FAS_forward.domains' % (phyloprofileDir), 'ab')
        finalPhyloprofile = open('%s/mode1.phyloprofile' % (phyloprofileDir), 'a')
        finalLen = open('%s/length.phyloprofile' % (phyloprofileDir), 'a')
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
                        # delete single fdog out
                        # shutil.rmtree('%s/%s' % (refDir, groupID))
        mergedFaFile.close()
        # calculate fas scores for merged _all.extended.fa using fdogFAS
        fdogFAS = 'fdogFAS -i %s -w %s --cores %s' % (mergedFa, annoDir, cpus)
        try:
            subprocess.run([fdogFAS], shell=True, check=True)
        except:
            print('\033[91mProblem occurred while running fdogFAS for \'%s\'\033[0m\n%s' % (mergedFa, fdogFAS))
    # move to phyloprofile output dir
    if not mode == 0:
        # phyloprofile file
        groupScoreFwd = {}
        groupScoreRev = {}
        groupOrtho = {}
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
        for groupID in groupOrtho:
            # calculate mean fas score for ortholog
            groupIDmod = '_'.join(groupID.split('_')[1:])
            groupOrthoMod = '_'.join(groupOrtho[groupID].split('_')[1:])
            newline = '%s\t%s\t%s\t%s\n' % (groupIDmod, 'ncbi' + str(queryID.split('@')[1]), groupOrthoMod, statistics.mean((statistics.mean(groupScoreFwd[groupID]), statistics.mean(groupScoreRev[groupID]))))
            finalPhyloprofile.write(newline)
            # append profile of core sequences
            meanCoreFile = '%s/core_orthologs/%s/%s/fas_dir/cutoff_dir/2.cutoff' % (coreDir, coreSet, groupIDmod)
            for tax in readFile(meanCoreFile):
                if not tax.split('\t')[0] == 'taxa':
                    if not tax.split('\t')[0] == groupRefspec[groupIDmod]:
                        ppCore = '%s\t%s\t%s|1\t%s\n' % (groupIDmod, 'ncbi' + str(tax.split('\t')[0].split('@')[1]), tax.split('\t')[2].strip(), tax.split('\t')[1])
                        finalPhyloprofile.write(ppCore)
        finalPhyloprofile.close()
        # length phyloprofile file and final fasta file
        for s in SeqIO.parse(mergedFa, 'fasta'):
            idMod = '_'.join(s.id.split('_')[1:])
            if not idMod.split('|')[1] == groupRefspec[idMod.split('|')[0]]:
                finalFa.write('>%s\n%s\n' % (idMod, s.seq))
                ppLen = '%s\t%s\t%s\t%s\n' % (idMod.split('|')[0], 'ncbi' + str(idMod.split('|')[1].split('@')[1]), idMod, len(s.seq))
                finalLen.write(ppLen)
        finalFa.close()
        finalLen.close()
        # join domain files
        shutil.copyfileobj(open('%s/%s_all_forward.domains' % (fdogOutDir, queryID), 'rb'), finalFwdDomain)
        finalFwdDomain.close()
        finalDomain = open('%s/FAS.domains' % (phyloprofileDir), 'w')
        for domains in readFile('%s/FAS_forward.domains' % (phyloprofileDir)):
            tmp = domains.split('\t')
            mGroup = '_'.join(tmp[0].split('#')[0].split('_')[1:])
            mQuery = '_'.join(tmp[0].split('#')[1].split('_')[1:])
            mSeed = '_'.join(tmp[1].split('_')[1:])
            domainLine = '%s\t%s\t%s\t%s\t%s\t%s\tNA\tN\n' % (mGroup+'#'+mQuery, mSeed, tmp[2], tmp[3], tmp[4], tmp[5])
            finalDomain.write(domainLine)
        finalDomain.close()
        os.remove('%s/FAS_forward.domains' % (phyloprofileDir))

def calcFAScmd(args):
    (seed, seedIDs, query, anno, out, name) = args
    if not os.path.exists('%s/%s.tsv' % (out, name)):
        cmd = 'calcFAS -s %s --seed_id %s -q %s -a %s -o %s --cpus 1 -n %s --domain > /dev/null 2>&1' % (seed, seedIDs, query, anno, out, name)
        try:
            subprocess.run([cmd], shell=True, check=True)
        except:
            print('\033[91mProblem occurred while running calcFAS\033[0m\n%s' % (cmd))

def calcFAScons(coreDir, outDir, coreSet, queryID, annoDir, cpus, force):
    # output files
    (mode, phyloprofileDir) = outputMode(outDir, coreSet, queryID, force, 'other')
    if mode == 1 or mode == 3:
        finalPhyloprofile = open('%s/mode4.phyloprofile' % (phyloprofileDir), 'w')
        finalPhyloprofile.write('geneID\tncbiID\torthoID\tFAS\n')
    elif mode == 2:
        finalPhyloprofile = open('%s/mode4.phyloprofile' % (phyloprofileDir), 'a')
    # parse single fdog output
    missing = []
    fdogOutDir = '%s/fcatOutput/%s/%s/fdogOutput' % (outDir, coreSet, queryID)
    calcFASjob = []

    annoDirTmp = '%s/fcatOutput/%s/%s/tmp/anno' % (outDir, coreSet, queryID)
    Path(annoDirTmp).mkdir(parents=True, exist_ok=True)
    fasDirOutTmp = '%s/fcatOutput/%s/%s/tmp/fasOut' % (outDir, coreSet, queryID)
    Path(fasDirOutTmp).mkdir(parents=True, exist_ok=True)

    fdogOutDir = '%s/fcatOutput/%s/%s/fdogOutput' % (outDir, coreSet, queryID)
    out = os.listdir(fdogOutDir)
    missing = []
    for refSpec in out:
        if os.path.isdir(fdogOutDir + '/' + refSpec):
            # make calcFAS job for founded orthologs and consensus seq of each group
            refDir = fdogOutDir + '/' + refSpec
            groups = os.listdir(refDir)
            groupFa = '%s/%s.fa' % (annoDirTmp, refSpec)
            groupFaFile = open(groupFa, 'w')
            for groupID in groups:
                if os.path.isdir(refDir + '/' + groupID):
                    # get seed and query fasta
                    singleFa = '%s/%s/%s.extended.fa' % (refDir, groupID, groupID)
                    if os.path.exists(singleFa):
                        consFa = '%s/core_orthologs/%s/%s/fas_dir/annotation_dir/cons.fa' % (coreDir, coreSet, groupID)
                        consFaLink = '%s/cons_%s.fa' % (annoDirTmp, groupID)
                        checkFileExist(consFa, '')
                        try:
                            os.symlink(consFa, consFaLink)
                        except FileExistsError:
                            os.remove(consFaLink)
                            os.symlink(consFa, consFaLink)
                        seedID = []
                        for s in SeqIO.parse(singleFa, 'fasta'):
                            if queryID in s.id:
                                idTmp = s.id.split('|')
                                seedID.append(idTmp[-2])
                                groupFaFile.write('>%s\n%s\n' % (idTmp[-2], s.seq))
                        # get annotations for seed and query
                        consJson = '%s/core_orthologs/%s/%s/fas_dir/annotation_dir/cons.json' % (coreDir, coreSet, groupID)
                        consJsonLink = '%s/cons_%s.json' % (annoDirTmp, groupID)
                        checkFileExist(consJson, '')
                        try:
                            os.symlink(consJson, consJsonLink)
                        except FileExistsError:
                            os.remove(consJsonLink)
                            os.symlink(consJson, consJsonLink)
                        # tmp fas output
                        calcFASjob.append([groupFa, ' '.join(seedID), consFaLink, annoDirTmp, fasDirOutTmp, groupID])
                    else:
                        missing.append(groupID)
            groupFaFile.close()
            # get annotation for orthologs
            if not os.path.exists('%s/%s.json' % (annoDirTmp, refSpec)):
                extractAnnoCmd = 'annoFAS -i %s -o %s -e -a %s/%s.json -n %s > /dev/null 2>&1' % (groupFa, annoDirTmp, annoDir, queryID, refSpec)
                try:
                    subprocess.run([extractAnnoCmd], shell=True, check=True)
                except:
                    print('\033[91mProblem occurred while running extracting annotation for \'%s\'\033[0m\n%s' % (seedFa, extractAnnoCmd))
    # do FAS calculation
    pool = mp.Pool(cpus)
    calcFASout = []
    for _ in tqdm(pool.imap_unordered(calcFAScmd, calcFASjob), total=len(calcFASjob)):
        calcFASout.append(_)
    # parse fas output into phyloprofile
    for tsv in os.listdir(fasDirOutTmp):
        if os.path.isfile('%s/%s' % (fasDirOutTmp, tsv)):
            for line in readFile('%s/%s' % (fasDirOutTmp, tsv)):
                if not line.split('\t')[0] == 'Seed':
                    groupID = tsv.split('.')[0]
                    ncbiID = 'ncbi' + str(queryID.split('@')[1])
                    orthoID = line.split('\t')[0]
                    fas = roundTo4(float(line.split('\t')[2].split('/')[0]))
                    finalPhyloprofile.write('%s\t%s\t%s\t%s\n' % (groupID, ncbiID, orthoID, fas))
    finalPhyloprofile.close()

def checkResult(fcatOut, force):
    if force:
        if os.path.exists(fcatOut):
            shutil.rmtree(fcatOut)
        return(0)
    else:
        if not os.path.exists('%s/phyloprofileOutput/mode1.phyloprofile' % fcatOut):
            if os.path.exists('%s/fdogOutput.tar.gz' % fcatOut):
                return(1)
            else:
                if os.path.exists(fcatOut):
                    shutil.rmtree(fcatOut)
                return(0)
        else:
            return(2)

def searchOrtho(args):
    coreDir = os.path.abspath(args.coreDir)
    coreSet = args.coreSet
    checkFileExist(coreDir + '/core_orthologs/' + coreSet, '')
    refspecList = str(args.refspecList).split(",")
    if len(refspecList) == 0:
        sys.exit('No refefence species given! Please specify reference taxa using --refspecList option!')
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
    keep = args.keep

    # check annotation of query species and get query ID
    doAnno = checkQueryAnno(annoQuery, annoDir)
    queryID = parseQueryFa(query, taxid, outDir, doAnno, annoDir, cpus)
    if doAnno == False:
        os.rename(annoDir+'/query.json', annoDir+'/'+queryID+'.json')

    # check old output files
    fcatOut = '%s/fcatOutput/%s/%s' % (outDir, coreSet, queryID)
    status = checkResult(fcatOut, force)

    print('Preparing...')
    groupRefspec = {}
    if status == 0:
        (fdogJobs, ignored, groupRefspec) = prepareJob(coreDir, coreSet, queryID, refspecList, outDir, blastDir, annoDir, annoQuery, force, cpus)
        print('Searching orthologs...')
        pool = mp.Pool(cpus)
        fdogOut = []
        for _ in tqdm(pool.imap_unordered(runFdog, fdogJobs), total=len(fdogJobs)):
            fdogOut.append(_)
        # write ignored groups and refspec for each group based on given refspec list
        if len(ignored) > 0:
            # print('\033[92mNo species in %s found in core set(s): %s\033[0m' % (refspecList, ','.join(ignored)))
            ignoredFile = open('%s/fcatOutput/%s/%s/ignored.txt' % (outDir, coreSet, queryID), 'w')
            ignoredFile.write('\n'.join(ignored))
            ignoredFile.write('\n')
            ignoredFile.close()
        if len(groupRefspec) > 0:
            refspecFile = open('%s/fcatOutput/%s/%s/last_refspec.txt' % (outDir, coreSet, queryID), 'w')
            for g in groupRefspec:
                refspecFile.write('%s\t%s\n' % (g, groupRefspec[g]))
            refspecFile.close()
    elif status == 1:
        # untar old fdog output to create phyloprofile files
        shutil.unpack_archive('%s/fdogOutput.tar.gz' % fcatOut, fcatOut + '/', 'gztar')

    if not status == 2:
        if len(groupRefspec) == 0:
            if os.path.exists('%s/last_refspec.txt' % fcatOut):
                groupRefspec = readRefspecFile('%s/last_refspec.txt' % fcatOut)
        print('Calculating pairwise FAS scores between query orthologs and sequences of refspec...')
        missing = calcFAS(coreDir, outDir, coreSet, queryID, annoDir, cpus, force)
        print('Calculating FAS scores between query orthologs and all sequences in each core group...')
        calcFASall(coreDir, outDir, coreSet, queryID, annoDir, cpus, force, groupRefspec)
        print('Calculating FAS scores between query orthologs and consensus sequence in each core group...')
        calcFAScons(coreDir, outDir, coreSet, queryID, annoDir, cpus, force)
        # remove tmp folder
        if os.path.exists('%s/tmp' % fcatOut):
            shutil.rmtree('%s/tmp' % fcatOut)
        # write missing groups
        if len(missing) > 0:
            missingFile = open('%s/fcatOutput/%s/%s/missing.txt' % (outDir, coreSet, queryID), 'w')
            missingFile.write('\n'.join(missing))
            missingFile.write('\n')
            missingFile.close()

    if os.path.exists('%s/fdogOutput' % fcatOut):
        try:
            make_archive('%s/fdogOutput' % fcatOut, '%s/fdogOutput.tar.gz' % fcatOut, 'gztar')
        except:
            print('Cannot archiving fdog output!')

    if keep == False:
        print('Cleaning up...')
        if os.path.exists('%s/genome_dir/' % (outDir)):
            shutil.rmtree('%s/genome_dir/' % (outDir))
        if os.path.exists('%s/fdogOutput/' % (fcatOut)):
            shutil.rmtree('%s/fdogOutput/' % (fcatOut))
    print('Done! Check output in %s' % fcatOut)

def main():
    version = '0.0.1'
    parser = argparse.ArgumentParser(description='You are running searchOrtho version ' + str(version) + '.')
    required = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('optional arguments')
    required.add_argument('-d', '--coreDir', help='Path to core set directory, where folder core_orthologs can be found', action='store', default='', required=True)
    required.add_argument('-c', '--coreSet', help='Name of core set, which is subfolder within coreDir/core_orthologs/ directory', action='store', default='', required=True)
    required.add_argument('-r', '--refspecList', help='List of reference species', action='store', default='')
    required.add_argument('-q', '--querySpecies', help='Path to gene set for species of interest', action='store', default='')
    optional.add_argument('-o', '--outDir', help='Path to output directory', action='store', default='')
    optional.add_argument('-b', '--blastDir', help='Path to BLAST directory of all core species', action='store', default='')
    optional.add_argument('-a', '--annoDir', help='Path to FAS annotation directory', action='store', default='')
    optional.add_argument('--annoQuery', help='Path to FAS annotation for species of interest', action='store', default='')
    optional.add_argument('-i', '--taxid', help='Taxonomy ID of gene set for species of interest', action='store', default=0, type=int)
    optional.add_argument('--cpus', help='Number of CPUs used for annotation. Default = 4', action='store', default=4, type=int)
    optional.add_argument('--force', help='Force overwrite existing data', action='store_true', default=False)
    optional.add_argument('--keep', help='Keep temporary phyloprofile data', action='store_true', default=False)
    args = parser.parse_args()

    start = time.time()
    searchOrtho(args)
    ende = time.time()
    print('Finished in ' + '{:5.3f}s'.format(ende-start))

if __name__ == '__main__':
    main()
