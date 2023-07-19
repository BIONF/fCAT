# -*- coding: utf-8 -*-

#######################################################################
#  Copyright (C) 2022 Vinh Tran
#
#  Search orthologs for each core gene in the query species
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
import glob
import argparse
from pathlib import Path
from Bio import SeqIO
import subprocess
import multiprocessing as mp
import shutil
from tqdm import tqdm
import time
import datetime
from pkg_resources import get_distribution
import fcat.functions as fcatFn

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

def checkQueryAnno(annoQuery, annoDir, taxid, query):
    doAnno = True
    queryID = query.split('/')[-1].split('.')[0]
    queryIDtmp = queryID.split('@')
    if not (len(queryIDtmp) == 3 and isInt(queryIDtmp[1])):
        id = taxid
    else:
        id = queryIDtmp[1]
    if not annoQuery == '':
        annoQuery = os.path.abspath(annoQuery)
        fcatFn.checkFileExist(annoQuery, '')
        if not annoDir in annoQuery:
            try:
                os.symlink(annoQuery, '%s/query_%s.json' % (annoDir, id))
            except FileExistsError:
                os.remove( '%s/query_%s.json' % (annoDir, id))
                os.symlink(annoQuery,  '%s/query_%s.json' % (annoDir, id))
        doAnno = False
    return(doAnno, id)

def parseQueryFa(coreSet, query, annoQuery, taxid, outDir, doAnno, annoDir, cpus):
    queryID = query.split('/')[-1].split('.')[0]
    queryIDtmp = queryID.split('@')
    if not (len(queryIDtmp) == 3 and isInt(queryIDtmp[1])):
        if taxid == '0':
            sys.exit('Query taxon does not have suitable ID format (e.g. HUMAN@9606@3). Please provide its taxonomy ID additionaly using --taxid option!')
        else:
            addTaxon = 'fdog.addTaxon -f %s -i %s -o %s --replace' % (query, taxid, outDir)
            if doAnno == False:
                addTaxon = addTaxon + ' --noAnno'
            # else:
            #     if annoQuery == '':
            #         print('Annotation for %s not given! It will take a while for annotating...' % queryID)
            try:
                addTaxonOut = subprocess.run([addTaxon], shell=True, capture_output=True, check=True)
            except:
                sys.exit('Problem occurred while parsing query fasta file\n%s' % addTaxon)
            queryID = ''
            for line in addTaxonOut.stdout.decode().split('\n'):
                if "Species name" in line:
                    queryID = line.split('\t')[1]
                    # print('Query ID used by fCAT and fDOG: %s' % queryID)
            if len(queryID) == 0:
                sys.exit('Cannot identidy queryID!')
    else:
        Path('%s/fcatOutput/%s/%s/searchTaxa_dir/%s' % (outDir, coreSet, queryID, queryID)).mkdir(parents=True, exist_ok=True)
        shutil.copy(query, '%s/fcatOutput/%s/%s/searchTaxa_dir/%s/%s.fa' % (outDir, coreSet, queryID, queryID, queryID))
        os.chmod('%s/fcatOutput/%s/%s/searchTaxa_dir/%s/%s.fa' % (outDir, coreSet, queryID, queryID, queryID), 0o755)
        checkedFile = open('%s/fcatOutput/%s/%s/searchTaxa_dir/%s/%s.fa.checked' % (outDir, coreSet, queryID, queryID, queryID), 'w')
        now = datetime.datetime.now()
        checkedFile.write(now.strftime("%Y-%m-%d %H:%M:%S"))
        checkedFile.close()
        if doAnno:
            annoFAS = 'fas.doAnno -i %s -o %s --cpus %s > /dev/null 2>&1' % (query, outDir, cpus)
            try:
                subprocess.run([annoFAS], shell=True, check=True)
            except:
                print('\033[91mProblem occurred while running annoFAS for query protein set\033[0m\n%s' % annoFAS)
    # make link to new created annotation file (if needed)
    if doAnno:
        try:
            os.symlink('%s/annotation_dir/%s.json' % (outDir, queryID), '%s/%s.json' % (annoDir, queryID))
        except FileExistsError:
            os.remove( '%s/%s.json' % (annoDir, queryID))
            os.symlink('%s/annotation_dir/%s.json' % (outDir, queryID),  '%s/%s.json' % (annoDir, queryID))
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
    for line in fcatFn.readFile(refspecFile):
        groupRefspec[line.split('\t')[0]] = line.split('\t')[1]
    return(groupRefspec)

def prepareJob(coreDir, coreSet, queryID, refspecList, outDir, coreTaxa_dir, annoDir, annoQuery, force, cpus):
    fdogJobs = []
    ignored = []
    groupRefspec = {}
    hmmPath = coreDir + '/core_orthologs/' + coreSet
    groups = os.listdir(hmmPath)
    if len(groups) > 0:
        searchTaxa_path = '%s/fcatOutput/%s/%s/searchTaxa_dir' % (outDir, coreSet, queryID)
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
                        fdogJobs.append([groupFa, groupID, refspec, outPath, coreTaxa_dir, hmmPath, searchTaxa_path, force])
                    groupRefspec[groupID] = refspec
    else:
        sys.exit('No core group found at %s' % (coreDir + '/core_orthologs/' + coreSet))
    return(fdogJobs, ignored, groupRefspec)

def runFdog(args):
    (seqFile, jobName, refSpec, outPath, coreTaxa_path, hmmPath, searchTaxa_path, force) = args
    fdog = 'fdog.run --seqFile %s --jobName %s --refspec %s --outpath %s --corepath %s --hmmpath %s --searchpath %s --fasOff --reuseCore --cpu 1' % (seqFile, jobName, refSpec, outPath, coreTaxa_path, hmmPath, searchTaxa_path)
    if force:
        fdog += ' --force'
    fdog += ' > /dev/null 2>&1'
    try:
        subprocess.run([fdog], shell=True, check=True)
    except Exception as e:
        print(f'\033[91mProblem occurred while running fDOG for \'{jobName}\' core group:\033[0m\n{e}\nCommand used:\n{fdog}')


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
    annoToolsFcat = fcatFn.getAnnoToolFile()
    missing = []
    fdogOutDir = '%s/fcatOutput/%s/%s/fdogOutput' % (outDir, coreSet, queryID)
    out = os.listdir(fdogOutDir)
    for refSpec in out:
        if os.path.isdir(fdogOutDir + '/' + refSpec):
            # merge single extended.fa files for each refspec
            refDir = fdogOutDir + '/' + refSpec
            items = os.listdir(refDir)
            mergedFa = '%s/%s.extended.fa' % (refDir, refSpec)
            if not os.path.exists(mergedFa) or force:
                mergedFaFile = open(mergedFa, 'wb')
                for item in items:
                    if item.endswith('extended.fa'):
                        groupID = item.replace('.extended.fa','')
                        singleFa = '%s/%s.extended.fa' % (refDir, groupID)
                        if os.path.exists(singleFa):
                            with open(singleFa) as f:
                                if not queryID in f.read():
                                    missing.append(groupID)
                            shutil.copyfileobj(open(singleFa, 'rb'), mergedFaFile)
                        else:
                            missing.append(groupID)
                mergedFaFile.close()
                # calculate fas scores for merged extended.fa using fas.runFdogFas
                print('Calculating pairwise FAS scores between query orthologs and sequences of refspec...')
                fdogFAS = 'fas.runFdogFas -i %s -w %s --cores %s --redo_anno --featuretypes %s' % (mergedFa, annoDir, cpus, annoToolsFcat)
                try:
                    subprocess.run([fdogFAS], shell=True, check=True)
                except:
                    print('\033[91mProblem occurred while running fas.runFdogFas for \'%s\'\033[0m\n%s' % (mergedFa, fdogFAS))
    return(missing)


def calcFASall(coreDir, outDir, coreSet, queryID, annoDir, cpus, force):
    annoToolsFcat = fcatFn.getAnnoToolFile()
    fdogOutDir = '%s/fcatOutput/%s/%s/fdogOutput' % (outDir, coreSet, queryID)
    mergedFa = '%s/%s_all.extended.fa' % (fdogOutDir, queryID)
    count = {}
    if not os.path.exists(mergedFa) or force:
        mergedFaFile = open(mergedFa, 'w')
        out = os.listdir(fdogOutDir)
        for refSpec in out:
            if os.path.isdir(fdogOutDir + '/' + refSpec):
                refDir = fdogOutDir + '/' + refSpec
                items = os.listdir(refDir)
                for item in items:
                    groupID = item.replace('.extended.fa','')
                    # merge each ortholog seq in single extended.fa file with core group fasta file
                    # and write into mergedFaFile
                    groupFa = '%s/core_orthologs/%s/%s/%s.fa' % (coreDir, coreSet, groupID, groupID)
                    singleFa = '%s/%s.extended.fa' % (refDir, groupID)
                    if os.path.exists(singleFa) and not groupID == refSpec:
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
        mergedFaFile.close()
        if os.stat(mergedFa).st_size == 0:
            sys.exit('WARNING: No ortholog for any core gene found!')
        # calculate fas scores for merged _all.extended.fa using fas.runFdogFas
        print('Calculating FAS scores between query orthologs and all sequences in each core group...')
        fdogFAS = 'fas.runFdogFas -i %s -w %s --cores %s --redo_anno --featuretypes %s' % (mergedFa, annoDir, cpus, annoToolsFcat)
        # print(fdogFAS)
        try:
            subprocess.run([fdogFAS], shell=True, check=True)
        except:
            print('\033[91mProblem occurred while running fas.runFdogFas for \'%s\'\033[0m\n%s' % (mergedFa, fdogFAS))

# def calcFAScmd(args):
#     (seed, seedIDs, query, anno, out, name) = args
#     if not os.path.exists('%s/%s.tsv' % (out, name)):
#         cmd = 'calcFAS -s %s --seed_id %s -q %s -a %s -o %s --cpus 1 -n %s --domain > /dev/null 2>&1' % (seed, seedIDs, query, anno, out, name)
#         try:
#             subprocess.run([cmd], shell=True, check=True)
#         except:
#             print('\033[91mProblem occurred while running calcFAS\033[0m\n%s' % (cmd))
#
# def calcFAScons(coreDir, outDir, coreSet, queryID, annoDir, cpus, force):
#     # output files
#     (mode, phyloprofileDir) = outputMode(outDir, coreSet, queryID, force, 'other')
#     if mode == 1 or mode == 3:
#         finalPhyloprofile = open('%s/mode4.phyloprofile' % (phyloprofileDir), 'w')
#         finalPhyloprofile.write('geneID\tncbiID\torthoID\tFAS\n')
#     elif mode == 2:
#         finalPhyloprofile = open('%s/mode4.phyloprofile' % (phyloprofileDir), 'a')
#     # parse single fdog output
#     missing = []
#     fdogOutDir = '%s/fcatOutput/%s/%s/fdogOutput' % (outDir, coreSet, queryID)
#     calcFASjob = []
#
#     annoDirTmp = '%s/fcatOutput/%s/%s/tmp/anno' % (outDir, coreSet, queryID)
#     Path(annoDirTmp).mkdir(parents=True, exist_ok=True)
#     fasDirOutTmp = '%s/fcatOutput/%s/%s/tmp/fasOut' % (outDir, coreSet, queryID)
#     Path(fasDirOutTmp).mkdir(parents=True, exist_ok=True)
#
#     fdogOutDir = '%s/fcatOutput/%s/%s/fdogOutput' % (outDir, coreSet, queryID)
#     out = os.listdir(fdogOutDir)
#     missing = []
#     for refSpec in out:
#         if os.path.isdir(fdogOutDir + '/' + refSpec):
#             # make calcFAS job for founded orthologs and consensus seq of each group
#             refDir = fdogOutDir + '/' + refSpec
#             groups = os.listdir(refDir)
#             groupFa = '%s/%s.fa' % (annoDirTmp, refSpec)
#             groupFaFile = open(groupFa, 'w')
#             for groupID in groups:
#                 if os.path.isdir(refDir + '/' + groupID):
#                     # get seed and query fasta
#                     singleFa = '%s/%s/%s.extended.fa' % (refDir, groupID, groupID)
#                     if os.path.exists(singleFa):
#                         consFa = '%s/core_orthologs/%s/%s/fas_dir/annotation_dir/cons.fa' % (coreDir, coreSet, groupID)
#                         consFaLink = '%s/cons_%s.fa' % (annoDirTmp, groupID)
#                         fcatFn.checkFileExist(consFa, '')
#                         try:
#                             os.symlink(consFa, consFaLink)
#                         except FileExistsError:
#                             os.remove(consFaLink)
#                             os.symlink(consFa, consFaLink)
#                         seedID = []
#                         for s in SeqIO.parse(singleFa, 'fasta'):
#                             if queryID in s.id:
#                                 idTmp = s.id.split('|')
#                                 seedID.append(idTmp[-2])
#                                 groupFaFile.write('>%s\n%s\n' % (idTmp[-2], s.seq))
#                         # get annotations for seed and query
#                         consJson = '%s/core_orthologs/%s/%s/fas_dir/annotation_dir/cons.json' % (coreDir, coreSet, groupID)
#                         consJsonLink = '%s/cons_%s.json' % (annoDirTmp, groupID)
#                         fcatFn.checkFileExist(consJson, '')
#                         try:
#                             os.symlink(consJson, consJsonLink)
#                         except FileExistsError:
#                             os.remove(consJsonLink)
#                             os.symlink(consJson, consJsonLink)
#                         # tmp fas output
#                         calcFASjob.append([groupFa, ' '.join(seedID), consFaLink, annoDirTmp, fasDirOutTmp, groupID])
#                     else:
#                         missing.append(groupID)
#             groupFaFile.close()
#             # get annotation for orthologs
#             if not os.path.exists('%s/%s.json' % (annoDirTmp, refSpec)):
#                 extractAnnoCmd = 'fas.doAnno -i %s -o %s -e -a %s/%s.json -n %s > /dev/null 2>&1' % (groupFa, annoDirTmp, annoDir, queryID, refSpec)
#                 try:
#                     subprocess.run([extractAnnoCmd], shell=True, check=True)
#                 except:
#                     print('\033[91mProblem occurred while running extracting annotation for \'%s\'\033[0m\n%s' % (groupFa, extractAnnoCmd))
#     # do FAS calculation
#     pool = mp.Pool(cpus)
#     calcFASout = []
#     for _ in tqdm(pool.imap_unordered(calcFAScmd, calcFASjob), total=len(calcFASjob)):
#         calcFASout.append(_)
#     # parse fas output into phyloprofile
#     for tsv in os.listdir(fasDirOutTmp):
#         if os.path.isfile('%s/%s' % (fasDirOutTmp, tsv)):
#             for line in fcatFn.readFile('%s/%s' % (fasDirOutTmp, tsv)):
#                 if not line.split('\t')[0] == 'Seed':
#                     groupID = tsv.split('.')[0]
#                     ncbiID = 'ncbi' + str(queryID.split('@')[1])
#                     orthoID = line.split('\t')[0]
#                     fas = fcatFn.roundTo4(float(line.split('\t')[2].split('/')[0]))
#                     finalPhyloprofile.write('%s\t%s\t%s\t%s\n' % (groupID, ncbiID, orthoID, fas))
#     finalPhyloprofile.close()

def checkResult(fcatOut, force):
    if force:
        fcatFn.deleteFolder('%s/fdogOutput.tar.gz' % fcatOut)
        fcatFn.deleteFolder('%s/fdogOutput' % fcatOut)
        fcatFn.deleteFolder('%s/phyloprofileOutput' % fcatOut)
        return(0)
    else:
        if not os.path.exists('%s/phyloprofileOutput/mode1.phyloprofile' % fcatOut):
            if os.path.exists('%s/fdogOutput.tar.gz' % fcatOut):
                return(1)
            else:
                # if os.path.exists(fcatOut):
                #     shutil.rmtree(fcatOut)
                return(0)
        else:
            return(2)

def searchOrtho(args):
    coreDir = os.path.abspath(args.coreDir)
    coreSet = args.coreSet
    fcatFn.checkFileExist(coreDir + '/core_orthologs/' + coreSet, '')
    refspecList = str(args.refspecList).split(",")
    if len(refspecList) == 0:
        sys.exit('No refefence species given! Please specify reference taxa using --refspecList option!')
    query = args.querySpecies
    fcatFn.checkFileExist(os.path.abspath(query), '')
    query = os.path.abspath(query)
    taxid = str(args.taxid)
    outDir = args.outDir
    if outDir == '':
        outDir = os.getcwd()
    else:
        Path(outDir).mkdir(parents=True, exist_ok=True)
    outDir = os.path.abspath(outDir)
    coreTaxa_dir = args.coretaxadb
    if coreTaxa_dir == '':
        coreTaxa_dir = '%s/coreTaxa_dir' % coreDir
    coreTaxa_dir = os.path.abspath(coreTaxa_dir)
    fcatFn.checkFileExist(coreTaxa_dir, 'Please set path to blastDB using --coreTaxa_dir option.')
    annoDir = args.annoDir
    if annoDir == '' or annoDir == '%s/annotation_dir' % os.path.abspath(args.coreDir) or annoDir == '%s/annotation_dir/' % os.path.abspath(args.coreDir):
        annoDir = '%s/fcatOutput/%s/annotation_dir' % (outDir, args.coreSet)
        Path(annoDir).mkdir(parents=True, exist_ok=True)
        # annoDir = '%s/annotation_dir' % coreDir
    annoDir = os.path.abspath(annoDir)
    fcatFn.checkFileExist(annoDir, 'Please set path to annotation directory using --annoDir option.')
    for annoFile in glob.glob('%s/annotation_dir/*.json' % args.coreDir):
        annoFileName = annoFile.split('/')[-1]
        if not os.path.exists('%s/%s' % (annoDir, annoFileName)):
            try:
                os.symlink('%s/annotation_dir/%s' % (args.coreDir, annoFileName), '%s/%s' % (annoDir, annoFileName))
            except FileExistsError:
                os.remove('%s/%s' % (annoDir, annoFileName))
                os.symlink('%s/annotation_dir/%s' % (args.coreDir, annoFileName), '%s/%s' % (annoDir, annoFileName))
    annoQuery = args.annoQuery

    cpus = args.cpus
    if cpus >= mp.cpu_count():
        cpus = mp.cpu_count()-1
    force = args.force
    keep = args.keep

    currDir = os.getcwd()

    # check annotation of query species and get query ID
    (doAnno, queryTaxId) = checkQueryAnno(annoQuery, annoDir, taxid, query)
    queryID = parseQueryFa(coreSet, query, annoQuery, taxid, outDir, doAnno, annoDir, cpus)
    if doAnno == False:
        if os.path.exists( '%s/query_%s.json' % (annoDir, queryTaxId)):
            os.rename('%s/query_%s.json' % (annoDir, queryTaxId), annoDir+'/'+queryID+'.json')
    # empty old directory if force is specified
    if force:
        fcatFn.deleteFolder(f'{outDir}/fcatOutput/{coreSet}/{queryID}')
        # if os.path.exists(f'{outDir}/fcatOutput/{coreSet}/{queryID}'):
        #     shutil.rmtree(f'{outDir}/fcatOutput/{coreSet}/{queryID}')
    # move searchTaxa_dir into fcatOutput/coreSet/query folder
    src = '%s/searchTaxa_dir/%s' % (outDir, queryID)
    dest = '%s/fcatOutput/%s/%s/searchTaxa_dir/%s' % (outDir, coreSet, queryID, queryID)
    if not os.path.exists(dest):
        shutil.move(src, dest)
    # check old output files
    fcatOut = '%s/fcatOutput/%s/%s' % (outDir, coreSet, queryID)
    status = checkResult(fcatOut, force)
    print('Preparing...')
    groupRefspec = {}
    if status == 0:
        (fdogJobs, ignored, groupRefspec) = prepareJob(coreDir, coreSet, queryID, refspecList, outDir, coreTaxa_dir, annoDir, annoQuery, force, cpus)
        print('Searching orthologs...')
        os.chdir('%s/fcatOutput/%s/%s' % (outDir, coreSet, queryID))
        fdogOut = []
        if cpus == 1:
            for job in fdogJobs:
                tmp = runFdog(job)
                fdogOut.append(tmp)
        else:
            pool = mp.Pool(cpus)
            for _ in tqdm(pool.imap_unordered(runFdog, fdogJobs), total=len(fdogJobs)):
                fdogOut.append(_)
            pool.close()
            pool.join()
        # write ignored groups and refspec for each group based on given refspec list
        ignoredFile = open('%s/fcatOutput/%s/%s/ignored.txt' % (outDir, coreSet, queryID), 'w')
        if len(ignored) > 0:
            # print('\033[92mNo species in %s found in core set(s): %s\033[0m' % (refspecList, ','.join(ignored)))
            ignoredFile.write('\n'.join(ignored))
            ignoredFile.write('\n')
        ignoredFile.close()
        if len(groupRefspec) > 0:
            refspecFile = open('%s/fcatOutput/%s/%s/last_refspec.txt' % (outDir, coreSet, queryID), 'w')
            for g in groupRefspec:
                refspecFile.write('%s\t%s\n' % (g, groupRefspec[g]))
            refspecFile.close()

    # if not status == 2:
        # Calculating FAS scores between query orthologs and all sequences in each core group
        calcFASall(coreDir, outDir, coreSet, queryID, annoDir, cpus, force)
        # Calculating pairwise FAS scores between query orthologs and sequences of refspec
        missing = calcFAS(coreDir, outDir, coreSet, queryID, annoDir, cpus, force)
        # print('Calculating FAS scores between query orthologs and consensus sequence in each core group...')
        # calcFAScons(coreDir, outDir, coreSet, queryID, annoDir, cpus, force)
        # remove tmp folder
        if os.path.exists('%s/tmp' % fcatOut):
            shutil.rmtree('%s/tmp' % fcatOut)
        # write missing groups
        missingFile = open('%s/fcatOutput/%s/%s/missing.txt' % (outDir, coreSet, queryID), 'w')
        if len(missing) > 0:
            missingFile.write('\n'.join(missing))
            missingFile.write('\n')
        missingFile.close()

    if os.path.exists('%s/fdogOutput' % fcatOut):
        try:
            make_archive('%s/fdogOutput' % fcatOut, '%s/fdogOutput.tar.gz' % fcatOut, 'gztar')
        except:
            print('Cannot archiving fdog output!')

    if keep == False:
        fcatFn.deleteFolder('%s/fcatOutput/%s/%s/searchTaxa_dir' % (outDir, coreSet, queryID))
        fcatFn.deleteFolder('%s/searchTaxa_dir' % outDir)
        # # print('Cleaning up...') ### no idea why rmtree not works :(
        # if os.path.exists('%s/fcatOutput/%s/%s/searchTaxa_dir' % (outDir, coreSet, queryID)):
        #     shutil.rmtree('%s/fcatOutput/%s/%s/searchTaxa_dir' % (outDir, coreSet, queryID))
    os.chdir(currDir)
    print('Done! Check output in %s' % fcatOut)

def main():
    version = get_distribution('fcat').version
    parser = argparse.ArgumentParser(description='You are running fcat version ' + str(version) + '.')
    required = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('optional arguments')
    required.add_argument('-d', '--coreDir', help='Path to core set directory, where folder core_orthologs can be found', action='store', default='', required=True)
    required.add_argument('-c', '--coreSet', help='Name of core set, which is subfolder within coreDir/core_orthologs/ directory', action='store', default='', required=True)
    required.add_argument('-r', '--refspecList', help='List of reference species', action='store', default='', required=True)
    required.add_argument('-q', '--querySpecies', help='Path to gene set for species of interest', action='store', default='', required=True)
    optional.add_argument('-o', '--outDir', help='Path to output directory', action='store', default='')
    optional.add_argument('--coretaxadb', help='Path to BLAST directory of all core species', action='store', default='')
    optional.add_argument('--annoDir', help='Path to FAS annotation directory', action='store', default='')
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
