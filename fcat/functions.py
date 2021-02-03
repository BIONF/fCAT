# functions used for fCAT
import sys
import os
import shutil

def checkFileExist(file, msg):
    if not os.path.exists(os.path.abspath(file)):
        sys.exit('%s not found! %s' % (file, msg))

def readFile(file):
    if os.path.exists(file):
        with open(file, 'r') as f:
            lines = f.read().splitlines()
            f.close()
            return(lines)
    else:
        sys.exit('%s not found' % file)

def removeEmptyLine(file):
    if os.path.exists(file):
        nf = open(file+'.temp','w')
        with open(file,'r+') as f:
            for line in f:
                if not len(line.strip()) == 0:
                    nf.write(line)
        nf.close()
        os.replace(file+'.temp', file)


def removeDupEmpty(file):
    if os.path.exists(file):
        lines_seen = set() # holds lines already seen
        outfile = open(file+'.temp', 'w')
        for line in open(file, 'r'):
            if line not in lines_seen: # not a duplicate
                if not len(line.strip()) == 0: # not empty
                    outfile.write(line)
                    lines_seen.add(line)
        outfile.close()
        os.replace(file+'.temp', file)


def roundTo4(number):
    return("%.4f" % round(number, 4))

def deleteFolder(folder):
    if os.path.exists(folder):
        if os.path.isfile(folder):
            os.remove(folder)
        else:
            shutil.rmtree(folder)
