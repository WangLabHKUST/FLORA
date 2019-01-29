# select genes by gene type
# after selection, use the selected genes to clean bam file

import sys
import os.path
from subprocess import call
from multiprocessing import Pool
from functools import partial
from handleGTF import getGtfComponent


# this function will generate a temporal file to store the selected genes
def selectTranscriptsByType(typeList, gtfFile, outputDirectory):
    if not os.path.isdir(outputDirectory):
        os.mkdir(outputDirectory)
    
    if outputDirectory[-1] != '/':
        outputDirectory = outputDirectory + '/'

    # generate the correct name for gtf.selected output
    gtfFileName = gtfFile[gtfFile.rfind('/')+1:]
    outputFile = open(outputDirectory + gtfFileName + '.selected', 'w')
    outputFile.close()
    outputFile = open(outputDirectory + gtfFileName+'.selected', 'a')

    holdGeneId = ''
    with open(gtfFile, 'r') as inputFile:
        for line in inputFile:
            if line[0] == '#':
                continue
            else:
                lineElements = line.strip().split('\t')
                if lineElements[2] == 'gene':
                    # check type
                    geneType = getGtfComponent(line, 'gene_type')
                    if geneType in typeList:
                        holdGeneId = getGtfComponent(line, 'gene_id')
                        outputFile.write(line)
                else:
                    lineId = getGtfComponent(line, 'gene_id')
                    if lineId == holdGeneId:
                        outputFile.write(line)
    outputFile.close()
    return outputDirectory + gtfFileName+'.selected'
    
    # start bedtools intersect with gtfFile.selected

def bedtoolsClean(bamlist, gtfFile, nThread, outputDirectory):
    if not os.path.isdir(outputDirectory):
        os.mkdir(outputDirectory)
    with open(bamlist, 'r') as p:
        paths = p.read().splitlines()
    paths = [x for x in paths if x != '']

    if outputDirectory[-1] != '/':
        outputDirectory = outputDirectory + '/'

    # save commands to a list
    commands = list()
    for i in paths:
        outputFilePath = outputDirectory + i[i.rfind('/')+1:] + '.clean.bam'
        command = 'bedtools intersect -a ' + i + ' -b ' + gtfFile + ' -v' + ' > ' + outputFilePath
        commands.append(command)
    
    pool = Pool(nThread)
    for i, returncode in enumerate(pool.imap_unordered(partial(call, shell=True), commands)):
        if returncode != 0:
            print("%d command failed: %d" % (i, returncode))

    return 0
            


