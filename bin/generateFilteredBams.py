#!/usr/bin/env python

import os
import sys
import argparse as ag

# generate clean bams for transcript assembly

from lncmodule import selectTranscriptsByType
from lncmodule import bedtoolsClean

def generateFilteredBams(gtfFile, bamList, removeTypes, nThread, outputDir):
    # create a temporal gtf file holding the regions to be removed
    print >> sys.stdout, 'Generate temporal GTF file used in bedtools intersect.'
    tGtfPath = selectTranscriptsByType(removeTypes, gtfFile, outputDir)

    # start bedtools intersect
    print >> sys.stdout, 'Start bedtools intersect.'
    success = bedtoolsClean(bamList, tGtfPath, nThread, outputDir)

    if success == 0:
        print >> sys.stdout, 'Filtered BAM files are created successfully!'


def main():
    parser = ag.ArgumentParser(prog='generateFilteredBams', description='Generate filtered bam files that does not contain reads in selected regions.')
    parser.add_argument('-g', dest='inputGtf', required=True, default='', type=str, help='The path to reference gene annotation file in GTF format. It can be downloaded from gencode website.')
    parser.add_argument('-t', dest='types',required=False, default='protein_coding', nargs='+', type=str, help='The gene types to be removed from BAM files.')
    parser.add_argument('-o', dest= 'outputDir',required=False, default='FLORA', type=str, help='Output directory for the output files')
    parser.add_argument('-n', dest='nThread', required=False, default=1, type=int, help='Number of threads to be used for bedtools intersect.')
    parser.add_argument("inputBams", nargs=1, type=str, help='The path to a txt file that contains paths to input bam files separated by the newline character.')
    
    # process parameters
    args = parser.parse_args()
    bams = args.inputBams[0]

    if not os.path.isfile(bams):
        sys.exit('Path to BAM list file is not valid!')
    if not os.path.isfile(args.inputGtf):
        sys.exit('Path to GTF file is not valid!')
    
    nThread = args.nThread
    if nThread <= 0:
        nThread = 1

    generateFilteredBams(args.inputGtf, bams, args.types, nThread, args.outputDir)
    print >> sys.stdout, "FLORA is finished!"




    

if __name__ == '__main__':
    main()
