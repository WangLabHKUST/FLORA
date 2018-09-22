#!/usr/bin/env python

# annotate transcripts with reference annotation from GENCODE or RefSeq
# annotation files can be: 1. GFF file from RefSeq 2. GFF file from GENCODE


import os
import sys
import argparse as ag

from lncmodule import organizeGencodeGFF
from lncmodule import organizeRefSeqGFF
from lncmodule import organizeInput
from lncmodule import annotateGenes

def annotateTranscripts(referenceType, referencePath, refseqIdentifier, inputGTF, output, distance, number):
    # organize reference file
    if referenceType == 'refseq':
        organizedReference = organizeRefSeqGFF(referencePath, refseqIdentifier)
    if referenceType == 'gencode':
        organizedReference = organizeGencodeGFF(referencePath)
    
    # organize inputGTF GTF file
    organizedInput = organizeInput(inputGTF)

    # generate annotation file
    result = annotateGenes(organizedReference, organizedInput, output, distance, number)

    return result

def main():
    parser = ag.ArgumentParser(prog='annotateTranscripts', description='Annotated assembled transcripts with reference annotation from RefSeq or GENCODE')
    parser.add_argument('-r', dest='reference', required=True, default='refseq', type=str, help='Identify the source of reference annotation. Currently, it should be either refseq or gencode')
    parser.add_argument('-f', dest='path', required=True, default='', type=str, help='The path to the reference annotation file.')
    parser.add_argument('-i', dest='identifier', required=False, default='', type=str, help='The path to the RefSeq assembly report. The file is used to replace RefSeq sequence identifiers with USCS identifiers. Required if you use RefSeq annotation.')
    parser.add_argument('-d', dest='distance', required=False, default=10000, type=int, help='Specify the distance within which the nearby genes should locate. Default: 10,000 bp. Unit: bp.')
    parser.add_argument('-n', dest='number', required=False, default=1, type=int, help='Specify the number of nearby genes obtained. Default: 1.')
    parser.add_argument('-o', dest= 'output',required=False, default='FLORA_Annotation.txt', type=str, help='Output file name for the final annotation. Default: FLORA_annotation.txt')
    parser.add_argument("inputGTF", nargs=1, default='',type=str, help='Input transcriptome file in GTF format. The file could be output from Cufflinks or StringTie')
    
    # process parameters
    args = parser.parse_args()

    inputGTF = args.inputGTF[0]
    output = args.output
    reference = args.reference
    path = args.path
    identifier = args.identifier

    distance = args.distance
    number = args.number

    if not os.path.isfile(inputGTF):
        sys.exit('Path to input GTF file is not valid!')
    if not os.path.isfile(path):
        sys.exit('Path to reference file is not valid!')
    if reference != 'refseq' and reference != 'gencode':
        sys.exit('reference has to be either refseq or gencode!')
    if reference == 'refseq' and identifier == '':
        sys.exit('Please provide RefSeq assembly report if you use RefSeq annotation')
    if identifier != '' and not os.path.isfile(identifier):
        sys.exit('Path to RefSeq assembly report is not valid.')

    # start annotation
    result = annotateTranscripts(reference, path, identifier, inputGTF, output, distance, number)
    if result == 0:
        print >> sys.stdout, "FLORA annotationTranscritps is finished"
    return result


if __name__ == '__main__':
    main()
