#!/usr/bin/env python

# generate the keep transcript list
# create the gtf with the keep transcripts

import os
import sys
import argparse as ag

# import self-defined module
from lncmodule import gtfToCpatBed
from lncmodule import filterExonNumberAndLength
from lncmodule import filterCPAT
from lncmodule import getTranscriptId

# call functions in lncmodule to finish the workflow
def generateLncTranscriptome(inputFileName, exonNumberCutoff, lengthCutoff, cpatCutoff, cpatParameters, outputfileName):
    if not os.path.isdir('tmp/'):
        os.mkdir('tmp/')
    # convert the input gtf file (inputFileName) into bed file
    print >> sys.stdout, 'convert input gtf file into bed file'
    newName = inputFileName.split("/")[-1]
    bedFilePositive = 'tmp/'+ inputFileName + '.positive.bed'
    bedFileNegative = 'tmp/'+ inputFileName + '.negative.bed'
    gtfToCpatBed(newName, bedFilePositive, '+')
    gtfToCpatBed(newName, bedFileNegative, '-')
    # summarize the exon number and transcript length from bed file
    # generate the list that contains transcript id that fit exonNumberCutoff and lengthCutoff
    print >> sys.stdout, 'summarize the exon number and length of transcripts'
    exonSummarize = 'tmp/' + newName + '.exonAndLength.txt'
    exonFiltered = filterExonNumberAndLength(bedFilePositive, exonNumberCutoff, lengthCutoff, exonSummarize)

    # run CPAT
    print >> sys.stdout, 'start run CPAT'
    cpatPositiveFiltered = filterCPAT(bedFilePositive, cpatParameters, cpatCutoff)
    cpatNegativeFiltered = filterCPAT(bedFileNegative, cpatParameters, cpatCutoff)


    transcriptFiltered = set.intersection(set(exonFiltered), set(cpatPositiveFiltered), set(cpatNegativeFiltered))
    
    # keep transcripts in transcriptFiltered and write them into output
    print >> sys.stdout, 'output presumed long noncoding RNAs'
    output = open(outputfileName, 'w')
    output.close()
    output =  open(outputfileName, 'a')
    for line in open(inputFileName, 'r'):
        if line[0] == '#':
            continue
        transcriptId = getTranscriptId(line)
        if transcriptId in transcriptFiltered:
            output.write(line)
    output.close()

# python filterTranscripts.py -r reference.fa -e exon_number_cutoff -l transcript_length_cutoff -c cpat_cutoff -h hexamer_table -m logit_model -o output_name inputGTF
def main():
    parser = ag.ArgumentParser(prog='filterTranscript', description='Filter out non-lncRNA transcripts after assembly')
    parser.add_argument('-r', dest='reference', required=True, default='', type=str, help='The path to reference genome sequences in FASTA format. It will be indexed automatically by CPAT if .fai file is not present.')
    parser.add_argument('-e', dest='exon',required=False, default=2, type=int, help='The least number of exons a transcipt should have in order to be kept in the final transcriptome')
    parser.add_argument('-l', dest='length',required=False, default=200, type=int, help='The shortest transcript to be kept in the final transcriptome')
    parser.add_argument('-c', dest='cpat',required=False, default=0.364, type=float, help='CPAT cutoff used to filter transcripts')
    parser.add_argument('-x', dest='hexamer',required=True, default='', type=str, help='The path to hexamer table required by CPAT. Can be downloaded from CPAT website.')
    parser.add_argument('-m', dest='logit', required=True, default='', type=str, help='The path to logit model required by CPAT. Can be downloaded CPAT website.')
    parser.add_argument('-o', dest= 'outputGTF',required=False, default='FLORA', type=str, help='Output prefix for the final transcriptome GTF file')
    parser.add_argument("inputGTF", nargs=1, default='',type=str, help='Input transcriptome file in GTF format. Can be output from stringtie or cufflinks')
    
    # process parameters
    args = parser.parse_args()

    inputGTF = args.inputGTF[0]
    outputGTF = args.outputGTF + '.gtf'
    hexamer_table = args.hexamer
    logit_model = args.logit
    reference = args.reference

    if not os.path.isfile(inputGTF):
        sys.exit('Path to input GTF file is not valid!')
    if not os.path.isfile(hexamer_table):
        sys.exit('Path to hexamer table is not valid!')
    if not os.path.isfile(logit_model):
        sys.exit('Path to logit model is not valie!')
    if not os.path.isfile(reference):
        sys.exit('Path to reference is not valie!')

    cpat_parameters = [reference, hexamer_table, logit_model]

    generateLncTranscriptome(inputGTF, args.exon, args.length, args.cpat, cpat_parameters, outputGTF)
    print >> sys.stdout, "FLORA filterTranscripts is finished!"
    return 0

    

if __name__ == '__main__':
    main()

#input_file = 'TCGA_GC_v27_merged.gtf'
#exon_number_cutoff = 2
#length_cutoff = 300
#cpat_cutoff = 0.3
#cpat_parameters = ['GRCh38.d1.vd1.fa', '/Users/hongyushi/Downloads/NORI/inst/extdata/Human_Hexamer.tab', '/Users/hongyushi//Downloads/NORI/inst/extdata/Human_train.RData']
#output_file = 'TCGA_GC_v27_lncRNA_filtered.gtf'

#generateLncTranscriptome(input_file, exon_number_cutoff, length_cutoff, cpat_cutoff, cpat_parameters, output_file)

#import subprocess
#subprocess.call('date')
