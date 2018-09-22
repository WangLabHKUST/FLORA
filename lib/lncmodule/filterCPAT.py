# run CPAT
# summarize CPAT result

import sys
import subprocess

# cpatList contains 4 elements
# cpatList: fasta referece path, hexamer path, logit model path
def filterCPAT(inputFileName, cpatList, cutoff):
    outputList = list()
    outputName = inputFileName + '.cpat'
    try:
        cpatReturn = subprocess.call(['cpat.py', '-g', inputFileName, '-o', outputName, '-x', cpatList[1], '-d', cpatList[2], '-r', cpatList[0]])
    except:
        sys.exit('RUN CPAT failed! Please check the input parameters')
    if cpatReturn == 0:
        with open(outputName, 'r') as cpatResult:
            for line in cpatResult:
                if line[0:9] == 'mRNA_size':
                    continue
                else:
                    lineElements = line.strip().split('\t')
                    prob = float(lineElements[5])
                    if prob <= cutoff:
                        outputList.append(lineElements[0])
    return outputList
