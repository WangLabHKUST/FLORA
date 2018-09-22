# record the transcript length and exon number of the transcript
# read transcript information from bed file

# output a tuple of filtered transcript id
def filterExonNumberAndLength(inputFileName, exonNumberCutoff, lengthCutoff, outputSummarize):
    outputId = list()
    output = open(outputSummarize,'w')
    output.close()
    output = open(outputSummarize, 'a')

    with open(inputFileName, 'r') as input:
        for line in input:
            if(line[0:3] == 'chr'):
                lineElements = line.strip().split('\t')
                exon = int(lineElements[9])
                lengthList = lineElements[10].strip().split(',')
                lengthList = lengthList[0:-1]
                length = 0
                for i in lengthList:
                    length = length + int(i)
                transcriptId = lineElements[3]
                output.write(transcriptId + '\t' + str(exon) + '\t' + str(length) + '\n')
                if(exon >= exonNumberCutoff and length >= lengthCutoff):
                    outputId.append(transcriptId)
    output.close()
    return outputId
