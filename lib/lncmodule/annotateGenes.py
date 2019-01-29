# with organized reference GFF and organized inputGTF GTF, do annotation with reference for inputGTF
# annotation: find the overlapped genomic elements or nearby genomic elements for each gene in inputGTF GTF file

# the default value for distance is 1 Mbp, the default value for number of nearby elements reported is 1
def annotateGenes(reference, inputGTF, output, distanceLimit, numberLimit):
    # loop through the inputGTF to do annotation
    writer = open(output, 'w+')
    writer.write('#gene_id\tnearby_gene_name\tnearby_gene_distance\tnearby_gene_type\n')

    for key in inputGTF:
        try:
            currentReference = reference[key]
        except KeyError:
            print "Genomic Region Identifier " + key + " is not present in reference GFF file"
            continue
        for i, elem in enumerate(inputGTF[key]):
            # for each gene, check with currentReference
            # if currentReference is empty
            if not currentReference:
                print "Genomic Region identified by " + key + " does not contain any annotation"
                continue
            result = findNearbyGenomicElements(elem, currentReference, distanceLimit, numberLimit)
            # print results to output
            printToOutput(elem, result, writer)
    writer.close()
    return 0


def findNearbyGenomicElements(gene, reference, distance, number):
    # handle the cornered case: gene is the first & gene is the last
    firstT = len(reference) - 1
    secondT = 0
    # assume the reference is sorted
    # find the first transition point
    for i, elem in enumerate(reference):
        if elem[1] < gene[0]:
            # not yet the first transition
            continue
        else:
            firstT = i
            break
    # from the first transition point find the second transition point
    for i, elem in enumerate(reference):
        if elem[0] < gene[1]:
            continue
        else:
            secondT = i
            break
    # from the first transition point, check the genomic elements ahead
    result = []

    lowerBound = max(firstT-1-number, -1)
    for i in range(firstT-1, lowerBound, -1):
        if gene[0] - reference[i][1] <= distance:
            add = (reference[i][2], reference[i][3], reference[i][4], reference[i][1] - gene[0])
            result.insert(0, add)
        else:
            break
    # check the genomic elements between the two transition points
    for i in range(firstT, secondT):
        add = (reference[i][2], reference[i][3], reference[i][4], 0)
        result.append(add)
    
    # from the second transition point, check the following genomic elements 
    upperBound = min(secondT+number, len(reference))
    for i in range(secondT, upperBound):
        if reference[i][0] - gene[1] <= distance:
            add = (reference[i][2], reference[i][3], reference[i][4], reference[i][0] - gene[1])
            result.append(add)
    return result


def printToOutput(inputGene, outputData, outputFile):
    # summirize the nearby gene name, gene distance, gene type
    geneNameList = [i[1] for i in outputData]
    geneDistanceList = [str(i[3]) for i in outputData]
    geneTypeList = [i[2] for i in outputData]

    # generate output string
    output = inputGene[2] + '\t' + ';'.join(geneNameList) + '\t' + ';'.join(geneDistanceList) + '\t' + ';'.join(geneTypeList) + '\n'
    outputFile.write(output)


# test
# 
# from organizeGencodeGFF import organizeGencodeGFF
# from organizeRefSeqGFF import organizeRefSeqGFF
# from organizeAnnotationInput import organizeInput
# 
# refseq = organizeGencodeGFF('test/gencode.v27.long_noncoding_RNAs.gff3')
# inputFile = organizeInput('test/TCGA_GC_lncRNA_mainChr_novel.gtf')
# 
# annotateGenes(refseq, inputFile, 'outputtest.txt', 10000, 4)
# 
