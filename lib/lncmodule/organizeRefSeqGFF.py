# organize the RefSeq reference GFF file into a library containing multiple lists
# each list contains gene information from one chromosome
# information of each gene is stored in a tuple: (start, end, gene_id, gene_name, gene_type)


from handleGFF import getElementGFF

# summarize genomic elements with no parent
def organizeRefSeqGFF(referencePath, identifierPath):
    chromosomeMap = mapIdentifier(identifierPath)
    result = {}
    with open(referencePath, 'r') as reference:
        regionIdentifier = ('', '')
        for line in reference:
            # Pass if the line is comment
            if line[0] == '#':
                continue
            # Pass if the line is a child element
            elif line.find('Parent=') >= 0:
                continue
            else:
                elements = line.strip().split('\t')
                 # if the line contians region information
                if elements[2] == 'region':
                    # update the chromosome information
                    if elements[0] in chromosomeMap:
                        chromosome = chromosomeMap[elements[0]]
                        regionIdentifier = (chromosome, elements[0])
                        if regionIdentifier[0] not in result:
                            result[regionIdentifier[0]] = []
                    else:
                        continue
                # if the line is a gene element
                else:
                    # check whether the chromosome identifier match, if not match, continue this gene
                    if elements[0] != regionIdentifier[1]:
                        continue
                    else:
                        # get (start, end, gene_id, gene_name, gene_type)
                        start = int(elements[3])
                        end = int(elements[4])
                        gene_id = getElementGFF(line, 'ID')
                        gene_name = getElementGFF(line, 'Name')
                        gene_type = getElementGFF(line, 'gene_biotype')
                        if gene_id=='' or gene_name == '' or gene_type == '':
                            continue
                        toBeAdded = (start, end, gene_id, gene_name, gene_type)
                        result[regionIdentifier[0]].append(toBeAdded)

    for key in result:
        result[key].sort()
    return result

# map the RefSeq region identifier with UCSC region identifier
# eg. NC_000001.11 <=> chr1
def mapIdentifier(identifierFile):
    result = {}
    with open(identifierFile, 'r') as f:
        for line in f:
            if line[0] == '#':
                continue
            else:
                elements = line.strip().split('\t')
                if elements[9] in ('na', 'NA', 'Na'):
                    continue
                else:
                    result[elements[6]] = elements[9]
    return result


