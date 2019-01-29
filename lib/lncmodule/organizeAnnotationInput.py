# organize input GTF file for annotation
# The sample of GTF file is from the output of StringTie
# for each gene, store its smallest coordinate and the largest coordinate

from handleGTF import getGtfComponent

def organizeInput(path):
    result = {}
    geneHold = ''
    chromosomeHold = ''
    minCo = 0
    maxCo = 0
    with open(path, 'r') as input:
        for line in input:
            if line[0] == '#':
                continue
            else:
                gene_id = getGtfComponent(line, 'gene_id')
                if gene_id == '':
                    continue
                else:
                    elements = line.strip().split('\t')
                    if geneHold == '' or geneHold != gene_id:
                        # store the results on hold
                        if geneHold != '':
                            toBeAdded = (minCo, maxCo, geneHold)
                            if chromosomeHold not in result:
                                result[chromosomeHold] = []
                            result[chromosomeHold].append(toBeAdded)
                        geneHold = gene_id
                        chromosomeHold = elements[0]
                        minCo = int(elements[3])
                        maxCo = int(elements[4])
                    else:
                        minCo = min(minCo, int(elements[3]))
                        maxCo = max(maxCo, int(elements[4]))
    # finish the last round
    toBeAdded = (minCo, maxCo, geneHold)
    if chromosomeHold not in result:
        result[chromosomeHold] = []
    result[chromosomeHold].append(toBeAdded)
    return result

