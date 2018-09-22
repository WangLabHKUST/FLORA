# organize the GENCODE reference GFF3 file into a library containing multiple lists
# each list contains gene information from one chromosome
# information of each gene is stored in a tuple: (start, end, gene_id, gene_name, gene_type)

from handleGFF import getElementGFF

def organizeGencodeGFF(referencePath):
    result = {}
    with open(referencePath, 'r') as reference:
        for line in reference:
            # continue if the line is comment
            if line[0] == '#':
                continue
            else:
                elements = line.strip().split('\t')
                if elements[2] == 'gene':
                    start = int(elements[3])
                    end = int(elements[4])
                    gene_id = getElementGFF(line, 'ID')
                    gene_name = getElementGFF(line, 'gene_name')
                    gene_type = getElementGFF(line, 'gene_type')
                    toBeAdded = (start, end, gene_id, gene_name, gene_type)
                    if elements[0] not in result:
                        result[elements[0]] = []
                    result[elements[0]].append(toBeAdded)
    for key in result:
        result[key].sort()
    return result
