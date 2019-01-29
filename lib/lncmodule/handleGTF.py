def getGtfComponent(input, keyword):
    geneIdPosition = input.find(keyword)
    if geneIdPosition < 0:
        return ''
    length = len(keyword) + 2
    geneIdEndPosition = input[geneIdPosition+length:].find('"')
    return input[geneIdPosition+length:geneIdPosition+length+geneIdEndPosition]
