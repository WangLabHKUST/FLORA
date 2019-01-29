def getElementGFF(line, element):
    length = len(element) + 1
    position = line.find(element+'=')
    if position == -1:
        return ''
    else:
        endPosition = line[position+length:].find(';')
        if endPosition <= 0:
            endPosition = len(line[position+length:])
        return line[position+length:position+length+endPosition].strip()