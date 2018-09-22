
def gtfToCpatBed(inputfile, outputfile, arbitrary_strand):
    input = open(inputfile, 'r')
    output = open(outputfile, 'w')
    transcript_id = ""
    transcript_string = ""
    for line in input:
        #if the input line transcript_id is not equal to the transcript_id
        #process transcript_string
        #update transcript_id and transcript_string
        incoming_transcript_id = getTranscriptId(line)
        if incoming_transcript_id != transcript_id:
            if transcript_id == "":
                transcript_id = incoming_transcript_id
                transcript_string = line
            else:
                bed_line = processTranscriptString(transcript_string, transcript_id, arbitrary_strand)
                output.write(bed_line)
                transcript_id = incoming_transcript_id
                transcript_string = line
        else:
            transcript_string += line
    output.write(processTranscriptString(transcript_string, transcript_id, arbitrary_strand))
    output.close


def getTranscriptId(input):
    position = input.find("transcript_id")
    endPosition = input[position+15:].find('"')
    return input[position+15:position+15+endPosition]

def processTranscriptString(transcript_string, transcript_id, arbitrary_strand):
    transcript_list = transcript_string.strip().split('\n')
    exon_eg = transcript_list[0].strip().split('\t')
    chromosome = exon_eg[0]
    strand = arbitrary_strand
    exon_list = []
    for transcript in transcript_list[1:]:
        transcript_element = transcript.strip().split('\t')
        positions = (int(transcript_element[3]), int(transcript_element[4]))
        exon_list.append(positions)
    start = exon_list[0][0]
    end = exon_list[-1][1]
    length_string = ""
    exon_position_string = ""
    for exon in exon_list:
        length_string += str(exon[1]-exon[0]+1) + ','
        exon_position_string += str(exon[0] - start) + ','
    #construct output bed format line
    if chromosome[0:3] == "chr":
        outputstring = chromosome + '\t' + str(start) + '\t' + str(end) + '\t' + transcript_id + '\t0\t' + strand + '\t' + str(start)+'\t'+str(end) + '\t0,0,0\t'+str(len(transcript_list)-1)+'\t'+length_string+'\t'+exon_position_string+'\n'
        return outputstring
    else:
        return ""

def countTranscriptGTF(inputfile):
    idList = set()
    input = open(inputfile, 'r')
    for line in input:
        idList.add(getTranscriptId(line))
    return len(idList)

