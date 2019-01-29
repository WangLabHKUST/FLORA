# FLORA: Functional Long-noncoding RNA Assembly Workflow

# Ownership
[Wang Lab at HKUST](http://wang-lab.ust.hk)


# Status
Active Development


# Introduction
Fast Long-noncoding RNA Assembly Workflow (FLORA) is a collection of easy-to-use command line tools for fast lncRNA transcriptome assembly from RNA-seq BAM files. It was initially developed by [Hongyu Shi](https://github.com/AlexHelloWorld) and is still under active development by other members of [Wang Lab at HKUST](http://wang-lab.ust.hk).


## Prerequisites

• bedtools (https://bedtools.readthedocs.io/en/latest/)

• CPAT (http://rna-cpat.sourceforge.net/)

• (optional) CPAT prebuilt logit model and hexamer table files (https://sourceforge.net/projects/rna-cpat/files/v1.2.2/prebuilt_model/)

• GTF files from Gencode (https://www.gencodegenes.org/)

• GTF files and assembly report from Refseq (ftp://ftp.ncbi.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/)

• gProfilerR package in R (https://cran.r-project.org/package=gProfileR)

## Installing

Download/Clone the repository and change the working directory to the FLORA directory

```
python setup.py install
```


## Run

### generateFilteredBams.py

Generate BAM files with certain regions removed.

```
usage: generateFilteredBams [-h] -g INPUTGTF [-t TYPES [TYPES ...]]
                            [-o OUTPUTDIR] [-n NTHREAD]
                            inputBams

positional arguments:
  inputBams             The path to a txt file that contains paths to input
                        bam files separated by the newline character.

optional arguments:
  -h, --help            show this help message and exit
  -g INPUTGTF           The path to reference gene annotation file in GTF
                        format. It can be downloaded from gencode website.
  -t TYPES [TYPES ...]  The gene types to be removed from BAM files.
  -o OUTPUTDIR          Output directory for the output files
  -n NTHREAD            Number of threads to be used for bedtools intersect.
```

### filterTranscripts.py

Identify lncRNA transcripts from transcriptome assembled by Stringtie or Cufflinks.

```
usage: filterTranscript [-h] -r REFERENCE [-e EXON] [-l LENGTH] [-c CPAT] -x
                        HEXAMER -m LOGIT [-o OUTPUTGTF]
                        inputGTF

positional arguments:
  inputGTF      Input transcriptome file in GTF format. Can be output from
                stringtie or cufflinks

optional arguments:
  -h, --help    show this help message and exit
  -r REFERENCE  The path to reference genome sequences in FASTA format. It
                will be indexed automatically by CPAT if .fai file is not
                present.
  -e EXON       The least number of exons a transcipt should have in order to
                be kept in the final transcriptome
  -l LENGTH     The shortest transcript to be kept in the final transcriptome
  -c CPAT       CPAT cutoff used to filter transcripts
  -x HEXAMER    The path to hexamer table required by CPAT. Can be downloaded
                from CPAT website.
  -m LOGIT      The path to logit model required by CPAT. Can be downloaded
                CPAT website.
  -o OUTPUTGTF  Output prefix for the final transcriptome GTF file
```

### annotateTranscripts.py

Find overlapped and nearby genes of lncRNAs in reference annotation (RefSeq or GENCODE annotation in GFF format).

```
usage: annotateTranscripts [-h] -r REFERENCE -f PATH [-i IDENTIFIER]
                           [-d DISTANCE] [-n NUMBER] [-o OUTPUT]
                           inputGTF

positional arguments:
  inputGTF       Input transcriptome file in GTF format. The file could be
                 output from Cufflinks or StringTie

optional arguments:
  -h, --help     show this help message and exit
  -r REFERENCE   Identify the source of reference annotation. Currently, it
                 should be either refseq or gencode
  -f PATH        The path to the reference annotation file.
  -i IDENTIFIER  The path to the RefSeq assembly report. The file is used to
                 replace RefSeq sequence identifiers with USCS identifiers.
                 Required if you use RefSeq annotation.
  -d DISTANCE    Specify the distance within which the nearby genes should
                 locate. Default: 10,000 bp. Unit: bp.
  -n NUMBER      Specify the number of nearby genes obtained. Default: 1.
  -o OUTPUT      Output file name for the final annotation. Default:
                 FLORA_annotation.txt
```

### functionalPrediction


Gene regulatory network is constructed via "ARACNe-AP" (https://github.com/califano-lab/ARACNe-AP) based on expression data.
[Command Line]
```
#########  Example code for ARACNe-AP #########

# Perform 3 bootstrapping and network reconstruction
for i in {1..3}
do
  (
    java -Xmx120G -jar .../aracne.jar -e .../expression.matrix.txt  -o output/ --tfs .../lnc_tf.txt --pvalue 1E-8 --seed $i
  ) &
  wait
done

# Generating consensus network
java -Xmx120G -jar .../aracne.jar -o output/ --consolidate

```

LncRNA's function was predicted based on gene regulatory network.
[R]
```
getnetwork()
    usage: getnetwork(lnc.info, coding.info,         # lncRNA and coding genes' information (id and name)
                      network,                       # gene regulatory network (predicted by ARACNe-AP)
                      lnc.name)                      # name of target lncRNA   
    output: coding - target lncRNA  regulatory network

makePrediction()
    usage: makePrediction(lnc.name,           # name of lncRNA 
                          lnc.coding,         # coding - target lncRNA  regulatory network, generated in getnetwork() function
                          gotype)             # ["regulator","target","all"] use regulator/target/all genes in the netowrk to do the prediction
    output: list of GO terms

```

[Example]
```
######### Example code for predicting the function of LINC01614 (your lncRNA of interest) #########

[Command Line]
Rscript example.R ../bin/functionalPrediction.R lnc.info.txt coding.info.txt network.txt output_dir LINC01614 

outputs:
  LINC01614_lnc.coding.txt   # table of genes connected with your lncRNA of interest
  LINC01614_GO.txt           # table of significant GO terms associated with your lncRNA of interest
  LINC01614_GO.pdf           # figure of significant GO terms associated with your lncRNA of interest

```
![image](https://github.com/shuangat/FLORA/blob/master/data/LINC01614.png?raw=true)

