# FLORA: Functional Long-noncoding RNA Assembly Workflow


## Ownership
[Wang Lab at HKUST](http://wang-lab.ust.hk)


## Status
Active Development


## Introduction
Fast Long-noncoding RNA Assembly Workflow (FLORA) is a collection of easy-to-use command line tools for fast lncRNA transcriptome assembly from RNA-seq BAM files. It was initially developed by [Hongyu Shi](https://github.com/AlexHelloWorld) and is still under active development by other members of [Wang Lab at HKUST](http://wang-lab.ust.hk).


## Prerequisites

• bedtools (https://bedtools.readthedocs.io/en/latest/)

• CPAT (http://rna-cpat.sourceforge.net/)

• (optional) CPAT prebuilt logit model and hexamer table files (https://sourceforge.net/projects/rna-cpat/files/v1.2.2/prebuilt_model/)

• GTF files from Gencode (https://www.gencodegenes.org/)

• GTF files and assembly report from Refseq (ftp://ftp.ncbi.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/)

• gProfilerR package in R (https://cran.r-project.org/package=gProfileR)


## Install

Download/Clone the repository and change the working directory to the FLORA directory

```
python setup.py install
```


## Workflow

### Step 1: Generate Filtered BAM

Starting from bam files, **generateFilteredBams.py** generates prefiltered bam files by automatically removing reads with low mapping quality and reads mapped to ribosomal RNAs (built in for hg38, hg19, mm10, mm9, danRer11, danRer10, dm6). In addition, user-defined list of gene types (such as protein-coding, pseudogenes, immunoglobin and T-cell receptor genes) can also be filtered out from the bam files.

```
usage: generateFilteredBams [-h] -g INPUTGTF [-t TYPES [TYPES ...]]
                            [-o OUTPUTDIR] [-n NTHREAD]
                            inputBams

positional arguments:
  inputBams             The path to a txt file that contains paths to input
                        bam files separated by the newline character.

optional arguments:
  -h, --help            show this help message and exit.
  -g INPUTGTF           The path to reference gene annotation file in GTF
                        format. It can be downloaded from gencode website.
  -t TYPES [TYPES ...]  The gene types to be removed from BAM files.
  -o OUTPUTDIR          Output directory for the output files.
  -n NTHREAD            Number of threads to be used for bedtools intersect.

```
Example:
```
generateFilteredBams.py -g reference.gtf -t protein_coding inputBams.txt
```

### Step 2: Transcriptome Assembly and Merging

Transcriptome assembly can be achieved with **StringTie** from the preprocessed RNA-seq data. Since this step can be slow on large datasets, we recommend running StringTie on high-performance clusters in parallel (example code provided below).
Assembled transcripts in each sample were selected for merging by **Stringtie merge** with customized parameters (parameter settings in our manuscript: (A) longer than 200 nucleotides; (B) with expression level over 0.1 TPM and 0.1 FPKM; (C) account for over 10% of all isoforms from the same loci).

Example:
```
# StringTie
stringtie clean.bam -o clean.bam.gtf -G reference.long_noncoding_RNAs.gtf

# StringTie Merge: merge multiple clean.bam.gtf to one merge.gtf
stringtie --merge -m 200 -F 0.1 -T 0.1 -f 0.1 -o merge.gtf clean.bam.gtf.list
```

### Step 3: Filtering for LncRNAs

**filterTranscripts.py** identify prospective lncRNAs from the assembled transcriptome. Transcripts are commonly selected as prospective lncRNAs by the following criteria: (A) longer than 200 nucleotides; (B) containing two or more exons; (C) coding potential score larger than 0.364 as predicted by Coding Potential Assessment Tool (CPAT).

```
filterTranscripts.py -r reference.fa -x Human_Hexamer.tsv -m Human_logitModel.RData -o FLORA_out merge.gtf

usage: filterTranscript [-h] -r REFERENCE [-e EXON] [-l LENGTH] [-c CPAT] -x
                        HEXAMER -m LOGIT [-o OUTPUTGTF]
                        inputGTF

positional arguments:
  inputGTF      Input transcriptome file in GTF format. Can be output from
                stringtie or cufflinks.

optional arguments:
  -h, --help    show this help message and exit.
  -r REFERENCE  The path to reference genome sequences in FASTA format. It
                will be indexed automatically by CPAT if .fai file is not
                present.
  -e EXON       The least number of exons a transcipt should have in order to
                be kept in the final transcriptome. Default: 2.
  -l LENGTH     The shortest transcript to be kept in the final transcriptome.    
                Default: 200.
  -c CPAT       CPAT cutoff used to filter transcripts. Default: 0.364
  -x HEXAMER    The path to hexamer table required by CPAT. Can be downloaded
                from CPAT website.
  -m LOGIT      The path to logit model required by CPAT. Can be downloaded
                CPAT website.
  -o OUTPUTGTF  Output prefix for the final transcriptome GTF file.
```

Example:
```
filterTranscripts.py -r reference.fa -x Human_Hexamer.tsv -m Human_logitModel.RData -o lncRNA.gtf merge.gtf
```

### Step 4: LncRNA Annotation

**AnnotateNovelLncRNA.sh** utilizes cuffcompare to compare the prospective lncRNAs with reference annotation (downloadable from Gencode, RefSeq, Ensembl; allowed up to 3 GTF files in the comparison), and lncRNA-expressing loci with no overlap with known genomic features were defined as novel. A detailed report of each lncRNA loci and its overlapping genes, as well as a GTF containing only novel lncRNA loci are generated.

Example:
```
./AnnotateNovelLncRNA.sh lncRNA.gtf Gencode.gtf Ensembl.gtf RefSeq.gtf
```

### Step 5: Functional Prediction

FLORA predicted the potential functions of lncRNA(s) based on expression data.  **functionalPrediction.R** provides different correlation methods(Pearson, Spearman, Distance) to construct network, which can also further combines with the gene transcriptional network. Later, all coding genes with significant associations with the target lncRNA are selected to perform Gene Ontology (GO) enrichment analysis with g-profiler. This step will output the list of genes with significant associations with the target lncRNA (LINCX.net.txt) and the enriched GO terms of the genes associated with the target lncRNA (LINCX.GO.txt and LINCX.GO.pdf).

Optional, gene transcriptional network can be constructed by **ARACNe-AP** software (https://github.com/califano-lab/ARACNe-AP) based on normalized expression data. To obtain a stable network, 100 reproducible bootstraps were performed and consolidated (time-consuming step, we recommend network construction in parallel on high-performance clusters; example code provided below).

Example code for ARACNe-AP:
```
# Perform 100 bootstrapping and network reconstruction
for i in {1..100}
do
  (
    java -Xmx120G -jar .../aracne.jar -e .../expression.matrix.txt  -o output/ --tfs .../lnc_tf.txt --pvalue 1E-8 --seed $i
  ) &
  wait
done

# Generating consensus network
java -Xmx120G -jar .../aracne.jar -o output/ --consolidate
```

To analyze the genes associated with the lncRNAs in gastric cancer, our example code is provided below to automatically output the significantly associated genes and GO terms. Additionally, the gene transcription network constructed by ARACNe-AP based on the normalized expression of all coding genes and lncRNAs in gastric cancer is provided (example/data/ARACNe_network.txt).

Example code for predicting the function of "H19" and "LINC01614" (your lncRNA(s) of interest):
```
Rscript functionalPrediction.R -a example/data/lnc.info.txt -b example/data/coding.info.txt -c example/data/lnc.FPKM.txt -d example/data/coding.FPKM.txt -m spearman  -o output -i H19,LINC01614 -x 0.5 -n data/ARACNe_network.txt

Usage: functionalPrediction.R [-h] -a LncInfo -b CodingInfo -c LncExpr -d CodingExpr -o OutputDir -i LINC1,LINC2,LINC3,... -m CorMethod [-n ARACNeNetwork] [-x CorCutoff] [-y CorAdjpCutoff] [-z GOCutoff] [-p PorN]

arguments:
    -a|--LncInfo          information of lncRNAs, including gene.id and gene.name
    -b|--CodingInfo       information of coding genes, including gene.id and gene.name
    -c|--LncExpr          expression of lncRNAs
    -d|--CodingExpr       expression of coding genes
    -o|--OutputDir        output directory
    -i|--LINC             the name(s) of your lncRNA(s) of interest, separated by ","
    -m|--CorMethod        correlation method,  method = c("pearson", "spearman", "distance")

optional arguments:
    -h|--Help             help
    -n|--ARACNeNetwork    output network of ARACNe, optional
    -x|--CorCutoff        cutoff for correlation coefficient. by default use = 0
    -y|--CorAdjpCutoff    cutoff for correlation adjust p-value. by default use = 0.001
    -z|--GOCutoff         cutoff for GO terms p-value. by default use = 0.001
    -p|--PorN             use only "positive", "negative" or "all" related coding genes to do the prediction. by default use = "positive"

outputs:
  LINCX.net.txt          # table of genes connected with LINCX
  LINCX.GO.txt           # table of significant GO terms associated with LINCX
  LINCX.GO.pdf           # figure of significant GO terms associated with LINCX
```
<div align=center><img width="500" height="450" src="https://github.com/WangLabHKUST/FLORA/blob/shuangat/example/output/LINC0164.GO.png"/></div>



21 Jan 2020
