import pandas as pd
import numpy as np
from sys import argv
import csv

# read merged cuffcompare reports of lncRNA transcripts
report = argv[1]
inputgtf = argv[2]
data = pd.read_table(report)

# annotate each transcript from lncRNA-expressing genes, mark overlapped genes
cs = [i for i in data.columns if i.find('ref_gene_id')!=-1]
lncgenes = sorted(list(set(data['cuff_gene_id'])))

anno = np.zeros((len(lncgenes), len(cs)))
agene = np.empty([len(lncgenes), len(cs)], dtype=object)

for i in range(len(lncgenes)):
    select = data[data['cuff_gene_id']==lncgenes[i]]
    for j in range(len(cs)):
        if list(set(select.loc[:,cs[j]]))==["-"]:
            agene[i,j]="NA"
        else:
            agene[i,j]=",".join([s for s in set(select.ref_gene_id) if s!="-"])
            anno[i,j] = 1 

df1 = pd.DataFrame(agene, columns=["GeneName_ref"+str(i+1) for i in range(len(cs))])
df2 = pd.DataFrame(anno, columns=["Overlap_ref"+str(i+1) for i in range(len(cs))])
df2['novelLoci'] = 0
for i in range(df2.shape[0]):
    if sum(df2.iloc[i,:-1])==0:
        df2.iloc[i,-1] = "Novel lncRNA"
    else:
        df2.iloc[i,-1] = "Overlap with annotated genes"

# generate a report of novel lncRNAs and lncRNAs overlapped with annotations
df = pd.concat([df1, df2], axis=1)
df.insert(0,"lncRNA_ID",pd.Series(lncgenes))
df.insert(0,"lncRNA_order",pd.Series([int(i.split(".")[1]) for i in lncgenes]))
df = df.sort_values(by="lncRNA_order")

df2 = df.drop(['lncRNA_order']+[i for i in df.columns if i.find("Overlap")!=-1], axis=1)
df2.to_csv('%s.detailed.txt' % report, index=False, sep="\t")

# generate a gtf file with novel lncRNAs only
novel = list(df[df.novelLoci=="Novel lncRNA"].lncRNA_ID)
gtf = pd.read_table(inputgtf, header=None)
gtf['geneid']=gtf[8].str.split("\"", 2, expand=True)[1]
novelgtf = gtf[gtf['geneid'].isin(novel)]
print ("number of novel lncRNA genes:", len(novel))
print ("number of novel lncRNA genes found in the GTF:", len(set(novelgtf.geneid)))
novelgtf2 = novelgtf.drop(['geneid'],axis=1)
novelgtf2.to_csv('%s.novelLnc.gtf' % inputgtf, header=None, index=False, sep="\t", 
                 quoting=csv.QUOTE_NONE, doublequote=False)

