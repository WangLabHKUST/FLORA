# Rscript example.R .../functionalPrediction.R .../lnc.info.txt .../coding.info.txt .../network.txt output_dir LINC01614 
args<-commandArgs(T)

dir.create(args[5])
output_dir <- args[5]
source(args[1])

# Load information of all coding genes and lncRNAs
lnc.info <- read.table(args[2], sep = "\t", head=T, row.names=1)
coding.info <- read.table(args[3], sep = "\t", head=T, row.names=1)

# Load gene regulatory network constructed by ARACNe-AP
network <- read.table(args[4], sep = "\t", head=T)

# Input the name of your lncRNA of interest
lnc.name <- args[6]

# Extract genes connected with the lncRNA and perform functional prediction
lnc.coding <- getnetwork(lnc.info, coding.info, network, lnc.name)
results <- makePrediction(lnc.name, lnc.coding)
gene.use <- results$gene.use
GO <- subset(results$GO, domain %in% c("BP", "CC", "MF"))
write.table(lnc.coding, paste(output_dir,"/",lnc.name,"_lnc.coding.txt",sep=""), sep="\t", row.names = FALSE)
write.table(GO, paste(output_dir,"/",lnc.name,"_GO.txt",sep=""), sep="\t", row.names = FALSE)

# Plotting significant GO terms associated with the networks of the lncRNA
library(ggplot2)
pdf( file=paste(output_dir,"/",lnc.name,"_GO.pdf",sep = "")) 
GO$term.name = factor(GO$term.name, levels = GO$term.name)
ggplot(aes(x = term.name, y = -log10(GO$p.value), fill = domain), data = GO) + 
  geom_bar(stat="identity",position="dodge") +
  labs(x = "", y = "-log10_Pvalue", title= lnc.name )  + coord_flip() +
  theme_classic() + scale_fill_manual(values = c('#8470FF','#87CEFA','#FFC125')) +
  theme(legend.title=element_blank(), axis.title.x = element_text(size=9) )
dev.off()
