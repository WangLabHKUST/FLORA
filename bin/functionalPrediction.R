# Rscript functionalPrediction.R .../lnc.info.txt .../coding.info.txt .../lnc.FPKM.txt .../coding.FPKM.txt  cor.method  [.../ARACNe_network(noARACNe)] output_dir LINC [CorCutoff] [CorAdjpCutoff] [GOCutoff] [PorN]

if (!require('getopt')) install.packages('getopt')
if (!require('foreach')) install.packages('foreach')
if (!require('doParallel')) install.packages('doParallel')
if (!require('energy')) install.packages('energy')
if (!require('gProfileR')) install.packages('gProfileR')
if (!require('ggplot2')) install.packages('ggplot2')

library(getopt)
library(foreach)
library(doParallel)
library(energy)
library(gProfileR)
library(ggplot2)

command <- matrix(c( 
  'Help', 'h', 0, 'loical', 'help',
  'LncInfo', 'a', 1, 'character', 'information of lncRNAs, including gene.id and gene.name',
  'CodingInfo', 'b', 1, 'character', 'information of coding genes, including gene.id and gene.name',
  'LncExpr', 'c', 1, 'character', 'expression of lncRNAs',
  'CodingExpr', 'd', 1, 'character', 'expression of coding genes',
  'CorMethod', 'm', 1, 'character', 'correlation method,  method = c("pearson", "spearman", "distance")',
  'ARACNeNetwork', 'n', 2, 'character', 'output network of ARACNe, optional',
  'OutputDir', 'o', 1, 'character', 'output directory',
  'LINC', 'i', 1, 'character', 'the name(s) of your lncRNA(s) of interest, separated by ","',
  'CorCutoff', 'x', 2, 'numeric', 'cutoff for correlation coefficient. by default use = 0',
  'CorAdjpCutoff', 'y', 2, 'numeric', 'cutoff for correlation adjust p-value. by default use = 0.001',
  'GOCutoff', 'z', 2, 'numeric', 'cutoff for GO terms p-value. by default use = 0.001',
  'PorN', 'p', 2, 'character', 'use only "positive", "negative" or "all" related coding genes to do the prediction. by default use = "positive"'), byrow=T, ncol=5)
args <- getopt(command)

if (!is.null(args$Help) || is.null(args$LncInfo) || is.null(args$CodingInfo) || is.null(args$LncExpr) || is.null(args$CodingExpr) ||
   ((args$CorMethod %in% c("pearson", "spearman", "distance")) == F) || is.null(args$OutputDir) || is.null(args$LINC)  ) {
  cat(paste(getopt(command, usage = T), "\n"))
  q(status=1)
}


# Load gene information
lnc.info <- read.table(args$LncInfo, sep = "\t", header = T)
coding.info <- read.table(args$CodingInfo, sep = "\t", header = T)

# load expression data
lnc.expr <- read.table(args$LncExpr, sep = "\t", header = T, row.names = 1)
coding.expr <- read.table(args$CodingExpr, sep = "\t", header = T, row.names = 1)
lnc.expr <- log2(lnc.expr + 1)
coding.expr <- log2(coding.expr + 1)

# correlation method,  method = c("pearson", "spearman", "distance")
cor.method <- args$CorMethod

# Load gene regulatory network constructed by ARACNe-AP, optional
ARACNe_network <- args$ARACNeNetwork
if (!is.null(ARACNe_network)){
  ARACNe <- read.table(ARACNe_network, sep = "\t", header = T)
}

# output_dir
dir.create(args$OutputDir)
output_dir <- args$OutputDir

# Input the name(s) of your lncRNA(s) of interest
lnc.name.pool <- strsplit(args$LINC, split = ",")[[1]]
lnc.name.pool <- as.data.frame(lnc.name.pool)
colnames(lnc.name.pool) <- "name"
for (i in 1:nrow(lnc.name.pool)){
  lnc.name.pool[i, "id"] <- lnc.info[which(lnc.info$name == as.character(lnc.name.pool[i, "name"]))[1], "id"]
}

lnc.info <- subset(lnc.info, lnc.info$id %in% as.character(lnc.name.pool$id))
lnc.expr <- lnc.expr[as.character(lnc.name.pool$id), ]

# set cutoff 
if (is.null(args$CorCutoff)) {
  CorCutoff <- 0
} else {
  CorCutoff <- args$CorCutoff
}
if (is.null(args$CorAdjpCutoff)) {
  CorAdjpCutoff <- 0.001
} else {
  CorAdjpCutoff <- args$CorAdjpCutoff
}
if (is.null(args$GOCutoff)) {
  GOCutoff <- 0.001
} else {
  GOCutoff <- args$GOCutoff
}

# use only "positive", "negative" or "all" related coding genes
if (is.null(args$PorN)){
  PorN <- "positive"
} else {
  PorN <- args$PorN
}

#########################################################################################################
if (cor.method %in% c("pearson", "spearman")){
  # calculate correlation cor
  cores <- detectCores(logical=T)
  cl <- makeCluster(cores)
  registerDoParallel(cl,cores = cores)
  cor <- foreach(i=1:nrow(coding.expr), .combine = 'rbind', .export = c("coding.expr", "lnc.expr", "cor.method") ) %dopar%
  { 
    out <- matrix(nrow = 1, ncol = nrow(lnc.expr))
    out <- data.frame(out)
    rownames(out) <- rownames(coding.expr)[i]
    colnames(out) <- rownames(lnc.expr)
    for (j in 1:nrow(lnc.expr)){
      out[1,j] <- cor(as.numeric(coding.expr[i,]), as.numeric(lnc.expr[j,]), method = cor.method)
    }
    out
  }
  stopImplicitCluster()
  stopCluster(cl)
  
  # calculate correlation cor.p
  cores <- detectCores(logical=T)
  cl <- makeCluster(cores)
  registerDoParallel(cl,cores = cores)
  cor.p <- foreach(i=1:nrow(coding.expr), .combine = 'rbind', .export = c("coding.expr", "lnc.expr", "cor.method") ) %dopar%
  { 
    out <- matrix(nrow = 1, ncol = nrow(lnc.expr))
    out <- data.frame(out)
    rownames(out) <- rownames(coding.expr)[i]
    colnames(out) <- rownames(lnc.expr)
    for (j in 1:nrow(lnc.expr)){
      out[1,j] <- cor.test(as.numeric(coding.expr[i,]), as.numeric(lnc.expr[j,]), method = cor.method)$p.value
    }
    out
  }
  stopImplicitCluster()
  stopCluster(cl)
  
} else if (cor.method == "distance"){
  # calculate correlation cor
  cores <- detectCores(logical=T)
  cl <- makeCluster(cores)
  registerDoParallel(cl,cores = cores)
  cor <- foreach(i=1:nrow(coding.expr), .combine = 'rbind', .export = c("coding.expr", "lnc.expr", "cor.method"), .packages = "energy" ) %dopar%
  { 
    out <- matrix(nrow = 1, ncol = nrow(lnc.expr))
    out <- data.frame(out)
    rownames(out) <- rownames(coding.expr)[i]
    colnames(out) <- rownames(lnc.expr)
    for (j in 1:nrow(lnc.expr)){
      out[1,j] <- dcor(as.numeric(coding.expr[i,]), as.numeric(lnc.expr[j,]))
    }
    out
  }
  stopImplicitCluster()
  stopCluster(cl)
  
  # calculate correlation cor.p
  cores <- detectCores(logical=T)
  cl <- makeCluster(cores)
  registerDoParallel(cl,cores = cores)
  cor.p <- foreach(i=1:nrow(coding.expr), .combine = 'rbind', .export = c("coding.expr", "lnc.expr", "cor.method"), .packages = "energy" ) %dopar%
  { 
    out <- matrix(nrow = 1, ncol = nrow(lnc.expr))
    out <- data.frame(out)
    rownames(out) <- rownames(coding.expr)[i]
    colnames(out) <- rownames(lnc.expr)
    for (j in 1:nrow(lnc.expr)){
      out[1,j] <- dcor.ttest(as.numeric(coding.expr[i,]), as.numeric(lnc.expr[j,]), distance = F)$p.value
    }
    out
  }
  stopImplicitCluster()
  stopCluster(cl)
  
}

#######################################################
for (l in 1:nrow(lnc.name.pool)){
  lnc.name <- as.character(lnc.name.pool[l, "name"])
  lnc.id <- as.character(lnc.name.pool[l, "id"])
  
  mydata <- matrix(nrow = nrow(cor), ncol = 7)
  mydata <- as.data.frame(mydata)
  colnames(mydata) <- c("lnc.id", "lnc.name", "coding.id", "coding.name", "cor", "cor.p", "cor.adjp")
  mydata$lnc.id <- lnc.id
  mydata$lnc.name <- lnc.name
  mydata$coding.id <- rownames(cor)
  mydata$cor <- cor[, lnc.id]
  mydata$cor.p <- cor.p[, lnc.id]
  mydata$cor.adjp <- mydata$cor.p * nrow(mydata)
  
  mydata <- subset(mydata, abs(mydata$cor) > CorCutoff)
  mydata <- subset(mydata, mydata$cor.adjp < CorAdjpCutoff)
  if (PorN == "negative"){
    mydata <- subset(mydata, mydata$cor < 0)
  } else if (PorN == "all"){
    mydata <- mydata
  } else {
    mydata <- subset(mydata, mydata$cor > 0)
  }
  
  ## ARACNe_network
  if (!is.null(ARACNe_network)){
    ARACNe.tmp <- subset(ARACNe, (ARACNe$Gene1.id == lnc.id) | (ARACNe$Gene2.id == lnc.id))
    ARACNe.tmp <- subset(ARACNe.tmp, (ARACNe.tmp$Gene1.id %in% as.character(mydata$coding.id)) | (ARACNe.tmp$Gene2.id %in% as.character(mydata$coding.id)))
    mydata <- subset(mydata, mydata$coding.id %in% c(as.character(ARACNe.tmp$Gene1.id), as.character(ARACNe.tmp$Gene2.id)))
    if (nrow(ARACNe.tmp) > 0){
      for (i in 1:nrow(ARACNe.tmp)) {
        pcg.id <- setdiff(c(as.character(ARACNe.tmp[i, "Gene1.id"]), as.character(ARACNe.tmp[i, "Gene2.id"])), lnc.id)
        mydata[which(mydata$coding.id == pcg.id), "MI"] <- ARACNe.tmp[i, "MI"]
        mydata[which(mydata$coding.id == pcg.id), "ARACNe.pvalue"] <- ARACNe.tmp[i, "pvalue"]
      }
    }
  }
  
  ## GO enrichment
  GO <- NULL
  if (nrow(mydata) > 0) {
    for (i in 1:nrow(mydata)){
      mydata[i, "coding.name"] <- as.character(coding.info[which(coding.info$id == as.character(mydata[i, "coding.id"])), "name"])
    }
    if (is.null(ARACNe_network)) {
      mydata <- mydata[order(abs(mydata$cor), decreasing = T), ]
    } else {
      mydata <- mydata[order(mydata$MI, decreasing = T), ]
    }
    GO <- gprofiler(as.character(mydata$coding.name), organism = "hsapiens", ordered_query = T, hier_filtering="strong")
    GO <- GO[order(GO[,"p.value"]),]
    GO <- subset(GO, p.value < GOCutoff)
  }
  
  write.table(mydata, paste0(output_dir, "/", lnc.name, ".net.txt"), sep="\t", row.names = F, quote = F)
  write.table(GO, paste0(output_dir, "/", lnc.name, ".GO.txt"), sep="\t", row.names = F, quote = F)
  
  ## draw plot
  if (!is.null(GO)){
    GO <- subset(GO, domain %in% c("BP", "CC", "MF", "keg", "rea", "mir", "cor"))
    if (nrow(GO) > 0) {
      GO[which(GO$domain == "BP"), "tmp.rank"] <- 1
      GO[which(GO$domain == "CC"), "tmp.rank"] <- 2
      GO[which(GO$domain == "MF"), "tmp.rank"] <- 3
      GO[which(GO$domain == "keg"), "tmp.rank"] <- 4
      GO[which(GO$domain == "rea"), "tmp.rank"] <- 5
      GO[which(GO$domain == "mir"), "tmp.rank"] <- 6
      GO[which(GO$domain == "cor"), "tmp.rank"] <- 7
      
      GO.clean <- NULL
      for (i in unique(GO$term.name)){
        GO.clean <- rbind(GO.clean, subset(GO, GO$term.name == i)[1,])
      }
      GO.clean <- GO.clean[order(GO.clean[,"tmp.rank"]),]
      GO.clean$term.name = factor(GO.clean$term.name, levels = rev(GO.clean$term.name))
      GO.clean$domain = factor(GO.clean$domain, levels = c("BP", "CC", "MF", "keg", "rea", "mir", "cor"))
      
      p <- ggplot(aes(x = GO.clean$term.name, y = -log10(GO.clean$p.value), fill = domain), data = GO.clean) + 
        geom_bar(stat="identity", position="dodge") +
        labs(x = "", y = "-log10_Pvalue", title= lnc.name )  + coord_flip() +
        theme_classic() + scale_fill_manual(values = c("BP" = '#8470FF', "CC" = '#87CEFA',  "MF" = '#FFC125', 
                                                       "keg" = "#66c2a5",  "rea"= "#fc8d62", "mir" = "#8da0cb",
                                                       "cor" = "#e78ac3")) +
        theme(legend.title=element_blank(), axis.title.x = element_text(size=9), plot.title = element_text(color = "black", face = "italic") )
      ggsave(file = paste0(output_dir, "/", lnc.name, ".GO.pdf"), p)
    }
  }
  
}


