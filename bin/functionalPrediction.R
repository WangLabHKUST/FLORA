getnetwork <- function(lnc.info, coding.info, network, lnc.name) {
  j <- which(lnc.info$name == lnc.name)
  lnc.id <- as.character(lnc.info[j,1]) 
  
  tmp <- subset(network, (Regulator == lnc.id) | (Target == lnc.id) )
  tmp <- data.frame(tmp)
  
  tmp.name <- matrix(nrow = nrow(tmp), ncol = ncol(tmp))
  for (i in 1: nrow(tmp)) {
    for (j in 1:2) {
      k <- which(lnc.info$id == as.character(tmp[i,j]) )
      if (length(k)>0) {
        tmp.name[i,j] <- as.character(lnc.info[k,2]) 
      } else {
        k <- which(coding.info$id == as.character(tmp[i,j]) )
        tmp.name[i,j] <- as.character(coding.info[k,2]) 
      }
    }
  }
  tmp.name[,3] <- tmp[,3]
  tmp.name[,4] <- tmp[,4]
  colnames(tmp.name) <- colnames(network)
  
  lnc.coding <- tmp.name
  for ( i in 1: nrow(lnc.coding)){
    k1 <- which(lnc.info$name == as.character(lnc.coding[i,1]) )
    if (length(k1)>0){
      k2 <- which(lnc.info$name == as.character(lnc.coding[i,2]) )
      if (length(k2) >0 ){
        lnc.coding[i,3] <- "-"
      }
    }
  }
  od <- data.frame(order(lnc.coding[,3],decreasing = T))
  lnc.coding <-  lnc.coding[od[,1],]
  lnc.coding <-  subset(lnc.coding, lnc.coding[,3] != "-")
  
  return(lnc.coding)
}


makePrediction <- function(lnc.name, lnc.coding, gotype="all") {
  library(gProfileR)
  
  gene.regulator <- setdiff(lnc.coding[,1], lnc.name)
  gene.target <- setdiff(lnc.coding[,2], lnc.name)
  gene.all <- setdiff(union(lnc.coding[,1], lnc.coding[,2]), lnc.name)
  
  if (gotype == "regulator") {
    gene.use <- gene.regulator
  } else if (gotype == "target") {
    gene.use <- gene.target
  } else {
    gene.use <- gene.all
  }
  
  GO <- gprofiler(gene.use, organism = "hsapiens",ordered_query = T,hier_filtering="strong")
  GO <- GO[order(GO[,3]),]
  
  BP <- subset(GO, domain == "BP")
  CC <- subset(GO, domain == "CC")
  MF <- subset(GO, domain == "MF")
  
  return(list(gene.use=gene.use, GO=GO, BP=BP, CC=CC, MF=MF))
}

