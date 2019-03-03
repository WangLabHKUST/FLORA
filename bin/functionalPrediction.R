getnetwork <- function(lnc.info, coding.info, network, lnc.name) {
  j <- which(lnc.info$name == lnc.name)
  lnc.id <- as.character(lnc.info[j,1]) 
  
  tmp <- subset(network, (Gene1 == lnc.id) | (Gene2 == lnc.id) )
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


functionalPrediction <- function(lnc.name, lnc.coding) {
  library(gProfileR)
  
  lnc.coding <- data.frame(lnc.coding)
 
  gene.use <- matrix(nrow = nrow(lnc.coding), ncol = 1)
  for (i in 1:nrow(lnc.coding)) {
    gene.use[i] <- setdiff(c(as.character(lnc.coding[i,1]),as.character(lnc.coding[i,2])), lnc.name)
  }
  gene.use <- as.character(gene.use)
  
  GO <- gprofiler(gene.use, organism = "hsapiens",ordered_query = T,hier_filtering="strong", max_p_value=0.01)
  GO <- GO[order(GO[,3]),]
  
  BP <- subset(GO, domain == "BP")
  CC <- subset(GO, domain == "CC")
  MF <- subset(GO, domain == "MF")
  
  return(list(gene.use=gene.use, GO=GO, BP=BP, CC=CC, MF=MF))
}
