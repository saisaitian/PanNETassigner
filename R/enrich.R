enrich <-function(labels1,labels2,outfile)
{
  #filename1 = "data/label_miR_subtypes_normal_removed.txt"
  #filename2 = "data/label_mRNA_subtypes_normal_removed.txt"
  # outfile PROVIDE YOUR CHOICE OF NAME
  
  ############################
  ### Reading the datasets ###
  ############################
  cor.map<-labels1
  cat("labels1 dataset has", ncol(cor.map),"genes and with",rownames(cor.map),"as row names.\n")
  
  class<-labels2
  cat("Subtype dataset has", ncol(class),"genes.\n")
  
  # Creating storage for results and setting the working directory
  dir.create("PanNETassignerResults", showWarnings = FALSE)
  dir.create("PanNETassignerResults/Enrichment", showWarnings = FALSE)
  dir<-getwd()  
  setwd(paste0(dir,"/PanNETassignerResults/Enrichment"))
  
  #####################################################
  ### Enrichment analysis with hypergeometric test for genes. ###
  #####################################################
  m3<-match(rownames(cor.map),class[,1])
  w3<-which(!is.na(m3))
  z3<-m3[w3]
  idcor.map<-cor.map[,ncol(cor.map)]
  idclass<-class[z3,ncol(class)]
  if(length(z3)<length(m3))
  {
    m33<-match(class[z3,1],rownames(cor.map))
    w33<-which(!is.na(m33))
    z33<-m33[w33]
    idcor.map<-cor.map[z33,ncol(cor.map)]
  }
  b<-table(idclass,idcor.map)
  
  pval<-matrix(NA,length(table(idclass)),length(table(idcor.map)))
  for(i in 1:length(table(idcor.map))) {
    for(j in 1:length(table(idclass))) {
      pval[j,i]<-phyper((b[j,i] - 1),sum(b[,i]),sum(b)-sum(b[,i]),sum(b[j,]),lower.tail=FALSE)
    }
  }
  adj.pval<-matrix(p.adjust(c(pval),"fdr"),length(table(idclass)),length(table(idcor.map)))
  rownames(adj.pval)<-rownames(pval)<-rownames(b)
  colnames(adj.pval)<-colnames(pval)<-colnames(b)
  
  #####################################################
  ### Enrichment analysis with heatmap. ###
  #####################################################
  file3<-unlist(strsplit(outfile,".txt"))[1]
  pdf(paste(file3,"_enrichment_heatmap1.pdf",sep=""), height=13, width=13)
  heatmap(adj.pval,col=rainbow(1000,start=(0+.7*min(adj.pval)),end=.7*max(adj.pval)),Rowv=NA,Colv=NA,scale="none",
          margins=c(10,10),cexRow=1.3,cexCol=1.3,main="Hypergeometric test Enrichment analysis",xlab="Gene subtypes",ylab="miR subtypes")
  dev.off() 
  setwd(dir)
}
