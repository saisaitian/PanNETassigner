PanNETsubtypesmicRNAVal <-function(Indata,IndSAMdata)
{
  # Selecting non-redudant genes
  res1<-screenExpr(Indata,0)
  cat("Data after screenExpr has", ncol(res1)-1, "samples", "\n", append = T)
  sdv<-paste("_sd0.txt",sep="")
  InputFileName_sd<-paste("241210-Weidenmann-49-samples-data-genesonly_1_sd0_removed_non_PNET_few_Mets",sdv,sep="")
  
  # Creating storage for results and setting the working directory
  dir.create("PanNETassignerResults", showWarnings = FALSE)
  dir.create("PanNETassignerResults/micRNAVal", showWarnings = FALSE)
  dir<-getwd()   
  setwd(paste0(dir,"/PanNETassignerResults/micRNAVal"))
  
  # Selecting SAM genes  sam<-IndSAMdata
  sam<-IndSAMdata
  sam_mat<-match(sam[,1],res1[,1])
  sam_w<-which(!is.na(sam_mat))
  data_sam<-res1[sam_mat[sam_w],]
  InputFileName_sam<-paste((strsplit(InputFileName_sd,".txt")[1]),"_SAM.txt",sep="")
  write.table(data_sam,InputFileName_sam, row.names=FALSE,sep="\t",quote=FALSE)
  
  dataUnique<-apply(data.matrix(data_sam[,2:dim(res1)[2]]),2,as.numeric)
  rownames(dataUnique)<-data_sam[,1]
  
  # Median centering data
  rowMed <- apply(dataUnique,1,median, na.rm=TRUE);
  dataUnique<-dataUnique-rowMed
  dataUnique_1<-cbind(rownames(dataUnique),dataUnique)
  InputFileName_sam_med<-paste((strsplit(InputFileName_sam,".txt")[1]),"_data_cdt.txt",sep="")
  write.table(dataUnique_1,InputFileName_sam_med, row.names=FALSE,sep="\t",quote=FALSE)
  
  # NMF analysis
  cat("Performing NMF analysis...")
  nmfconsensus(input.ds=InputFileName_sam_med,k.init=2,k.final=4,num.clusterings=20,maxniter=500, error.function="euclidean",data=dataUnique)
  
  k=1:5
  # Ordering samples by NMF consensus k=3 cluster
  con<-read.delim("consensus.k.3.gct", stringsAsFactors=FALSE)
  m<-match(con[,1],colnames(dataUnique_1))
  #m
  
  InputFileName_sam_med_NMF<-paste((strsplit(InputFileName_sam_med,".txt")[1]),"_NMF_ordered_3.txt",sep="")
  data_2<-cbind(dataUnique_1[,1],dataUnique_1[,m])
  write.table(data_2,InputFileName_sam_med_NMF, row.names=FALSE,sep="\t",quote=FALSE)
  setwd(dir)
}
