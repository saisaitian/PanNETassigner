PanNETsubtypesmicRNA <-function(Indata)
{
  # Selecting genes with SD > sdvalue
  sdvalue<-0.5
  res1<-screenExpr(Indata,sdvalue)
  cat("Data after screenExpr has", ncol(res1)-1, "samples", "\n", append = T)
  dataUnique<-apply(data.matrix(res1[,2:dim(res1)[2]]),2,as.numeric)
  rownames(dataUnique)<-res1[,1]
  
  # Median centering data
  rowMed <- apply(dataUnique,1,median, na.rm=TRUE);
  dataUnique<-round(dataUnique-rowMed, digits=6)
  
  # Creating storage for results and setting the working directory
  dir.create("PanNETassignerResults", showWarnings = FALSE)
  dir.create("PanNETassignerResults/micRNA", showWarnings = FALSE)
  dir<-getwd()   
  setwd(paste0(dir,"/PanNETassignerResults/micRNA"))
  
  # NMF analysis
  InputFileName_sd_med<-paste("mirna-matched-samples-mRNA_tumor_only",sdvalue,"_data_cdt.txt",sep="")
  cat("Performing NMF analysis...")
  nmfconsensus(input.ds=InputFileName_sd_med,k.init=2,k.final=5,num.clusterings=20,maxniter=500, error.function="euclidean",data=dataUnique)
  
  # Ordering samples by NMF consensus k=3 cluster
  con<-read.delim("consensus.k.3.gct", stringsAsFactors=FALSE)
  m<-match(con[,1],colnames(dataUnique))
  
  InputFileName_sd_med_NMF<-paste((strsplit(InputFileName_sd_med,".txt")[[1]])[1],"_NMF_ordered_3.txt",sep="")
  data_2<-cbind(rownames(dataUnique),dataUnique[,m])
  write.table(data_2,InputFileName_sd_med_NMF, row.names=FALSE,sep="\t",quote=FALSE)
  
  data_2_m<-matrix(as.numeric(data_2[,2:dim(data_2)[2]]),nrow=dim(data_2)[1],ncol=dim(data_2)[2]-1)
  rownames(data_2_m)<-data_2[,1]
  
  # SAM analysis
  cat("Performing SAM analysis...")
  sam.out<-sam(data_2_m,con[,2],rand=123,gene.names=rownames(data_2_m))
  write.table(sam.out@mat.fdr,"samoutput.txt",row.names=FALSE,sep="\t",quote=FALSE)
  
  # Selecting SAM genes with delta value "1.3", which provides 30 miRs
  g<-summary(sam.out, 1.3)@row.sig.genes
  data_2_sam<-data_2[g,]
  
  outsam<-paste("_SAM_",length(g),".txt",sep="")
  OutputFileName_sd_med_NMF_SAM<-paste((strsplit(InputFileName_sd_med_NMF,".txt")[[1]])[1],outsam,sep="")
  write.table(data_2_sam,OutputFileName_sd_med_NMF_SAM,row.names=FALSE,sep="\t",quote=FALSE)
  sessionInfo()
  save.image("PanNETsubtypesmicRNA.RData")
  setwd(dir)
}
