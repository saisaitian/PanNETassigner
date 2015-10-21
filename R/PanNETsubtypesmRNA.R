PanNETsubtypesmRNA <-function(Indata)
{
  # Selecting genes with SD > 0.8
  res1<-screenExpr(Indata,0.8)
  cat("Data after screenExpr has", ncol(res1)-1, "samples", "\n", append = T)
  dataUnique<-apply(data.matrix(res1[,2:dim(res1)[2]]),2,as.numeric)
  rownames(dataUnique)<-res1[,1]
  
  # Median centering data
  rowMed <- apply(dataUnique,1,median, na.rm=TRUE);
  dataUnique<-dataUnique-rowMed
  
  # Creating storage for results and setting the working directory
  dir.create("PanNETassignerResults", showWarnings = FALSE)
  dir.create("PanNETassignerResults/mRNA", showWarnings = FALSE)
  dir<-getwd()   
  setwd(paste0(dir,"/PanNETassignerResults/mRNA"))
  
  # NMF analysis
  InputFileName_sd_med<-"mrna-1_tumor_and_islets_only_removed_replicates_islet_sd0.8_data_cdt.txt"
  cat("Performing NMF analysis...")
  nmfconsensus(input.ds=InputFileName_sd_med,k.init=2,k.final=10,num.clusterings=20,maxniter=500, error.function="euclidean",data=dataUnique)
  
  # Ordering samples by NMF consensus k=5 cluster
  con<-read.delim("consensus.k.5.gct", stringsAsFactors=FALSE)
  m<-match(con[,1],colnames(dataUnique))
  
  InputFileName_sd_med_NMF<-paste((strsplit(InputFileName_sd_med,".txt")[1]),"_NMF_ordered_5.txt",sep="")
  data_2<-cbind(rownames(dataUnique),dataUnique[,m])
  write.table(data_2,InputFileName_sd_med_NMF, row.names=FALSE,sep="\t",quote=FALSE)
  
  data_2_m<-matrix(as.numeric(data_2[,2:dim(data_2)[2]]),nrow=dim(data_2)[1],ncol=dim(data_2)[2]-1)
  rownames(data_2_m)<-data_2[,1]
  
  # SAM analysis
  cat("Performing SAM analysis...")
  sam.out<-sam(data_2_m,con[,2],rand=123,gene.names=rownames(data_2_m))
  write.table(sam.out@mat.fdr,"samoutput.txt",row.names=FALSE,sep="\t",quote=FALSE)
  
  # Selecting SAM genes with delta value "11.1", which provides 221 genes
  g<-summary(sam.out,11.1)@row.sig.genes
  data_2_sam<-data_2[g,]
  
  OutputFileName_sd_med_NMF_SAM<-paste((strsplit(InputFileName_sd_med_NMF,".txt")[1]),"_SAM_221_genes.txt",sep="")
  write.table(data_2_sam,OutputFileName_sd_med_NMF_SAM,row.names=FALSE,sep="\t",quote=FALSE)
  setwd(dir)
}
