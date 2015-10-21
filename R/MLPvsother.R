MLPvsother <-function(Indata,Indlabdata)
{
  
  data<-Indata
  cat("Input data has", ncol(data)-1, "samples", "\n", append = T)
  dataUnique<-data.matrix(data[,2:dim(data)[2]])
  rownames(dataUnique)<-data[,1]
  dim(dataUnique)
  
  # Creating storage for results and setting the working directory
  dir.create("PanNETassignerResults", showWarnings = FALSE)
  dir.create("PanNETassignerResults/MLPvsother", showWarnings = FALSE)
  dir<-getwd()   
  setwd(paste0(dir,"/PanNETassignerResults/MLPvsother"))
  
  # Median centering data
  rowMed <- apply(dataUnique,1,median, na.rm=TRUE);
  dataUnique<-dataUnique-rowMed
  dataUnique_1<-cbind(rownames(dataUnique),dataUnique)
  
  # NMF subtypes
  #Mapping MLP vs other subtypes
  con<-Indlabdata
  c1<-which(con[,2] == 1)
  c3<-which(con[,2] == 3)
  c2<-which(con[,2] == 2)
  c4<-which(con[,2] == 4)
  c13<-c(c1,c3)
  c24<-c(c2,c4)
  con_1<-con[c(c13,c24),]
  
  # matching the samples between NMF subtypes and the data file
  m1<-match(con_1[,1],colnames(dataUnique_1))
  data_3<-cbind(dataUnique_1[,1],dataUnique_1[,m1])
  
  InputFileName_to<-"mrna-1_tumor_and_islets_only_removed_replicates_islet_sd0.8_data_cdt_NMF_ordered_5_MLPvsothersubtypes.txt"
  InputFileName_to<-gsub(".*/", "", InputFileName_to, ignore.case=TRUE)
  write.table(data_3,InputFileName_to, row.names=FALSE,sep="\t",quote=FALSE)
  
  data_3_m<-matrix(as.numeric(data_3[,2:dim(data_3)[2]]),nrow=dim(data_3)[1],ncol=dim(data_3)[2]-1)
  rownames(data_3_m)<-data_3[,1]
  
  cl<-rep(c(0,1),c(length(c13),length(c24)))
  
  # SAM analysis
  sam.out1<-sam(data_3_m,cl,rand=123,gene.names=rownames(data_3_m))
  write.table(sam.out1@mat.fdr,"samoutput.txt",row.names=FALSE,sep="\t",quote=FALSE)
  
  g1<-summary(sam.out1,1.7)@row.sig.genes
  data_3_sam<-data_3[g1,]
  
  InputFileName_to_SAM<-paste((strsplit(InputFileName_to,".txt")[1]),"_SAM_49_genes.txt",sep="")
  
  write.table(data_3_sam,"mrna-1_tumor_and_islets_only_sd0.8_removed_replicates_data_cdt_NMF_ordered_5__MLPvsothersubtypes_SAM_49_genes.txt",row.names=FALSE,sep="\t",quote=FALSE)
  setwd(dir)
}
