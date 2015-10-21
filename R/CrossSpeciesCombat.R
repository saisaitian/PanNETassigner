CrossSpeciesCombat <-function(human, mouse, ordered )
{
  #SET YOUR CHOICE OF DIRECTORY
  # human<- "data/mrna-1_tumor_and_islets_only_removed_replicates_islet_normal_like.txt"
  #mouse<-"022910-Rip-Tag-stages-data-log2-genesonly_IT_IC2_only.txt"
  #sdvalue
  sdvalue=0.8
  # Creating storage for results and setting the working directory
  dir.create("PanNETassignerResults", showWarnings = FALSE)
  dir.create("PanNETassignerResults/CrossSpeciesCombat", showWarnings = FALSE)
  dir<-getwd()   
  setwd(paste0(dir,"/PanNETassignerResults/CrossSpeciesCombat"))
  
  # Selecting non-redudant genes
  res1<-screenExpr(human, sdvalue)
  res2<-screenExpr(mouse, sdvalue)
  
  # Reading mouse data
  ta<-res2
  dim(ta)
  
  # Reading human data
  hum<-res1
  dim(hum)
  
  # matching human and mouse data
  ta[,1] <- sub(" ", "", ta[,1])  ## removing spaces between quotes, do this twice.
  ta[,1] <- sub(" ", "", ta[,1]) 
  m<-match(toupper(ta[,1]),hum[,1])
  w<-which(!is.na(m))
  length(w)
  ta_1<-ta[w,]
  hum_1<-hum[m[w],]
  ta_hum<-cbind(hum_1,ta_1)  #[,2:dim(ta_1)[2]])
  dim(ta_hum)
  write.table(ta_hum,"20150105-human_NO_normal_like_mouse_IT_MLP_sd0.8.txt",row.names=FALSE,sep="\t",quote=FALSE)
  
  # Combat analysis
  ta_hum<-cbind(hum_1,ta_1[,2:dim(ta_1)[2]])
  ta_h<-apply(ta_hum[,2:dim(ta_hum)[2]],2,as.numeric)
  ta_hum_m<-data.matrix(ta_h)
  dim(ta_hum)
  bat<-rep(c(1,2),c((dim(hum)[2]-1),(dim(ta_1)[2]-1)))
  combat<-ComBat(dat=ta_hum_m,batch=bat,mod=NULL, par.prior=TRUE, prior.plots=FALSE)
  rownames(combat)<-ta_hum[,1]
  OutFile<-"20150105-human_NO_normal_like_mouse_IT_MLP_sd0.8_ComBat.txt"
  write.table(combat,OutFile, sep="\t",quote=FALSE)
  
  # providing information for Eorder according to subtype information for hierarchical clustering
  dim(ordered)
  com<-read.delim(OutFile, stringsAsFactors=FALSE)
  
  dim(com)
  
  mc<-match(colnames(ordered),colnames(com))
  
  wmc<-which(!is.na(mc))
  
  com_1<-cbind(com[,1],com[,mc[wmc]])
  
  dim(com_1)
  
  com_1_m<-data.matrix(com_1[,2:ncol(com_1)])
  
  med<-apply(com_1_m,1,median)
  
  com_2<-cbind(as.character(com_1[,1]),com_1_m-med)
  
  colnames(com_2)[1]<-"G"
  
  com_3<-rbind(ordered,com_2)
  
  write.table(com_3,"20150105-human_NO_normal_like_mouse_IT_MLP_sd0.8_ComBat_ordered_Eorder.txt", row.names=FALSE,sep="\t",quote=FALSE)
  setwd(dir)
}
