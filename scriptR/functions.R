


#### LOAD LIBRARIES

library(BSgenome.Hsapiens.UCSC.hg19)
library(BSgenome.Hsapiens.UCSC.hg19.masked)
library(GenomicRanges)
library(gkmSVM)
library(Biostrings)
library(JASPAR2020)
library(TFBSTools)
library(motifmatchr)
library(kebabs)
library(pROC)
library(glmnet)
library(doParallel)
library(ranger)
library(e1071)
library(tensorflow)
library(keras)
library(reticulate)
library(ggseqlogo)
library(motifStack)
library(data.table)


# Function to convert to one-hot encoding
convertOneHot<-function(peak.seq){
 peak.char=strsplit(as.character(peak.seq),'')
 peak.mat=do.call(rbind,peak.char)

 DNAletters=c("A","C","G","T")
 listMat=list()
 for(i in 1:length(DNAletters)){
  mat=matrix(0,nrow(peak.mat),ncol(peak.mat))
  mat[peak.mat==DNAletters[i]]=1
  listMat[[i]]=mat
 }
 arrayout=array(unlist(listMat), dim = c(nrow(listMat[[1]]), ncol(listMat[[1]]), length(listMat)))
 return(arrayout)
}

# Function to convert kernel filters to Position Frequency Matrices (PFMs)
convertToPWMs<-function(conv1d_output,peakPos.seq,activationThreshold){

 conv1d_posmax=apply(conv1d_output,c(1,3),function(x){if(max(x)>activationThreshold){which.max(x)}else{NA}})
 print(dim(conv1d_posmax))
 kernelPFMList=list()
 for(i in 1:ncol(conv1d_posmax)){
  if(sum(!is.na(conv1d_posmax[,i]))>0){
   idxSeqActivatedi=which(!is.na(conv1d_posmax[,i]))
   seqActivatedi.seq=peakPos.seq[idxSeqActivatedi]
   seqActivatedi.char=as.character(seqActivatedi.seq)
   idxStarti=conv1d_posmax[idxSeqActivatedi,i]
   idxEndi=conv1d_posmax[idxSeqActivatedi,i]+kernelSize-1
   subSeqActivatedi.char=as.character(sapply(1:length(seqActivatedi.char),function(x){substring(seqActivatedi.char[x],idxStarti[x],idxEndi[x])}))
   subSeqActivatedi.mat=do.call(rbind,strsplit(subSeqActivatedi.char,""))
   PFMmati=apply(subSeqActivatedi.mat,2,function(x){table(factor(x,levels=c("A","C","G","T")))})
   kernelPFMList[[paste0("kernel",i)]]=PFMatrix(ID=paste0("kernel",i),profileMatrix=PFMmati)
   #print(i)
  }
 }
 kernelPFMList=do.call(PFMatrixList,kernelPFMList)
 return(kernelPFMList)
}

# Function to trim and filter out motifs based on information content
motifTrimming<-function(PFML,ICthreshold){
 PFMTrimmedList=list()
 for(i in 1:length(kernelPFMList)){
  ICmotifi=colSums(toICM(kernelPFMList[[i]]))
  if(sum(ICmotifi>ICthreshold)>0){
   idxLeft = min(which(ICmotifi>ICthreshold))
   idxRight = max(which(ICmotifi>ICthreshold))
   if((idxRight-idxLeft)>5){
    PFMTrimmedList[[ID(kernelPFMList[[i]])]]=PFMatrix(ID = ID(kernelPFMList[[i]]), profileMatrix=as.matrix(kernelPFMList[[i]])[,idxLeft:idxRight])
   }
  }
 }
 PFMTrimmedList=do.call(PFMatrixList,PFMTrimmedList)
 return(PFMTrimmedList)
}




