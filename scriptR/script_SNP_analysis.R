

# conda activate /home/mourad/.local/share/r-miniconda/envs/r-reticulate
# install_tensorflow(method = 'conda', envname = '/home/raphael/.local/share/r-miniconda/envs/r-reticulate')


#### LOAD LIBRARIES

library(BSgenome.Hsapiens.UCSC.hg19)
library(GenomicRanges)
library(Biostrings)
#library(tensorflow)
#library(keras)
library(reticulate)
library(data.table)
library(ggplot2)
library(gplots)


#### SETUP PROJECT FOLDER

#setwd("/media/mourad/diskSave/MCF_Toulouse/enseignement/atelier_INSERM/")
setwd("/media/raphael/SSD2/atelier_INSERM")


#### LOAD FUNCTIONS

source("scriptR/functions.R")

use_condaenv('r-reticulate')
tf$constant("Hello Tensorflow")


#### SOME PARAMETERS

peakSize=201
kpeaks=4000
expe="DNase"
DNAletters=c("A","C","G","T")


#### INPUT FILE NAMES

fileBedPos=paste0("data/bed/",expe,"_GM12878_hg19_",kpeaks,"_pos.bed")
fileBedNeg=paste0("data/bed/",expe,"_GM12878_hg19_",kpeaks,"_neg.bed")
fileFastaPos=paste0("data/fasta/",expe,"_GM12878_hg19_",kpeaks,"_pos.fa")
fileFastaNeg=paste0("data/fasta/",expe,"_GM12878_hg19_",kpeaks,"_neg.fa")




#### DEEP LEARNING PREDICTIONS OF SNP EFFECTS

# Load CNN1 model
file_model=paste0("results/model/CNN1_model_",expe,".hdf5")
model=load_model_hdf5(file_model)

# rs907091 and rs9303277 were found in a GWAS of childhood acute lymphoblastic leukemia
# Article: https://www.nature.com/articles/s41467-017-02596-9 (supplementary file)
# rs907091 rs9303277 are also eQTLs in blood (https://www.gtexportal.org/home/snp/rs9303277)
SNP.GR1=GRanges("chr17",IRanges(c(37921742,37976469,38029120)),ref=c("C","C","C"),alt=c("T","T","G"),strand='+',
	rs=c("rs907091","rs9303277","rs12936231"),OR_GWAS=c(1.17,0.86,0.86), Effect_eQTL=c(0.32,-0.33,-0.34))

# rs7329174 was found from GWAS catalog. It's associated with Crohn's disease and located in intron.
SNP.GR2=GRanges("chr13",IRanges(c(41558110)),ref=c("A"),alt=c("G"), strand='+', rs=c(" rs7329174"),OR_GWAS=c(1.27))

# rs142811167 and rs10821936 were found from GWAS catalog. It's associated with leukemia and located in intron.
SNP.GR3=GRanges(c("chr11","chr10"),IRanges(c(7719585,63723577)),ref=c("C","T"),alt=c("T","C"), strand='+', rs=c("rs142811167","rs10821936"),OR_GWAS=c(87,1.8))

SNP.GR=c(SNP.GR1,SNP.GR2,SNP.GR3)

# Extract sequences surrounding the SNPs
genomeVersion=BSgenome.Hsapiens.UCSC.hg19
region.GRr=resize(SNP.GR,peakSize,fix="center")
region.seq=as.character(getSeq(genomeVersion, names=seqnames(region.GRr),
     start=start(region.GRr), end=end(region.GRr)))

# Make sequences with reference alleles
SNPposrel=ceiling((peakSize+1)/2)
region.seqRef=region.seq
substring(region.seqRef,SNPposrel,SNPposrel)=as.character(region.GRr$ref)
names(region.seqRef)=paste0("seq",1:length(region.seqRef))
#writeXStringSet(DNAStringSet(region.seqRef),filepath=paste0("data/fasta/sequences_SNPref_hg19.fa"))

# Make sequences with alternative alleles
region.seqAlt=region.seq
substring(region.seqAlt,SNPposrel,SNPposrel)=as.character(region.GRr$alt)
names(region.seqAlt)=paste0("seq",1:length(region.seqAlt))
#writeXStringSet(DNAStrinregionMaxProbgSet(region.seqAlt),filepath=paste0(pathfasta,"/sequences_SNPalt_hg19.fa"))

# Convert sequences to one hot encoding
oneHotRef=convertOneHot(region.seqRef)
oneHotAlt=convertOneHot(region.seqAlt)

# Predict SNP effect on the predicted probability using the deep learning model
predProbRef=predict(model,oneHotRef)
predProbAlt=predict(model,oneHotAlt)
deltaProbSNP=predProbAlt-predProbRef
print(deltaProbSNP)

# Predict SNP effects on the second to last layer (just before logistic transformation) using the deep learning model
secondlast_layer_model <- keras_model(inputs = model$input, outputs = get_layer(model, '2ndlastlayer_cnn1')$output)
secondlast_outputRef <- predict(secondlast_layer_model, oneHotRef)
secondlast_outputAlt <- predict(secondlast_layer_model, oneHotAlt)
delta2ndLastSNP=secondlast_outputAlt-secondlast_outputRef
print(delta2ndLastSNP)

# Mutation map around the 1st SNP
region.seq1=region.seq[1]
mutationMap=matrix(NA,4,nchar(region.seq1))
valRef=predict(secondlast_layer_model, convertOneHot(region.seq1))
for(i in 1:nchar(region.seq1)){
 for(j in 1:4){
  region.seqMut=region.seq1
  substring(region.seqMut,i,i)=DNAletters[j]
  oneHotMut=convertOneHot(region.seqMut)
  secondlast_outputMut=predict(secondlast_layer_model, oneHotMut)
  mutationMap[j,i]=secondlast_outputMut - valRef
 }
}
rownames(mutationMap)=DNAletters

# Plot mutation maps (heatmaps)
file_mutationMap_large=paste0("results/SNP/mutationMap_",peakSize,"b_",expe,".pdf")
pdf(file_mutationMap_large,30,4)
heatmap.2(mutationMap, trace="none", scale="none", dendrogram="none", Rowv=FALSE, Colv=FALSE, col=colorRampPalette(colors = c("blue","white","red")))
dev.off()

file_mutationMap_large=paste0("results/SNP/mutationMap_",41,"b_",expe,".pdf")
pdf(file_mutationMap_large,10,4)
mapRange=(ceiling(peakSize/2)-20):(ceiling(peakSize/2)+20)
heatmap.2(mutationMap[,mapRange], trace="none", scale="none", dendrogram="none", Rowv=FALSE, Colv=FALSE, col=colorRampPalette(colors = c("blue","white","red")))
dev.off()






