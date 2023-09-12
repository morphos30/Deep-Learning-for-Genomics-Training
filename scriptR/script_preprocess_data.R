


### SETUP DIRECTORY FOR R LIBRARIES
# For teaching at ENSAT
#.libPaths("/usr/local/lib/R/x86_64-pc-linux-gnu-library/4.2")


#### LOAD LIBRARIES

library(BSgenome.Hsapiens.UCSC.hg19)
library(BSgenome.Hsapiens.UCSC.hg19.masked)
library(GenomicRanges)
library(gkmSVM)
library(Biostrings)


#### SETUP PROJECT FOLDER

#setwd("/media/mourad/diskSave/MCF_Toulouse/enseignement/atelier_INSERM/")
setwd("/media/raphael/SSD2/atelier_INSERM")


#### LOAD FUNCTIONS

source("scriptR/functions.R")


#### CREATE FOLDERS

dir.create("data")
dir.create("data/bed")
dir.create("data/fasta")
dir.create("data")
dir.create("results")
dir.create("results/model")
dir.create("results/motif")
dir.create("results/SNP")


#### SOME PARAMETERS

peakSize=201
kpeaks=4000
expe="CTCF"
DNAletters=c("A","C","G","T")


#### INPUT FILE NAMES

fileBedPos=paste0("data/bed/",expe,"_GM12878_hg19_",kpeaks,"_pos.bed")
fileBedNeg=paste0("data/bed/",expe,"_GM12878_hg19_",kpeaks,"_neg.bed")
fileFastaPos=paste0("data/fasta/",expe,"_GM12878_hg19_",kpeaks,"_pos.fa")
fileFastaNeg=paste0("data/fasta/",expe,"_GM12878_hg19_",kpeaks,"_neg.fa")



#### LOAD, PROCESS AND SAVE PROCESSED DATA

# Chromosome information for hg19
Genome=BSgenome.Hsapiens.UCSC.hg19
Chr.V=c(paste0("chr",1:22),"chrX")
SeqInfo=seqinfo(BSgenome.Hsapiens.UCSC.hg19)[Chr.V]

# Load ChIP-seq/DNase-seq peaks
file_peaks=paste0("data/narrowPeak/",expe,"_GM12878_hg19.narrowPeak")
dataPeaks=read.table(file_peaks,sep="\t",header=F)

# Resize peaks
peaks.GR=GRanges(dataPeaks[,1],IRanges(dataPeaks[,2],dataPeaks[,3]),score=dataPeaks[,7])
peaks.GR=resize(peaks.GR, width=peakSize, fix="center")

# Select k best peaks and export to bed
peaksSorted.GR=peaks.GR[order(peaks.GR$score,decreasing=T)]
dataPeaksSorted=as.data.frame(peaksSorted.GR)
dataPeaksBest=dataPeaksSorted[1:kpeaks,]
write.table(dataPeaksBest[,1:3],fileBedPos,sep='\t',col.names=F,row.names=F,quote=F)

# Generate fasta sequences from the experimental peaks
# as well as control fasta sequences
# GC content and repeat distribution similar to the experimental peaks.
# The number of control peaks drawn should be similar to experimental peaks (balanced dataset).
genNullSeqs(inputBedFN=fileBedPos,nMaxTrials=10,xfold=1.5,genomeVersion="hg19",
	outputPosFastaFN=fileFastaPos,outputBedFN=fileBedNeg,outputNegFastaFN=fileFastaNeg,length_match_tol=0)

# genNullSeqs generates sequences of length 200b instead of 201b
# We want to correct this.
dataPeaksNeg=read.table(fileBedNeg,sep="\t",header=F)
peaksNeg.GR=GRanges(dataPeaksNeg[,1],IRanges(dataPeaksNeg[,2],dataPeaksNeg[,3]+1))
peaksNeg.GR=peaksNeg.GR[nchar(as.character(seqnames(peaksNeg.GR)))<6]
write.table(as.data.frame(peaksNeg.GR)[,1:3],fileBedNeg,sep='\t',col.names=F,row.names=F,quote=F)
peaksNeg.seq=getSeq(BSgenome.Hsapiens.UCSC.hg19, names=seqnames(peaksNeg.GR),
     start=start(peaksNeg.GR), end=end(peaksNeg.GR))
names(peaksNeg.seq)=paste0(seqnames(peaksNeg.GR),"_",1:length(peaksNeg.GR))
writeXStringSet(peaksNeg.seq,fileFastaNeg)
peaksPos.GR=GRanges(dataPeaksBest[,1],IRanges(dataPeaksBest[,2],dataPeaksBest[,3]))
peaksPos.seq=getSeq(BSgenome.Hsapiens.UCSC.hg19, names=seqnames(peaksPos.GR),
     start=start(peaksPos.GR), end=end(peaksPos.GR))
names(peaksPos.seq)=paste0(seqnames(peaksPos.GR),"_",1:length(peaksPos.GR))
writeXStringSet(peaksPos.seq,fileFastaPos)





