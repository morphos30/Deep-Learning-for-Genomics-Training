

# conda activate /home/mourad/.local/share/r-miniconda/envs/r-reticulate
# install_tensorflow(method = 'conda', envname = '/home/raphael/.local/share/r-miniconda/envs/r-reticulate')


#### LOAD LIBRARIES

library(GenomicRanges)
library(Biostrings)
library(pROC)
library(tensorflow)
library(keras)
library(reticulate)


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
expe="H3K4me3"
DNAletters=c("A","C","G","T")
vocab_size=length(DNAletters)


#### INPUT FILE NAMES

fileBedPos=paste0("data/bed/",expe,"_GM12878_hg19_",kpeaks,"_pos.bed")
fileBedNeg=paste0("data/bed/",expe,"_GM12878_hg19_",kpeaks,"_neg.bed")
fileFastaPos=paste0("data/fasta/",expe,"_GM12878_hg19_",kpeaks,"_pos.fa")
fileFastaNeg=paste0("data/fasta/",expe,"_GM12878_hg19_",kpeaks,"_neg.fa")


#### LOAD PROCESSED DATA

# Load ChIP-seq and control fasta sequences
peakPos.seq=readDNAStringSet(fileFastaPos)
peakNeg.seq=readDNAStringSet(fileFastaNeg)

# Bind positive and negative sequences, and make label
peakAll.seq=c(peakPos.seq,peakNeg.seq)
label=c(rep(1,length(peakPos.seq)),rep(0,length(peakNeg.seq)))

# Shuffle sequence indices
idxS=sample(1:length(peakAll.seq))
peakAllS.seq=peakAll.seq[idxS]
labelS=label[idxS]

# Split train and test indices
percTrain=0.7
idxTrain=1:(ceiling(length(labelS)*percTrain))
idxTest=(length(idxTrain)+1):length(labelS)
labelTrain=labelS[idxTrain]
labelTest=labelS[idxTest]


#### DEEP LEARNING METHODS


# One hot encoding
oneHotTrain=convertOneHot(peakAllS.seq[idxTrain,])
oneHotTest=convertOneHot(peakAllS.seq[idxTest,])

# CNN 1 (simple convolution + 10 dense). 
# Simple yet very efficient model.
kernelSize=16
model <- keras_model_sequential()
model %>% 
  layer_conv_1d(filters = 128, kernel_size = kernelSize, activation = 'relu',
                input_shape = c(peakSize,vocab_size), name="conv1d_cnn1") %>% 
  layer_global_max_pooling_1d() %>% 
  layer_dropout(rate = 0.5) %>% 
  layer_dense(units = 10, activation = "relu", name="10dense_cnn1") %>%
  layer_dense(units = 1, activation = 'linear', name="2ndlastlayer_cnn1") %>%
  layer_dense(units = 1, activation = 'sigmoid', name="lastlayer_cnn1") 

# CNN 2 (parallel convolution). 
# Useful to model different DNA motif sizes.
if(F){
inputs <- layer_input(shape = c(peakSize,vocab_size)) 
conv1 <- inputs %>%
  layer_conv_1d(filters = 256, kernel_size = 8, activation = 'relu') %>% 
  layer_global_max_pooling_1d() %>% 
  layer_dropout(rate = 0.5) 
conv2 <- inputs %>%
  layer_conv_1d(filters = 256, kernel_size = 16, activation = 'relu') %>% 
  layer_global_max_pooling_1d() %>% 
  layer_dropout(rate = 0.5) 
conv3 <- inputs %>%
  layer_conv_1d(filters = 256, kernel_size = 24, activation = 'relu') %>% 
  layer_global_max_pooling_1d() %>% 
  layer_dropout(rate = 0.5) 
outputs <- layer_concatenate(c(conv1,conv2,conv3)) %>%
  layer_dense(units = 10, activation = "relu") %>%
  layer_dense(units = 1, activation = 'sigmoid') 
model <- keras_model(inputs = inputs, outputs = outputs)
}

# CNN 3 (double convolution + 10 dense). 
# Useful for computer vision, but not for DNA sequences.
if(F){
model <- keras_model_sequential()
model %>% 
  layer_conv_1d(filters = 256, kernel_size = 8, activation = 'relu',
                input_shape = c(peakSize,vocab_size)) %>% 
  layer_conv_1d(filters = 256, kernel_size = 3, activation = 'relu',
                input_shape = c(peakSize,vocab_size)) %>% 
  layer_global_max_pooling_1d() %>% 
  layer_dropout(rate = 0.5) %>% 
  layer_dense(units = 10, activation = "relu") %>%
  layer_dense(units = 1, activation = 'sigmoid') 
}

# RNN with LSTM (no CNN). Bad predictions.
# Bad results. We need to include a convolution layer before the lstm layer!
# Slow (not made for efficient parallel computing)
if(F){
model <- keras_model_sequential()
model %>% 
  layer_lstm(128,input_shape = c(peakSize,vocab_size)) %>%
  layer_dropout(rate = 0.5) %>% 
  layer_dense(units = 10, activation = "relu") %>%
  layer_dense(1) %>%
  layer_activation("sigmoid")
}

# CNN + LTSM. 
# Simple yet very efficient model.
if(F){
model <- keras_model_sequential()
model %>% 
  layer_conv_1d(filters = 128, kernel_size = 16, activation = 'relu',
                input_shape = c(peakSize,vocab_size)) %>% 
  layer_dropout(rate = 0.5) %>% 
  layer_lstm(12) %>% 
  layer_dropout(rate = 0.5) %>% 
  layer_dense(units = 20, activation = "relu") %>%
  layer_dense(units = 1, activation = 'sigmoid') 
}

# Model summary
model %>% summary()

# Compile the model
# SGD and Adagrad = bad results
# Adam and RMSprop = good results
model %>% compile(
  optimizer = 'rmsprop', # adam or rmsprop are great here
  loss = 'binary_crossentropy',
  metrics = list('accuracy')
)

# Train the model
history <- model %>% fit(
  oneHotTrain,
  labelTrain,
  epochs = 10, 
  batch_size = 128,
  validation_split = 0.2,
  verbose=1
)
plot(history)

# Prediction
predCNN=predict(model,oneHotTest)
rocCNN_motif=pROC::roc(as.factor(labelTest),predCNN,ci=T)
plot(rocCNN_motif,main=paste0("AUROC: ", round(pROC::auc(rocCNN_motif),3)))

# Accuracy and loss
acc_loss <- model %>% evaluate(oneHotTest, labelTest)
acc_loss

# Save model
file_model=paste0("results/model/CNN1_model_",expe,".hdf5")
save_model_hdf5(model,file_model)





