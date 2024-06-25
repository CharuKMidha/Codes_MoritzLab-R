
##PEPTIDE#######################################################################
library(limma) 
library(tidyverse)
library(edgeR) 
library(sva)
library(psych)# we use this for the pairwise correlation plots



#### FILE LOADING 35GB
############NOW OPEN ME NEXT TIME long4col.csv at line 100 or so ##################################
data_start2 <- as.data.frame(read.csv("//cello/proteomics/mhoopman/data/Longevity/Human/phaseI.May2024.pep.quant.txt", sep = '\t'))

###THIS CODE SECTION ESTIMATES THE SUMMARY OF ALL THE COLUMNS WE ARE INTERESTED IN ##
#names(data_start2) <- c("Charge", "Proteins", "Run", "Scan", "Pool", "Fraction",
#                        "TMTChannel", "Abundance", "Probability", "SampleID", "Cohort")
# pepU<- unique(data_start2$Peptide)
# chargeU<- unique(data_start2$Charge)
# protU<- unique(data_start2$Proteins)
# runU<- unique(data_start2$Run)
# scanU<- unique(data_start2$Scan)
#length(scanU)
# poolU<- unique(data_start2$Pool)
# fractionU<- unique(data_start2$Fraction)
# channelU<- unique(data_start2$TMTChannel)
# minp<- min(data_start2$Probability)
# maxp<- max(data_start2$Probability)
# sampleidU<- unique(data_start2$SampleID)
# cohortU<- unique(data_start2$Cohort)

data<- data_start2 #MADE A COPY OF original data
dim(data_start2)
head( data)

data$Fraction<- NULL #REMOVINF EXTRA COLUMNS
data$Probability<- NULL #REMOVINF EXTRA COLUMNS
data$Run<- NULL #REMOVINF EXTRA COLUMNS
data$Charge<- NULL #REMOVINF EXTRA COLUMNS

head(data)

#REMOVING DECOYS, CONTAMS, KEEPING PROTEOTYPIC only
decoy<- data[!grepl("_HUMAN", data$Proteins), ]
data<-  data[grepl("_HUMAN", data$Proteins), ] #removed decoys kept mixed ones
dataproteingp<- data[grepl(";", data$Proteins), ] #proteotypic psms
data<- data[!grepl(";", data$Proteins), ] #removed proteotypic
data[grepl("DECOY", data$Proteins), ] #removed proteotypic should b 0

dim(decoy)
dim(data)
dim(dataproteingp)
dim(data_start2)
head (data)

#JOINING all details 
data$sampledetail<- paste0(data$SampleID,"#",data$Cohort, "@", data$Peptide, "!", data$Proteins)
data$poolscan<- paste0(data$Pool, "^" ,data$Scan)
data$SampleID <- NULL #REMOVINF EXTRA COLUMNS
data$Cohort<- NULL #REMOVINF EXTRA COLUMNS
data$Pool<- NULL #REMOVINF EXTRA COLUMNS
data$Scan<- NULL #REMOVINF EXTRA COLUMNS
data$Peptide<- NULL #REMOVINF EXTRA COLUMNS
data$Proteins<- NULL #REMOVINF EXTRA COLUMNS
head(data)

dim(data)
rownames(data)<- paste0(data$sampledetail, "%", data$poolscan)
data$sampledetail<- NULL #REMOVINF EXTRA COLUMNS
head(data)
write.csv(data, "long4col.csv") #4column file sample details TMTChannel Abundance poolscan

head(data)
# TMTChannel Abundance poolscan
# ID4662#Cent@n[305.21]HGHGHGK[432.30]!sp|P01042|KNG1_HUMAN%1^1664          1   1261.59   1^1664
# ID0002#MrOS@n[305.21]HGHGHGK[432.30]!sp|P01042|KNG1_HUMAN%1^1664          2   2745.59   1^1664
# ID0003#MrOS@n[305.21]HGHGHGK[432.30]!sp|P01042|KNG1_HUMAN%1^1664          3   1600.31   1^1664
# ID0004#MrOS@n[305.21]HGHGHGK[432.30]!sp|P01042|KNG1_HUMAN%1^1664          4   1903.98   1^1664
# ID0005#MrOS@n[305.21]HGHGHGK[432.30]!sp|P01042|KNG1_HUMAN%1^1664          5   3585.63   1^1664
# ID0006#MrOS@n[305.21]HGHGHGK[432.30]!sp|P01042|KNG1_HUMAN%1^1664          6  11713.82   1^1664




##PEPTIDE#######################################################################
library(limma) 
library(tidyverse)
library(edgeR) 
library(sva)
library(psych)# we use this for the pairwise correlation plots
library(reshape2) 
library(reshape) 
data<- read.csv("long4col.csv")
head(data)
sortabundance<- order(data$Abundance, decreasing = TRUE)
head(sortabundance)
rownames(data)<- data$X; data$X<- NULL
length(unique(row.names(data)))
dim(data)
head(data)

#MAKING CHANNEL WISE MATRIX 18 channels column and pool^scan as rows

finaldata = as.data.frame(data) %>%
  spread(key = TMTChannel ,
         value = Abundance)

head(finaldata)
dim(finaldata) #1386708      19

length(rownames(data))
head(rownames(data))
rownames(finaldata)<- finaldata$poolscan
finaldata$poolscan<- NULL
head(finaldata)
poolscan_channelmatrix<- finaldata #keeping finaldata safe it has entries with 0value in channel 18
#write.csv(finaldata, "channelfile.csv") #18 channels column and pool^scan as rownames

###########TESTING#######################
#3pools for testing

#sample(1:201,15, replace=F)
#[1] 177  87 187  47  20  56 182 115 141 130  31 144  22  30  33
# pool1<- finaldata[grep("^1\\^", rownames(finaldata)), ]
# pool120<- finaldata[grep("^120\\^", rownames(finaldata)), ]
# pool171<- finaldata[grep("^171\\^", rownames(finaldata)), ]
# pool177<- finaldata[grep("^177\\^", rownames(finaldata)), ]
# pool87<- finaldata[grep("^87\\^", rownames(finaldata)), ]
# 
# pool187<- finaldata[grep("^187\\^", rownames(finaldata)), ]
# pool47<- finaldata[grep("^47\\^", rownames(finaldata)), ]
# pool20<- finaldata[grep("^20\\^", rownames(finaldata)), ]
# pool56<- finaldata[grep("^56\\^", rownames(finaldata)), ]
# pool182<- finaldata[grep("^182\\^", rownames(finaldata)), ]
# 
# 
# pool115<- finaldata[grep("^115\\^", rownames(finaldata)), ]
# pool130<- finaldata[grep("^130\\^", rownames(finaldata)), ]
# pool31<- finaldata[grep("^131\\^", rownames(finaldata)), ]
# pool22<- finaldata[grep("^22\\^", rownames(finaldata)), ]
# pool144<- finaldata[grep("^144\\^", rownames(finaldata)), ]
# 
# threepools<- rbind.data.frame(pool1, pool120, pool171, pool177, pool87,
#                               pool187, pool47, pool20, pool56, pool182,
#                               pool115, pool130, pool31, pool22, pool144)
threepools<- finaldata
head(threepools)
tail(threepools)
threepools$sampledetail<- rownames(threepools)
rownames(threepools)<- NULL
threepools_long<- melt(threepools)
colnames(threepools_long)<- c("sampledetail", "TMTChannel", "Raw_abundance")
threepools_long$TMTChannel<- paste0(threepools_long$sampledetail,"_", threepools_long$TMTChannel,"!")
head(threepools_long)
sortabundance<- sort(threepools_long$Raw_abundance,decreasing = TRUE)
tail(sortabundance)
head(sortabundance)
threepools_long$sampledetail<- NULL

data<- datacp
head(data)
data$sampledetail<- row.names(data)
row.names(data)<- NULL

data$TMTChannel<- paste0(data$poolscan,"_",data$TMTChannel,"!")
head(data)
head(threepools_long)


data_raw_longmerged <- left_join(threepools_long, 
                                data, by = "TMTChannel")
head(data_raw_longmerged)
data_raw_longmerged$Abundance<- NULL
data_raw_longmerged$poolscan<- NULL
data_raw_longmerged$sampledetail<- gsub("%.*","", data_raw_longmerged$sampledetail)
data_raw_longmerged$sampledetail<- gsub("@.*!", "@",data_raw_longmerged$sampledetail)
head(data_raw_longmerged)
length(unique(data_raw_longmerged$sampledetail))
length(unique(data_raw_longmerged$TMTChannel))

#DISSOLVED SCAN NUMBER
data_raw_longmerged$TMTChannel<- gsub("\\^.*_","^",data_raw_longmerged$TMTChannel)
data_raw_longmerged$TMTChannel<- gsub("!","",data_raw_longmerged$TMTChannel)


data_raw_longmerged$sampledetail<- 
  paste0(data_raw_longmerged$TMTChannel,"%",
         data_raw_longmerged$sampledetail)
sortabundance<- sort(data_raw_longmerged$Raw_abundance,decreasing = TRUE)
tail(sortabundance)
head(sortabundance)

tail(data_raw_longmerged)

data_raw_longmerged$TMTChannel<- NULL

data_raw_longprotein_sum<- aggregate(x = data_raw_longmerged[c("Raw_abundance")], 
                                          by = data_raw_longmerged[c("sampledetail")], 
                                          FUN = sum)
dim(data_raw_longmerged)
dim(data_raw_longprotein_sum)
head(data_raw_longprotein_sum)
data_raw_longprotein_sum$protein<- gsub(".*@","",data_raw_longprotein_sum$sampledetail)
data_raw_longprotein_sum$sampledetail<- gsub("@.*", "",data_raw_longprotein_sum$sampledetail)

sortabundance<- sort(data_raw_longprotein_sum$Raw_abundance,decreasing = TRUE)
tail(sortabundance)
head(sortabundance)
head(data_raw_longprotein_sum)

#MAKING PROTEIN MATRIX
protein_raw_matrix = as.data.frame(data_raw_longprotein_sum) %>%
  spread(key = sampledetail ,
         value = Raw_abundance)

head(protein_raw_matrix)
dim(protein_raw_matrix) #1386708      19

rownames(protein_raw_matrix)<- protein_raw_matrix$protein
protein_raw_matrix$protein<- NULL
head(protein_raw_matrix)
protein_raw_matrix<- na.omit(protein_raw_matrix)
dim(protein_raw_matrix)
#protein_raw_matrix[sapply(protein_raw_matrix, is.infinite)] <- NA
write.csv(protein_raw_matrix,"protein_raw_matrixcpproteins.csv")

poolcols<- paste0("!",colnames(protein_raw_matrix))
poolno<- unique(gsub("\\^.*", "", colnames(protein_raw_matrix)))

pdf("testingfulldataphase1completedata.pdf", width = 30, height = 15)

boxplot(log(protein_raw_matrix,10), main= "Raw")
boxplot(log(protein_raw_matrix[,grepl("QC", colnames(protein_raw_matrix))],10), main = "Raw QC")
dim(protein_raw_matrix)
#dev.off()####1) SL
subset_protein<- list()
subset_sl_protein<- list()

# matched_cols <- grep(paste0("!", 1,"^"), poolcols, fixed = TRUE)
# subset_protein[1] <- protein_raw_matrix[,matched_cols ]
# head(subset_protein[1])
# columnsum<- colSums(subset_protein[1], na.rm = TRUE)
# NF <- columnsum/columnsum[18]
# subset_sl_protein[1]<- sweep(subset_protein[1] ,2, NF, "/")
# head(subset_sl_protein[1])
# colSums(subset_sl_protein[1], na.rm = TRUE)#should all be same

target<- mean(colSums(protein_raw_matrix, na.rm = TRUE))

# Find col indices that match the current pattern
#colnames(VSN_normalized_protein_matrix)<- paste0("!",colnames(VSN_normalized_protein_matrix))
for (i in poolno){
  matched_cols <- grep(paste0("!", i,"^"), poolcols, fixed = TRUE)
  subset_protein[[i]] <- protein_raw_matrix[,matched_cols ]
  head(subset_protein[[i]])
  #columnsum<- colSums(subset_protein[[i]], na.rm = TRUE)
  NF <- target/colSums(subset_protein[[i]], na.rm = TRUE)
  subset_sl_protein[[i]]<- sweep(subset_protein[[i]] ,2, NF, "*")
  head(subset_sl_protein[[i]])
  colSums(subset_sl_protein[[i]], na.rm = TRUE)#should all be same
}

SL_normalized_protein_matrix<- do.call(cbind, subset_sl_protein)
#boxplot(log(SL_normalized_protein_matrix,10), main= "SL")
#boxplot(log(SL_normalized_protein_matrix[,grepl("QC", colnames(SL_normalized_protein_matrix))],10), main = "SL QC")

#plot(colSums(subset_sl_protein[[3]], na.rm = TRUE) ,main = "sum of pool3")#should all be same
#plot(colSums(subset_sl_protein[[2]], na.rm = TRUE),main = "sum of pool2")#should all be same
write.csv(SL_normalized_protein_matrix,"SL_normalized_protein_matrix.csv")

####################IRS
####1) SL IRS
SL_normalized_protein_matrix
dim(SL_normalized_protein_matrix)
#nonaSL<- na.omit(SL_normalized_protein_matrix)

irs<- SL_normalized_protein_matrix[,grepl("QC", colnames(SL_normalized_protein_matrix))]
write.csv(irs, "QC_sl.csv")
irs$average<- rowMeans(irs, na.rm = TRUE)

subset_protein<- list()
sl_irs<- list()


# # Find col indices that match the current pattern
# # 
# matched_cols <- grep(paste0("!", 1,"^"), poolcols, fixed = TRUE)
# subset_protein[[1]] <- SL_normalized_protein_matrix[,matched_cols ]
# subset_protein[[1]]
# #scaling factor
# irsNF<- irs$average/ subset_protein[[1]][,grepl("QC", colnames(subset_protein[[1]]))]
# irsNF[sapply(irsNF, is.infinite)] <- NA
# irsNF[sapply(irsNF, is.na)] <- 1
# #irsNF<- as.numeric(unlist(irsNF))
# dim(subset_protein[[1]])
# length(irsNF)
# head(irsNF)
# sl_irs[[1]]<- subset_protein[[1]] * irsNF
# colSums(sl_irs[[1]])
#write.csv(SL_normalized_protein_matrix,"SL_normalized_protein_matrixnew.csv")
for (i in poolno) {
  # Find col indices that match the current pattern
  #colnames(SL_normalized_protein_matrix)<- paste0("!",colnames(SL_normalized_protein_matrix))
  matched_cols <- grep(paste0("!", i,"^"), poolcols, fixed = TRUE)
  subset_protein[[i]] <- SL_normalized_protein_matrix[,matched_cols ]
  subset_protein[[i]]
  #scaling factor
  #irsNF<- irs$average/ rowSums(subset_protein[[i]], na.rm = TRUE)#[grep(".*QC", colnames(subset_protein[[i]]))]/QC_mean
  irsNF<- irs$average/ subset_protein[[i]][,grepl("QC", colnames(subset_protein[[i]]))]
  irsNF[sapply(irsNF, is.infinite)] <- NA
  irsNF[sapply(irsNF, is.na)] <- 1 #making factor 1 for na val as we multiply
  sl_irs[[i]]<- subset_protein[[i]] * irsNF
  
  #irsNF<- as.numeric(unlist(irsNF))
  
  # subset_protein1[is.na(subset_protein1)]<- 1
  
  
  #sl_irs[[i]]<- subset_protein[[i]] / rep(irsNF, ncol(subset_protein[[i]]))
  #SL IRS DONE
}

sl_irs_protein<- do.call(cbind, sl_irs)
sl_irs_protein[sapply(sl_irs_protein, is.infinite)] <- NA
head(sl_irs_protein[1:5,1:5])
#boxplot(log(sl_irs_protein,10), main= "sl_irs_protein")
boxplot(log(SL_normalized_protein_matrix[,grepl("QC", colnames(SL_normalized_protein_matrix))],10), main = "sl_protein QC")
boxplot(log(sl_irs_protein[,grepl("QC", colnames(sl_irs_protein))],10), main = "sl_irs_protein QC")
#boxplot(log(protein_raw_matrix[,grepl("QC", colnames(protein_raw_matrix))],10), main = "raw_protein QC")
#plot(colSums(sl_irs[[3]], na.rm = TRUE), main = "sl_irs colsum pool3")#should all be same
#plot(colSums(sl_irs[[2]], na.rm = TRUE), main = "sl_irs colsum pool2")#should all be same
write.csv(sl_irs_protein, "sl_irs_protein.csv")
#dev.off()

###############################################################



#2) MEDIAN
head(protein_raw_matrix)
dim(protein_raw_matrix)
library(proBatch)
median_normalized_protein_matrixlog2 = normalize_data_dm(log(protein_raw_matrix,2),
                                                        normalize_func = 'medianCentering')

boxplot(median_normalized_protein_matrixlog2, main= "median")
#boxplot(median_normalized_protein_matrixlog2[,grepl("QC", colnames(median_normalized_protein_matrixlog2))], main = "median QC")
#plot(colSums(median_normalized_protein_matrixlog2, na.rm = TRUE), main = "median colsum")#should all be same
write.csv(median_normalized_protein_matrixlog2, "median_normalized_protein_matrixlog2.csv")
########MEDIAN IRS
irs<- median_normalized_protein_matrixlog2[,grepl("QC", colnames(median_normalized_protein_matrixlog2))]
write.csv(irs, "QC_median.csv")
irs$average<- rowMeans(irs, na.rm = TRUE)

subset_protein<- list()
median_irs<- list()

for (i in poolno) {
  # Find col indices that match the current pattern
  #colnames(SL_normalized_protein_matrix)<- paste0("!",colnames(SL_normalized_protein_matrix))
  matched_cols <- grep(paste0("!", i,"^"), poolcols, fixed = TRUE)
  subset_protein[[i]] <- median_normalized_protein_matrixlog2[,matched_cols ]
  subset_protein[[i]]
  #scaling factor
  irsNF<- irs$average/ subset_protein[[i]][,grepl("QC", colnames(subset_protein[[i]]))]
  #irsNF[sapply(irsNF, is.infinite)] <- NA
  #irsNF[sapply(irsNF, is.na)] <- 1 #making factor 1 for na val as we multiply
  median_irs[[i]]<- subset_protein[[i]] * irsNF
  
  #irsNF<- as.numeric(unlist(irsNF))
  
  # subset_protein1[is.na(subset_protein1)]<- 1
  
  
  #sl_irs[[i]]<- subset_protein[[i]] / rep(irsNF, ncol(subset_protein[[i]]))
  #SL IRS DONE
}

median_irs_protein<- do.call(cbind, median_irs)
median_irs_protein[sapply(median_irs_protein, is.infinite)] <- NA

head(median_irs_protein[1:5,1:5])
#boxplot(median_irs_protein, main= "median_irs_protein")
boxplot(median_normalized_protein_matrixlog2[,grepl("QC", colnames(median_normalized_protein_matrixlog2))], main = "median_protein QC")
boxplot(median_irs_protein[,grepl("QC", colnames(median_irs_protein))], main = "median_irs_protein QC")
#boxplot(log(protein_raw_matrix[,grepl("QC", colnames(protein_raw_matrix))],10), main = "raw_protein QC")
#plot(colSums(median_irs_protein, na.rm = TRUE), main = "median_irs colsum")#should all be same
write.csv(median_irs_protein, "median_irs_protein.csv")




#3) QUANTILE
#BiocManager::install("preprocessCore", configure.args="--disable-threading", force = TRUE)
quantile_normalized_protein_matrixlog2 = normalize_data_dm(log(as.matrix(protein_raw_matrix),2),
                                                          normalize_func = 'quantile')
quantile_normalized_protein_matrixlog2[sapply(quantile_normalized_protein_matrixlog2, is.infinite)] <- NA
boxplot(quantile_normalized_protein_matrixlog2, main= "Quantile")
#boxplot(quantile_normalized_protein_matrixlog2[,grepl("QC", colnames(quantile_normalized_protein_matrixlog2))], main = "Quantile QC")
#plot(colSums(quantile_normalized_protein_matrixlog2, na.rm = TRUE), main = "Quantile colsum")#should all be same

write.csv(quantile_normalized_protein_matrixlog2, "quantile_normalized_protein_matrixlog2.csv")

########quantile IRS
irs<- quantile_normalized_protein_matrixlog2[,grepl("QC", colnames(quantile_normalized_protein_matrixlog2))]
write.csv(irs, "QC_quantile.csv")
irs$average<- rowMeans(irs, na.rm = TRUE)

subset_protein<- list()
quantile_irs<- list()

for (i in poolno) {
  # Find col indices that match the current pattern
  #colnames(SL_normalized_protein_matrix)<- paste0("!",colnames(SL_normalized_protein_matrix))
  matched_cols <- grep(paste0("!", i,"^"), poolcols, fixed = TRUE)
  subset_protein[[i]] <- quantile_normalized_protein_matrixlog2[,matched_cols ]
  subset_protein[[i]]
  #scaling factor
  irsNF<- irs$average/ subset_protein[[i]][,grepl("QC", colnames(subset_protein[[i]]))]
  #irsNF[sapply(irsNF, is.infinite)] <- NA
  #irsNF[sapply(irsNF, is.na)] <- 1 #making factor 1 for na val as we multiply
  quantile_irs[[i]]<- subset_protein[[i]] * irsNF
  
  #irsNF<- as.numeric(unlist(irsNF))
  
  # subset_protein1[is.na(subset_protein1)]<- 1
  
  
  #sl_irs[[i]]<- subset_protein[[i]] / rep(irsNF, ncol(subset_protein[[i]]))
  #SL IRS DONE
}

quantile_irs_protein<- do.call(cbind, quantile_irs)
quantile_irs_protein[sapply(quantile_irs_protein, is.infinite)] <- NA
head(quantile_irs_protein[1:5,1:5])
#boxplot(quantile_irs_protein, main= "quantile_irs_protein")
boxplot(quantile_normalized_protein_matrixlog2[,grepl("QC", colnames(quantile_normalized_protein_matrixlog2))], main = "quantile_protein QC")
boxplot(quantile_irs_protein[,grepl("QC", colnames(quantile_irs_protein))], main = "quantile_irs_protein QC")
#boxplot(log(protein_raw_matrix[,grepl("QC", colnames(protein_raw_matrix))],10), main = "raw_protein QC")
#plot(colSums(quantile_irs_protein, na.rm = TRUE), main = "quantile_irs colsum")#should all be same
write.csv(quantile_irs_protein, "quantile_irs_protein.csv")




#4) VSN
VSN_normalized_protein_matrix<- normalizeVSN(protein_raw_matrix)
#boxplot(log(VSN_normalized_protein_matrix,10), main= "VSN")
#boxplot(log(VSN_normalized_protein_matrix[,grepl("QC", colnames(VSN_normalized_protein_matrix))],10), main = "VSN QC")

write.csv(VSN_normalized_protein_matrix, "VSN_normalized_protein_matrix.csv")
########vsn IRS
irs<- VSN_normalized_protein_matrix[,grepl("QC", colnames(VSN_normalized_protein_matrix))]
write.csv(irs, "QC_vsn.csv")
irs$average<- rowMeans(irs, na.rm = TRUE)

subset_protein<- list()
vsn_irs<- list()

for (i in poolno) {
  # Find col indices that match the current pattern
  #colnames(SL_normalized_protein_matrix)<- paste0("!",colnames(SL_normalized_protein_matrix))
  matched_cols <- grep(paste0("!", i,"^"), poolcols, fixed = TRUE)
  subset_protein[[i]] <- VSN_normalized_protein_matrix[,matched_cols ]
  subset_protein[[i]]
  #scaling factor
  irsNF<- irs$average/ subset_protein[[i]][,grepl("QC", colnames(subset_protein[[i]]))]
  #irsNF[sapply(irsNF, is.infinite)] <- NA
  #irsNF[sapply(irsNF, is.na)] <- 1 #making factor 1 for na val as we multiply
  vsn_irs[[i]]<- subset_protein[[i]] * irsNF
  
  #irsNF<- as.numeric(unlist(irsNF))
  
  # subset_protein1[is.na(subset_protein1)]<- 1
  
  
  #sl_irs[[i]]<- subset_protein[[i]] / rep(irsNF, ncol(subset_protein[[i]]))
  #SL IRS DONE
}

vsn_irs_protein<- do.call(cbind, vsn_irs)
vsn_irs_protein[sapply(vsn_irs_protein, is.infinite)] <- NA
head(vsn_irs_protein[1:5,1:5])
boxplot(vsn_irs_protein, main= "vsn_irs_protein")
boxplot(VSN_normalized_protein_matrix[,grepl("QC", colnames(VSN_normalized_protein_matrix))], main = "vsn_protein QC")
boxplot(vsn_irs_protein[,grepl("QC", colnames(vsn_irs_protein))], main = "vsn_irs_protein QC")
#boxplot(log(protein_raw_matrix[,grepl("QC", colnames(protein_raw_matrix))],10), main = "raw_protein QC")
#plot(colSums(vsn_irs_protein, na.rm = TRUE), main = "vsn_irs colsum")#should all be same
write.csv(vsn_irs_protein, "vsn_irs_protein.csv")
dev.off()




