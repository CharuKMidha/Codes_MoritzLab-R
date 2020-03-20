library(ggplot2)
#FILE 1

datafile15 = read.csv("Ecoli_15min_gradient_Ecoli_consensus_T6_100vw_peakview_Protein_Report.xls", header = TRUE, sep ="\t")

datafile15
dim(datafile15)
#REMOVE NA containing lines

data_na15=na.omit(datafile15)
dim(data_na15)

#Rmean = paste0("Rmean_",counter)
head (datafile15)
head(data_na15)

#FILE 2

datafile30 = read.csv("Ecoli_30min_gradient_Ecoli_consensus_T6_100vw_peakview_Protein_Report.xls", header = TRUE, sep ="\t")

datafile30
dim(datafile30)
#REMOVE NA containing lines

data_na30=na.omit(datafile30)

dim(data_na30)


#Rmean = paste0("Rmean_",counter)
head (datafile30)
head(data_na30)

#FILE 3
datafile60 = read.csv("Ecoli_60min_gradient_Ecoli_consensus_T6_100vw_peakview_Protein_Report.xls", header = TRUE, sep ="\t")

datafile60
dim(datafile60)
#REMOVE NA containing lines

data_na60=na.omit(datafile60)
dim(data_na60)


#Rmean = paste0("Rmean_",counter)
head (datafile60)
head(data_na60)

#FILE 4
datafile90 = read.csv("Ecoli_90min_gradient_Ecoli_consensus_T6_100vw_peakview_Protein_Report.xls", header = TRUE, sep ="\t")


datafile90
dim(datafile90)
#REMOVE NA containing lines

data_na90=na.omit(datafile90)
dim(data_na90)


#Rmean = paste0("Rmean_",counter)
head (datafile90)
head(data_na90)

Rmean15= rowMeans(data_na15[,2:6]) #6
Rmean30= rowMeans(data_na30[,2:6])
Rmean60= rowMeans(data_na60[,2:6])
Rmean90= rowMeans(data_na90[,2:6])

Sd15= apply(data_na15[,2:6], 1, sd)
Sd30= apply(data_na30[,2:6], 1, sd)
Sd60= apply(data_na60[,2:6], 1, sd)
Sd90= apply(data_na90[,2:6], 1, sd)

head (Sd15)

CV15 = Sd15*100/Rmean15
CV30 = Sd30*100/Rmean30
CV60 = Sd60*100/Rmean60
CV90 = Sd90*100/Rmean90



median(CV15)
median(CV30)
median(CV60)
median(CV90)

boxplot(CV15, CV30, CV60, CV90, names= c(15, 30, 60, 90), ylim = c(0, 45), border=c("blue","deepskyblue","powderblue", "lightseagreen"), cex = 2, lwd = 3)
