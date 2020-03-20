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

Rmean_log15= sort(log10(rowMeans(data_na15[,2:6])), decreasing = TRUE)
Rmean_log30= sort(log10(rowMeans(data_na30[,2:6])), decreasing = TRUE)
Rmean_log60= sort(log10(rowMeans(data_na60[,2:6])), decreasing = TRUE)
Rmean_log90= sort(log10(rowMeans(data_na90[,2:6])), decreasing = TRUE)

head (Rmean_log15)
plot(Rmean_log15)

write.csv(Rmean_log15, file = "Rmean_log15.txt")
write.csv(Rmean_log30, file = "Rmean_log30.txt")
write.csv(Rmean_log60, file = "Rmean_log60.txt")
write.csv(Rmean_log90, file = "Rmean_log90.txt")

data = read.csv("Rmean_log15_30_60_90.txt", header = TRUE, sep ="\t")


plot(data$X15min, type='p', col= "blue", cex= 2, lwd= 2, xlab="Number of Protein groups", ylab='Protein groups abudance (log10)', ylim = c(1, 7))
points(data$X30min, type='p', cex= 2, lwd= 2,col = "deepskyblue")
points(data$X60min, type='p', cex= 2, lwd= 2,col = "powderblue")
points(data$X90min, type='p', cex= 2, lwd= 2,col = "lightseagreen")  
