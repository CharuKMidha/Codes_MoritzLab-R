library(ggplot2)
#FILE 1

datafile = read.csv("K562inhouseGoodMergerall_SN_Protein_matrix.csv", header = TRUE, sep =",")

datafile
dim(datafile)
#REOMVE NA containing lines

data_na=na.omit(datafile)
dim(data_na)


#Rmean = paste0("Rmean_",counter)
newdata = data.frame (data_na = data_na[,1], Rmean= rowMeans(data_na[,2:6]))

write.csv(datafile, file ="a.txt") 
write.csv(data_na, file ="a-na.txt") 
write.csv(newdata, file ="a-na_mean.txt") 

#FILE 2

datafile2 = read.csv("K562inhouseRawMergerall_SN_Protein_matrix.csv", header = TRUE, sep =",")

datafile2
dim(datafile2)
#REOMVE NA containing lines

data_na2=na.omit(datafile2)
dim(data_na2)


#Rmean = paste0("Rmean_",counter)
newdata2 = data.frame (data_na2 = data_na2[,1], Rmean2= rowMeans(data_na2[,2:6]))

write.csv(datafile2, file ="a2_raw.txt") 
write.csv(data_na2, file ="a-na2_raw.txt") 
write.csv(newdata2, file ="a-na_mean2_raw.txt") 

##MERGE


data = merge(newdata, newdata2, by.x = "data_na", by.y="data_na2")

write.csv(data, file = "mergefiles_raw.txt")

#####Ploting 

library(ggplot2)
library(scales)
# Open a pdf file

#pdf("K562inhouserawMergerall_SN_Protein_matrix.pdf", width= 8, height = 8)
tiff("K562inhouserawMergerall_SN_Protein_matrixfinal.tif", units = "in", width = 5,height = 5, res = 300)
data = read.csv("mergefiles_raw.txt", header = TRUE, sep =",")
dim(data)
N=NROW(data)
# Default way to draw Scatter Plot
#ggplot(data, aes(x = X.2RT, y = X.3RT)) + geom_point(color = "midnightblue") + geom_smooth()
correlation = cor(data$Rmean, data$Rmean2, method = "pearson")
correlation = sprintf("%1f",correlation)
title = paste("n=", N, "   R2= ",correlation)
title
ggplot(data, aes(x = Rmean , y = Rmean2)) + 
  geom_point(shape = 1, size =1) + 
  geom_smooth( method = "lm", color = 'black',alpha=0, size = .3  )+
  
  #scale_color_gradientn(colors = c("#99FF33","#99ff99","#ffff00","#ff9933", "#ff0000", "#cc3399","#990099","#3399ff","#0000ff","#0099ff","#0099cc"))+
  theme_bw()+
  theme(
    plot.background = element_blank()
    ,panel.grid.major = element_blank()
    ,panel.grid.minor = element_blank()
    #,panel.border = element_blank()
    ,legend.position = "bottom"
  )+
  theme(axis.line = element_line(color = 'black', size = .3))+
  labs(title =title, x = "Good", y = "Raw")+
  scale_x_continuous(labels = scientific)+#, limits = c(0,5000))+
  scale_y_continuous(labels= scientific)# limits = c(0,5000))
dev.off()
