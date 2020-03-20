library(ggplot2)
#FILE 1

datafile = read.csv("Ecoli_90minproteotypic.txt", header = TRUE, sep ="\t")

datafile
dim(datafile)
#REMOVE NA containing lines

data_na=na.omit(datafile)
dim(data_na)


#Rmean = paste0("Rmean_",counter)
head (datafile)

newdata = data.frame (data_na = data_na[,1], Rmean= rowMeans(data_na[,2:6]))
head(newdata)
head(data_na)
write.csv(datafile, file ="a.txt") 
write.csv(data_na, file ="a-na.txt") 
write.csv(newdata, file ="a-na_mean.txt") 

#FILE 2

datafile2 = read.csv("Ecoli_60minproteotypic.txt", header = TRUE, sep ="\t")

datafile2
dim(datafile2)
#REMOVE NA containing lines

data_na2=na.omit(datafile2)
dim(data_na2)
head(data_na2)

#Rmean = paste0("Rmean_",counter)
newdata2 = data.frame (data_na2 = data_na2[,1], Rmean2= rowMeans(data_na2[,2:6]))
head (newdata2)
write.csv(datafile2, file ="a2_60min.txt") 
write.csv(data_na2, file ="a-na2_60min.txt") 
write.csv(newdata2, file ="a-na_mean2_60min.txt") 

##MERGE

data = merge(newdata, newdata2, by.x = "data_na", by.y="data_na2")

head (data)
write.csv(data, file = "mergefiles_60min.txt")

#####Ploting 

library(ggplot2)
library(scales)
# Open a pdf file

#pdf("K562inhouse60minMergerall_SN_Protein_matrix.pdf", width= 8, height = 8)
#tiff("90min-60minPHL_SN_Peptidesfinal.tif", units = "in", width = 5,height = 5, res = 300)
data = read.csv("mergefiles_60min.txt", header = TRUE, sep =",")
dim(data)
N=NROW(data)
# Default way to d60min Scatter Plot
#ggplot(data, aes(x = X.2RT, y = X.3RT)) + geom_point(color = "midnightblue") + geom_smooth()
correlation = cor(data$Rmean, data$Rmean2, method = "pearson")
correlation = sprintf("%1f",correlation)
title = paste("n=", N, "   R2= ",correlation)
title
ggplot(data, aes(x = Rmean , y = Rmean2)) + 
  geom_point(shape = 1, size =4, color = "#008080") + #2,8,red earlier values
  #geom_point(shape="\u25D6", colour="red", size=4) +
  #geom_point(size =4, color = "red") + 
  geom_smooth( method = "lm", color = 'black',alpha=0, size = 1  )+
  
  #scale_color_gradientn(colors = c("#99FF33","#99ff99","#ffff00","#ff9933", "#ff0000", "#cc3399","#990099","#3399ff","#0000ff","#0099ff","#0099cc"))+
  theme_bw()+
  theme(
    plot.background = element_blank()
    ,panel.grid.major = element_blank()
    ,panel.grid.minor = element_blank()
    #,panel.border = element_blank()
    ,legend.position = "bottom"
  )+
  theme(axis.line = element_line(color = 'black', size = .9),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    plot.title = element_text(size = 20),
  )+
  labs(title =title,  element_text(size = 20), x ="90min", y = "60min")+
  scale_x_continuous(labels = scientific, breaks = c(0,200000,400000,600000, 800000, 1000000, 1200000, 1400000))+#, limits = c(0,5000))+
  scale_y_continuous(labels = scientific)#, breaks = c(0,5000000,7000000,9000000,10000000,30000000, 50000000, 70000000,90000000))# limits = c(0,5000))
#dev.off()
