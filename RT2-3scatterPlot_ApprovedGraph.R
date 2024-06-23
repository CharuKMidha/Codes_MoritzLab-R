library(ggplot2)
#FILE 1

data = read.csv("UndepletedMPlasmaDIA.txt", header = TRUE, sep ="\t")

data[1:5,1:5]
dim(data)
row.names(data)<- data$PG.ProteinGroups
data$PG.ProteinGroups<- NULL
#data[is.na(data)] <- 0

library(zoo)
data[] <- t(na.aggregate(t(data)))
data[1:5,1:5]


data$rmean<- rowMeans(data)

data2<- cbind.data.frame(data[20],data[21], data[22])
colnames(data2)<-c("protein", "sample1","sample2", "mean")

data2<- na.omit(data2)
data2[1:5,]
dim(data2)
nonproteotypic <- data2[ grepl( ";" , data2$protein ), ] #extract all males
proteotypic <- data2[!grepl( ";" , data2$protein ), ] #extract all males


#nonproteotypic<- 
tail(nonproteotypic)
dim(nonproteotypic)

####Ploting 

library(ggplot2)
library(scales)
# Open a pdf file

#pdf("K562inhouserawMergerall_SN_Protein_matrix.pdf", width= 8, height = 8)
#tiff("CorrelationRT23peptidesPHL-SN.tif", units = "in", width = 5,height = 5, res = 600)

dim(data)
N=NROW(data)
# Default way to draw Scatter Plot
#ggplot(data, aes(x = X.2RT, y = X.3RT)) + geom_point(color = "midnightblue") + geom_smooth()
correlation = cor(data$X2RT , data$X3RT, method = "pearson")
correlation = sprintf("%1f",correlation)
title = paste0("Ecoli +2/+3 pair RT correlation, n=",N, ", R2=",correlation)
title
ggplot(data, aes(x = data$X2RT , y = data$X3RT)) + 
  geom_point(shape = 1, size =8, color = '#99FF33') + ## #008080 teal
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
  theme(axis.line = element_line(color = 'black', size = .9),
        axis.title.x = element_text(size = 40),
        axis.title.y = element_text(size = 40),
        axis.text.x = element_text(size = 40),
        axis.text.y = element_text(size = 40))+
  labs(title =title, x = "+2RT (mins)", y = "+3RT (mins)", size = 40)
#scale_x_continuous(labels = scientific)+#, limits = c(0,5000))+
#scale_y_continuous(labels= scientific)# limits = c(0,5000))
#dev.off()



