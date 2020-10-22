library(ggplot2)
#FILE 1

data = read.csv("2-3ChargeRT-densityPHL.txt", header = TRUE, sep ="\t")

data
dim(data)

#####Ploting 

library(ggplot2)
library(scales)
# Open a pdf file

#pdf("K562inhouserawMergerall_SN_Protein_matrix.pdf", width= 8, height = 8)
tiff("CorrelationRT23peptidesPHL-SN.tif", units = "in", width = 5,height = 5, res = 600)

dim(data)
N=NROW(data)
# Default way to draw Scatter Plot
#ggplot(data, aes(x = X.2RT, y = X.3RT)) + geom_point(color = "midnightblue") + geom_smooth()
correlation = cor(data$X.2RT , data$X.3RT, method = "pearson")
correlation = sprintf("%1f",correlation)
title = paste0("PHL +2/+3 pair RT correlation, n=",N, ", R2=",correlation)
title
ggplot(data, aes(x = data$X.2RT , y = data$X.3RT)) + 
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
  labs(title =title, x = "+2RT", y = "+3RT")
  #scale_x_continuous(labels = scientific)+#, limits = c(0,5000))+
  #scale_y_continuous(labels= scientific)# limits = c(0,5000))
dev.off()
