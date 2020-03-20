
library(ggplot2)

tiff("PHL_PV_q3bad0.05final.tif", units = "in", width = 8,height = 5, res = 300)

data = read.csv("K562_SWATH1_PHL_q3bad0.05_PV_score.txt", header = TRUE, sep= '\t')
#data = read.csv("a.txt", header = TRUE, sep= '\t')
head(data)
ggplot(data,aes(x=Score, fill = Assays)) + 
  geom_histogram(binwidth = .5, position= "dodge")+
  labs (x="Scores", y="Frequency")+
  theme(panel.background = element_rect(fill = 'white', colour = 'black'))+
  scale_x_continuous(breaks = c(-15,-10,-5,0,5,10,15,20))+#scale_x_continuous(limits = c(-15, 15))+
  scale_y_continuous(limits = c(0,25000))+
  geom_vline(xintercept = 2.769, linetype="dotted", color = "black", size=1)+ # 0.02=2.599, 0.05 = 2.769 , 0.01=2.518, Good= 2.44946
  theme(axis.line = element_line(color = 'black', size = .3))+
  labs(x = "Scores", y = "Frequency")
  
  
  theme(
    plot.background = element_blank()
    ,axis.line = element_blank()
  )
  
  

dev.off()
