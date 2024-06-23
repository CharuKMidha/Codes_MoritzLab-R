require(ISLR)
names(Smarket)
summary(Smarket)
?Smarket
Smarket[1:5,1:9]
pairs(Smarket,col=Smarket$Direction)
install.packages("devEMF")

penguins<- read.csv("scatterDIANN_SN3cells.txt", sep = "\t")
dim(penguins)
head(penguins)
tail(penguins)
penguins<- drop_na(penguins)
penguins$Protein= NULL
Direction<- penguins$Direction
penguins$Direction = NULL
penguins<-log10(penguins)
#pairs(penguins,col=Direction)
library(devEMF)
emf(file = "scaterplotHDbw.emf",emfPlus = FALSE)
?pairs
pairs(penguins,col="#ffa78c")
pairs(penguins)
dev.off()
