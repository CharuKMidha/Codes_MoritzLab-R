
rm(list = ls())
my.data = read.table("clipboard", head = T)
my.data
attach(my.data)

#Proteins
my.data<- data.frame(C3=c(3051, 50), C1= c(2225, 116), C9= c(1548, 211))

#row_name<- my.data$SingelCellsPrecursors
row.names(my.data)<- c("Lib-based","DirectDIA")
                       
data<- as.matrix(my.data)
#attach(my.data)
# create color palette:
library(RColorBrewer)
coul <- brewer.pal(3, "Pastel2") 
head(my.data)
library(devEMF)
#emf(file = "Identificationprecursor.emf",emfPlus = FALSE)
emf(file = "Identificationprotein.emf",emfPlus = FALSE)

# Make a stacked barplot--> it will be in %!
barplot(data, col=coul , border="white", xlab="group")
dev.off()

