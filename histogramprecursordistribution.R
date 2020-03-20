library(ggplot2)


#data = read.csv("precirsorHistogram.txt", header = TRUE, sep= '\t')

#data = read.csv("freqQ1Q3_freq.txt", header = TRUE, sep= '\t')
data = read.csv("PeptideLength.txt", header = TRUE, sep= '\t')


head(data$PeptideLength)

hist(data$PeptideLength)

hist(data$PeptideLength, breaks=50, xlim = c(5,50), ylim = c(0, 5000), 
     col = "royalblue", border ="cadetblue") #main="With breaks=100"

