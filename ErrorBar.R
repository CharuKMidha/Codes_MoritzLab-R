my.data = read.table("clipboard", head = T)
my.data
attach(my.data)
mean.worm = tapply(ValidCount, Columns, mean)
sd.worm = tapply(ValidCount, Columns, sd)
n.worm = tapply(ValidCount, Columns, length) #subjects per group
sem.worm = tapply(ValidCount, Columns, sd)/sqrt(tapply(ValidCount, Columns, length))#Standart error of the mean

mean.worm
sd.worm
n.worm
sem.worm

mids = barplot(mean.worm, xlab = "K562_Libraries",ylab= "# Peptides identified", ylim = c(0,40000), col = grey.colors(4,start=.3, end= .9, gamma = 2.2)) #col = rainbow(4)
arrows(mids, mean.worm - sem.worm, mids, mean.worm + sem.worm, code = 3, angle = 90, col = "blue")
#text(mids,2,paste("n= ",n.worm))
stripchart(ValidCount~Columns, data = my.data, vertical = TRUE, method= "jitter", pch = 21, col = 'red', bg = "yellow", add= TRUE)
