install.packages(c("tidyr", "rlang"))
library(tidyverse)
require(ggplot2)
#install.packages("cli")
#update.packages("cli")

df %>% 
  pivot_longer(cols = -User) %>%
  ggplot(aes(x= value, y= User, group = User, color = name)) +
  geom_line()+
  geom_point(size=4) +
  theme_classic()




df <- structure(list(User = c("A", "B", "C", "D", "E", "F"), Start = c(1L, 
                                                                       2L, 1L, 4L, 3L, 1L), End = c(3L, 5L, 3L, 5L, 4L, 5L)), 
                class = "data.frame", row.names = c(NA, -6L))




library(forestplot)
library(dplyr)
# Cochrane data from the 'rmeta'-package
cochrane_from_rmeta <- structure(list(mean  = c(NA, NA,-1.7202, -1.16955, -4.2275, -0.98505, -2.5082, -1.4535, -2.0459, -2.15645, -1.38365, -0.9646, -1.73975, -1.99395, -2.1198, 1.36255, -1.08125, -2.01415, -1.581625, -2.6978, -1.93785, -2.00035, -1.0402, -1.5441, -1.321, -2.1943, -1.7985, -1.52815, -2.2621, -1.3733, -0.88485), 
                                      lower = c(NA, NA, -3.2514, -2.7941, -5.055, -2.3111, -4.6074, -3.261, -3.1248, -4.6819, -2.8863, -2.4812, -3.8975, -3.1349, -3.3436, 2.1541, -3.2335, -5.0363, -3.55025, -4.6386, -3.0457, -3.7177, -2.8114, -2.7102, -2.297, -3.5036, -4.362, -3.6163, -3.9472, -3.2306, -2.4017),
                                      upper = c(NA, NA, -0.189, 0.455, -3.4, 0.341, -0.409, 0.354, -0.967, 0.369, 0.119, 0.552, 0.418, -0.853, -0.896, 0.571, 1.071, 1.008, 0.387, -0.757, -0.83, -0.283, 0.731, -0.378, -0.345, -0.885, 0.765, 0.56, -0.577, 0.484, 0.632)),
                                 .Names = c("mean", "lower", "upper"), 
                                 row.names = c(NA, -31L), 
                                 class = "data.frame")

tabletext <- cbind(c("", "Gene Name", "A1BG", "AMBP", "ANXA3", "APOH", "BCAM", "C4BPB", "CNDP1", "CRYAB", "DDAH2", "DKK3", "EFNB1", "FABP5", "FBLN5", "FEN1", "FLNC", "HBG1", "ITIH4", "ITIH5", "LAMA5", "LAMB2", "LBP", "LTBP4", "NFASC", "OTC", "PEBP4", "PTGDS", "PTPRM", "SAA2", "SLPI", "Gene Name", "A1BG", "AMBP", "ANXA3", "APOH", "BCAM", "C4BPB", "CNDP1", "CRYAB", "DDAH2", "DKK3", "EFNB1", "FABP5", "FBLN5", "FEN1", "FLNC", "HBG1", "ITIH4", "ITIH5", "LAMA5", "LAMB2", "LBP", "LTBP4", "NFASC", "OTC", "PEBP4", "PTGDS", "PTPRM", "SAA2", "SLPI", "Gene Name", "A1BG", "AMBP", "ANXA3", "APOH", "BCAM", "C4BPB", "CNDP1", "CRYAB", "DDAH2", "DKK3", "EFNB1", "FABP5", "FBLN5", "FEN1", "FLNC", "HBG1", "ITIH4", "ITIH5", "LAMA5", "LAMB2", "LBP", "LTBP4", "NFASC", "OTC", "PEBP4", "PTGDS", "PTPRM", "SAA2", "SLPI"))

cochrane_from_rmeta %>% 
  forestplot(labeltext = tabletext, 
             is.summary = c(rep(TRUE, 2), rep(FALSE, 8), TRUE),
             clip = c(0.1, 2.5), 
             xlog = TRUE, 
             col = fpColors(box = "royalblue",
                            line = "darkblue",
                            summary = "royalblue"))
# Create empty example plot
plot(0, 0, col = "white", xlab = "", ylab = "")       

# Create data frame with line-values
multiple_segments <- data.frame(x0 = c(2.1, 0.2, - 0.7, 0.4, - 0.8),   
                                y0 = c(0.8, 0.3, 0.5, - 0.4, 0.3))
                 #               x1 = c(0, 0.4, 0.5, - 0.5, - 0.7),
                  #              y1 = c(- 0.3, 0.4, - 0.5, - 0.7, 0.8))

# Draw multiple lines                       
segments(x0 = multiple_segments$x0,                   
         y0 = multiple_segments$y0)

#         x1 = multiple_segments$x1,
 #        y1 = multiple_segments$y1)

set.seed(1)
x <- c(2.1, 0.2, - 0.7, 0.4, - 0.8)
y <- c(0.8, 0.3, 0.5, - 0.4, 0.3),

# function
segmentInf <- function(xs, ys){
  fit <- lm(ys~xs)
  abline(fit)
}

plot(x,y)
segmentInf(x,y)
