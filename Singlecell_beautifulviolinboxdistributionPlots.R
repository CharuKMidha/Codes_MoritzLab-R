install.packages("ggstatsplot")
library(ggstatsplot)
library(palmerpenguins)
library(tidyverse)
#data("penguins", package = "palmerpenguins")
penguins<- read.csv("proteotypicDIANN_SN.txt", sep = "\t")
dim(penguins)
head(penguins)
penguins<- drop_na(penguins)
penguins<- penguins[grep(pattern = 'C3_SN', penguins$cell), ]

penguins$log10val<-log10(penguins$quantity)
penguins$log2val<-log2(penguins$quantity)
?ggbetweenstats
plt <- ggbetweenstats(
  data = penguins,
  x = cell,
  y = log10val
)

plt <- plt + 
  ggplot2::scale_y_continuous(
    limits = c(0, 10),
    breaks = seq(from = 0, to = 10, by = 1))+
  # Add labels and title
  labs(
    x = "Cell types",
    y = "Quantities",
    title = "Distribution of quantities in 3cells from SN and DIANN"
  ) + 
  # Customizations
  theme(
    #This is the new default font in the plot
    text = element_text(family = "Roboto", size = 8, color = "black"),
    plot.title = element_text(
      family = "Lobster Two",
      size = 20,
      face = "bold",
      color = "#2a475e"
    ),
    # Statistical annotations below the main title
    # plot.subtitle = element_text(
    #   family = "Roboto",
    #   size = 15,
    #   face = "bold",
    #   color="#1b2838"
    # ),
    plot.title.position = "plot", # slightly different from default
    axis.text = element_text(size = 10, color = "black"),
    axis.title = element_text(size = 12)
  )
library(devEMF)
emf(file = "violinboxdistplotHDlog10SN_C3.emf",emfPlus = FALSE, width = 2, height= 3)

plt

dev.off()
# install.packages("here")
# require(here)
# 
# plt
# ggsave(
#   filename = here::here("img", "fromTheWeb", "web-violinplot-with-ggstatsplot.png"),
#   plot = plt,
#   width = 8,
#   height = 8,
#   device = "png"
# )
