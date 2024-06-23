# Data Preprocessing
rm(list = ls())

require(dplyr)
require(tibble)
#remove.packages("ggplot2")
#install.packages('ggplot2', dependencies = TRUE)
require(ggplot2)
#install.packages('tidyverse', dependencies = TRUE)

library(tidyverse)#for str_replace()

# these are from Bioconductor
library(limma) 
library(edgeR) 
library(sva)

# we use this for the pairwise correlation plots
library(psych)
# Importing the dataset

#dataset <- read.csv("z:/cmidha/Longivity_imputation/MouseM005ck/Mouse_2021_10_Liver_protein_unnormalized_proteotypic-Phase2-append-repeat-samples.xls", header = TRUE, sep ="\t")    
#SPN16
dataset <- read.csv("//harp/proteomics/cmidha/Borrelia/4Helisa/Osm28/proteinMEAN2304x19Osm28.csv")    



head(dataset[1:5,1:5])
# peptide<- paste0(dataset$Peptide,"#",dataset$Protein)
# dataset$Peptide<- NULL
# dataset$Protein<- NULL
dim(dataset)#6005 proteins X 304 samples

head(dataset[1:5,1:5])
dataset_nocontam<- (dataset[!grepl("CONTAM.*", dataset$protein),])
dim(dataset_nocontam)
row.names(dataset_nocontam)<- dataset_nocontam$protein 
head (dataset_nocontam)
dataset_nocontam$protein<- NULL

# boxplot(dataset_nocontam, col = c(rep("red", 6), rep("blue",6), rep("green",3), rep("orange",3)), ylim= c(0,10000000),
#           names = c(rep("biotin", 6), rep("CONTROL",6), rep("HEKWTB",3), rep("HEKWTc",3)))
library(devEMF)
emf(file = "Rawlog10boxplot.emf")
boxplot(log(dataset_nocontam,10), 
        col = c(rep("red", 6), rep("blue",6), rep("green",3), rep("orange",3)),
        ylab  ="MS1 abundance in log 10" , 
        par(cex.lab=1.5), # is for y-axis
        par(cex.axis=1.5), # is for x-axis
        names = c(rep("biotin", 6), rep("CONTROL",6), rep("HEKWTB",3), rep("HEKWTc",3)))

dev.off()

biotin<- dataset_nocontam[ , grepl( "biotin.*" , names( dataset_nocontam ) ) ]
control1<-dataset_nocontam[ , grepl( "control.*" , names( dataset_nocontam ) ) ]
control2<- dataset_nocontam[ , grepl( "hekwt_b.*" , names( dataset_nocontam ) ) ]
control3<- dataset_nocontam[ , grepl( "hekwt_c.*" , names( dataset_nocontam ) ) ]

dim(biotin)
dim(control1)
dim(control2)
dim(control3)
head (biotin)
head(control1)




biotinvsc1<- cbind.data.frame(biotin, control1)
biotinvsc2<- cbind.data.frame(biotin, control2)
biotinvsc3<- cbind.data.frame(biotin, control3)


dim(biotinvsc1)
dim(biotinvsc2)
dim(biotinvsc3)



#Sample Loading Normalization:-Sample loading (SL) normalization is a grand 
#total scaling method. We will scale to column totals to something in the middle.
#That keeps the new values in the same ballpark as the original values. 
#Taking the arithmetic average of the column totals and doing that experiment-wide 
#(for all conditions). 


# function to do grand total (sample loading) normalization
SL_norm <- function(df, print_factors = TRUE) {
  df <- na.omit(df)
  # compute norm factors and scale columns
  norm_facs <- mean(colSums(df)) / colSums(df)
  df_sl  <- sweep(df, 2, norm_facs, FUN = "*")
  
  # print the normalization factors for QC check
  if (print_factors == TRUE) {
    cat(dim(df),"\n")
    cat("\nSample Loading Normalization Factors:\n ")
    cat(sprintf("%s - %0.3f\n", colnames(df), norm_facs))
  }
  
  df_sl # return normalized data
}



# SL norm the entire data subset


# biotin_sl <- SL_norm(data.frame(biotin))#, control1, control2, control3))
# control1_sl <- SL_norm(data.frame(control1))
# control2_sl <- SL_norm(data.frame(control2))
# control3_sl <- SL_norm(data.frame(control3))

biotinvsc1_sl <- SL_norm(data.frame(biotinvsc1))
biotinvsc2_sl <- SL_norm(data.frame(biotinvsc2))
biotinvsc3_sl <- SL_norm(data.frame(biotinvsc3))



# boxplot(c(log(biotinvsc1_sl,10),log(biotinvsc2_sl,10) ,log(biotinvsc3_sl,10)) ,
#         col = c(rep("red", 12), rep("blue",9), rep("green",9)),
#         ylab  ="MS1 abundance in log 10" , 
#         par(cex.lab=1.5), # is for y-axis
#         par(cex.axis=1), # is for x-axis
#         names = c(rep("biotinvsC1_SL", 12), rep("biotinvsC2_SL",9), rep("biotinvsC3_SL",9)))
# 
# 
# boxplot(c(log(biotinvsc1,10),log(biotinvsc2,10) ,log(biotinvsc3,10)) ,
#         col = c(rep("red", 12), rep("blue",9), rep("green",9)),
#         ylab  ="MS1 abundance in log 10" , 
#         par(cex.lab=1.5), # is for y-axis
#         par(cex.axis=1), # is for x-axis
#         names = c(rep("biotinvsC1", 12), rep("biotinvsC2",9), rep("biotinvsC3",9)))
# 
# 
# boxplot(c(biotinvsc1_sl,biotinvsc2_sl ,biotinvsc3_sl) ,
#         col = c(rep("red", 12), rep("blue",9), rep("green",9)),
#         ylab  ="MS1 abundance raw SL" , 
#         ylim = c(0, 10000000),
#         par(cex.lab=1.5), # is for y-axis
#         par(cex.axis=1), # is for x-axis
#         names = c(rep("biotinvsC1_SL", 12), rep("biotinvsC2_SL",9), rep("biotinvsC3_SL",9)))
# 
# 
# boxplot(c(biotinvsc1,biotinvsc2 ,biotinvsc3) ,
#         col = c(rep("red", 12), rep("blue",9), rep("green",9)),
#         ylab  ="MS1 abundance raw" ,
#         ylim = c(0, 10000000),
#         par(cex.lab=1.5), # is for y-axis
#         par(cex.axis=1), # is for x-axis
#         names = c(rep("biotinvsC1", 12), rep("biotinvsC2",9), rep("biotinvsC3",9)))



##RAW vs SL in log 10 comparison
emf(file = "Rawlog10vsSLboxplotC1.emf")
boxplot(c(log(biotinvsc1,10),log(biotinvsc1_sl,10)) ,
        col = c(rep("red", 6), rep("blue",6),rep("green", 6), rep("orange",6)),
        ylab  ="MS1 abundance in log 10" , 
        par(cex.lab=1.5), # is for y-axis
        par(cex.axis=1), # is for x-axis
        main = "Abundance log10: red(biotin),blue(c1), SL: green(biotin),orange(c1)",
        names = c(rep("biotin", 6), rep("Control1", 6),rep("biotin_SL",6), rep("control1_SL",6)))
dev.off()
emf(file = "Rawlog10vsSLboxplotC2.emf")
boxplot(c(log(biotinvsc2,10),log(biotinvsc2_sl,10)) ,
        col = c(rep("red", 6), rep("blue",3),rep("green", 6), rep("orange",3)),
        ylab  ="MS1 abundance in log 10" , 
        par(cex.lab=1.5), # is for y-axis
        par(cex.axis=1), # is for x-axis
        main = "Abundance log10: red(biotin),blue(c2), SL: green(biotin),orange(c2)",
        names = c(rep("biotin", 6), rep("Control2", 3),rep("biotin_SL",6), rep("control2_SL",3)))
dev.off()
emf(file = "Rawlog10vsSLboxplotC3.emf")
boxplot(c(log(biotinvsc3,10),log(biotinvsc3_sl,10)) ,
        col = c(rep("red", 6), rep("blue",3),rep("green", 6), rep("orange",3)),
        ylab  ="MS1 abundance in log 10" , 
        par(cex.lab=1.5), # is for y-axis
        par(cex.axis=1), # is for x-axis
        main = "Abundance log10: red(biotin),blue(c3), SL: green(biotin),orange(c3)",
        names = c(rep("biotin", 6), rep("Control3", 3),rep("biotin_SL",6), rep("control3_SL",3)))

dev.off()
format(round(colSums(biotinvsc1_sl), digits = 0), big.mark = ",") #ALL sum sould be same
format(round(colSums(biotinvsc2_sl), digits = 0), big.mark = ",")
format(round(colSums(biotinvsc3_sl), digits = 0), big.mark = ",")

#TMM NORMLAIZATION

sl_tmm <- calcNormFactors(biotinvsc1_sl)
biotinvsc1_sl_tmm <- sweep(biotinvsc1_sl, 2, sl_tmm, FUN = "/") # data after SL and TMM on original scale

sl_tmm <- calcNormFactors(biotinvsc2_sl)
biotinvsc2_sl_tmm <- sweep(biotinvsc2_sl, 2, sl_tmm, FUN = "/") # data after SL and TMM on original scale

sl_tmm <- calcNormFactors(biotinvsc3_sl)
biotinvsc3_sl_tmm <- sweep(biotinvsc3_sl, 2, sl_tmm, FUN = "/") # data after SL and TMM on original scale


write.csv(biotinvsc1_sl_tmm, "biotinvsc1_sl_tmm.csv")
write.csv(biotinvsc2_sl_tmm, "biotinvsc2_sl_tmm.csv")
write.csv(biotinvsc3_sl_tmm, "biotinvsc3_sl_tmm.csv")


##RAW vs SL & TMM in log 10 comparison
emf(file = "Rawlog10vsSLTMMboxplotC1.emf")
boxplot(c(log(biotinvsc1,10),log(biotinvsc1_sl_tmm,10)) ,
        col = c(rep("red", 6), rep("blue",6),rep("green", 6), rep("orange",6)),
        ylab  ="MS1 abundance in log 10" , 
        par(cex.lab=1.5), # is for y-axis
        par(cex.axis=1), # is for x-axis
        main = "Abundance log10: red(biotin),blue(c1), SL_TMM: green(biotin),orange(c1)",
        names = c(rep("biotin", 6), rep("Control1", 6),rep("biotin_SL_TMM",6), rep("control1_SL_TMM",6)))

dev.off()

emf(file = "Rawlog10vsSLTMMboxplotC2.emf")
boxplot(c(log(biotinvsc2,10),log(biotinvsc2_sl_tmm,10)) ,
        col = c(rep("red", 6), rep("blue",3),rep("green", 6), rep("orange",3)),
        ylab  ="MS1 abundance in log 10" , 
        par(cex.lab=1.5), # is for y-axis
        par(cex.axis=1), # is for x-axis
        main = "Abundance log10: red(biotin),blue(c2), SL_TMM: green(biotin),orange(c2)",
        names = c(rep("biotin", 6), rep("Control2", 3),rep("biotin_SL_TMM",6), rep("control2_SL_TMM",3)))
dev.off()

emf(file = "Rawlog10vsSLTMMboxplotC3.emf")
boxplot(c(log(biotinvsc3,10),log(biotinvsc3_sl_tmm,10)) ,
        col = c(rep("red", 6), rep("blue",3),rep("green", 6), rep("orange",3)),
        ylab  ="MS1 abundance in log 10" , 
        par(cex.lab=1.5), # is for y-axis
        par(cex.axis=1), # is for x-axis
        main = "Abundance log10: red(biotin),blue(c3), SL_TMM: green(biotin),orange(c3)",
        names = c(rep("biotin", 6), rep("Control3", 3),rep("biotin_SL_TMM",6), rep("control3_SL_TMM",3)))
dev.off()


dim(biotinvsc1_sl_tmm)
dim(biotinvsc2_sl_tmm)
dim(biotinvsc3_sl_tmm)

biotinvsc1_sl_tmm_biotin<- rowMeans(biotinvsc1_sl_tmm[,1:6], na.rm = TRUE)
biotinvsc1_sl_tmm_c1<- rowMeans(biotinvsc1_sl_tmm[,7:12], na.rm = TRUE)

biotinvsc2_sl_tmm_biotin<- rowMeans(biotinvsc2_sl_tmm[,1:6], na.rm = TRUE)
biotinvsc2_sl_tmm_c2<- rowMeans(biotinvsc2_sl_tmm[,7:9], na.rm = TRUE)

biotinvsc3_sl_tmm_biotin<- rowMeans(biotinvsc3_sl_tmm[,1:6], na.rm = TRUE)
biotinvsc3_sl_tmm_c3<- rowMeans(biotinvsc3_sl_tmm[,7:9], na.rm = TRUE)

LFbiotinvsc1_sl_tmm = log(biotinvsc1_sl_tmm_biotin/biotinvsc1_sl_tmm_c1, 2)
LFbiotinvsc2_sl_tmm = log(biotinvsc2_sl_tmm_biotin/biotinvsc2_sl_tmm_c2, 2)
LFbiotinvsc3_sl_tmm = log(biotinvsc3_sl_tmm_biotin/biotinvsc3_sl_tmm_c3, 2)





# log(10/2,2)
# log(10,2)-log(2,2)
#####################################
##DE proteins/genes biotin vs Control1 
#####################################

condition<- "control"
trimdata0_CON <- biotinvsc1_sl_tmm[ ,grepl(condition , colnames(biotinvsc1_sl_tmm) ) ]
dim(trimdata0_CON)

condition<- "biotin"
trimdata0_biotin <- biotinvsc1_sl_tmm[ ,grepl(condition , colnames(biotinvsc1_sl_tmm) ) ]
dim(trimdata0_biotin)

trimdata0convsbiotin<- cbind.data.frame(trimdata0_CON, trimdata0_biotin)



# CRvsC<- cbind.data.frame(UMHET3MCR, UMHET3MC)
# head (CRvsC)
# dim(CRvsC)


# dim(UMHET3MCR)
# dim(UMHET3MC)
# group = factor(rep(c("CR","C"),c(9,10)))

group_trimdata0 = factor(rep(c("CON","biotin"),c(6,6)))


# res = rowttests(as.matrix(CRvsC),factor(group))
# res$adjP = p.adjust(res$p.value,"BH")
# head(res)



######Assuming there is difference in variance in two groups

library(broom)
res_trimdata0 = apply(biotinvsc1_sl_tmm,1,function(i)tidy(t.test(i ~ group_trimdata0)))



res_trimdata0 = do.call(rbind,res_trimdata0)
res_trimdata0$adjP = p.adjust(res_trimdata0$p.value,"BH")




head(res_trimdata0)
#res1[,c("statistic","p.value","adjP")]

#row.names(res1)<- rownames(FCR)
rownames(res_trimdata0)<- rownames(biotinvsc1_sl_tmm)



# head(res1)
# 
# write.csv(res1,"welch2samplettestGHRCOvsCfemalec57bl6.csv", row.names = TRUE, col.names = TRUE)


#####################################
##DE proteins/genes CR vs C VOLVANO PLOT
#####################################
#DE_CRvsC<- cbind.data.frame(rownames(CRvsC),CRlog2FCratio, res$adjP)
DEtrimdata0df<-cbind.data.frame(rownames(biotinvsc1_sl_tmm),LFbiotinvsc1_sl_tmm ,res_trimdata0$adjP)


colnames(DEtrimdata0df)<- c("protein", "log2FC", "AdjPvalue")


write.csv(DEtrimdata0df,"DEbiotinvsc1_sl_tmm.csv")
-log(0.05,10)



# 
# if (!requireNamespace('BiocManager', quietly = TRUE))
#   install.packages('BiocManager')
# 
# BiocManager::install('EnhancedVolcano')
library(EnhancedVolcano)
#

#
# EnhancedVolcano(DEtrimdata0df,
#                 lab = DEtrimdata0df$protein ,
#                 x = 'log2FC',
#                 y = 'AdjPvalue')

proteingene<- DEtrimdata0df$protein

DEtrimdata0df$protein <- gsub("sp.", "", DEtrimdata0df$protein)
DEtrimdata0df$protein <- gsub("_HUMAN", "", DEtrimdata0df$protein)



VOLCANO<- function(DE,conditionname){
  EnhancedVolcano(DE,
                  lab = DE$protein,
                  x = 'log2FC',
                  y = 'AdjPvalue',
                  title = paste0(conditionname," versus Control"),
                  pCutoff = .05,
                  FCcutoff = 2,
                  pointSize = 5,
                  labSize = 6,
                  ylim = c(0, 6)
                  #xlim = c(-1, 1)
                  )
}
emf(file = "biotinC1SL_TMM_noimputVplot.emf")
VOLCANO(DEtrimdata0df,"Biotin-C1 SL & TMM LogFC >2, adj.pval <.05")
dev.off()




##################################################

#####################################
##DE proteins/genes biotin vs Control2 
#####################################

condition<- "hekwt_b"
trimdata0_CON <- biotinvsc2_sl_tmm[ ,grepl(condition , colnames(biotinvsc2_sl_tmm) ) ]
dim(trimdata0_CON)

condition<- "biotin"
trimdata0_biotin <- biotinvsc2_sl_tmm[ ,grepl(condition , colnames(biotinvsc2_sl_tmm) ) ]
dim(trimdata0_biotin)

trimdata0convsbiotin<- cbind.data.frame(trimdata0_CON, trimdata0_biotin)



# CRvsC<- cbind.data.frame(UMHET3MCR, UMHET3MC)
# head (CRvsC)
# dim(CRvsC)


# dim(UMHET3MCR)
# dim(UMHET3MC)
# group = factor(rep(c("CR","C"),c(9,10)))

group_trimdata0 = factor(rep(c("CON","biotin"),c(3,6)))


# res = rowttests(as.matrix(CRvsC),factor(group))
# res$adjP = p.adjust(res$p.value,"BH")
# head(res)



######Assuming there is difference in variance in two groups

library(broom)
res_trimdata0 = apply(biotinvsc2_sl_tmm,1,function(i)tidy(t.test(i ~ group_trimdata0)))



res_trimdata0 = do.call(rbind,res_trimdata0)
res_trimdata0$adjP = p.adjust(res_trimdata0$p.value,"BH")




head(res_trimdata0)
#res1[,c("statistic","p.value","adjP")]

#row.names(res1)<- rownames(FCR)
rownames(res_trimdata0)<- rownames(biotinvsc2_sl_tmm)



# head(res1)
# 
# write.csv(res1,"welch2samplettestGHRCOvsCfemalec57bl6.csv", row.names = TRUE, col.names = TRUE)


#####################################
##DE proteins/genes CR vs C VOLVANO PLOT
#####################################
#DE_CRvsC<- cbind.data.frame(rownames(CRvsC),CRlog2FCratio, res$adjP)
DEtrimdata0df<-cbind.data.frame(rownames(biotinvsc2_sl_tmm),LFbiotinvsc2_sl_tmm ,res_trimdata0$adjP)


colnames(DEtrimdata0df)<- c("protein", "log2FC", "AdjPvalue")


write.csv(DEtrimdata0df,"DEbiotinvsc2_sl_tmm.csv")




# 
# if (!requireNamespace('BiocManager', quietly = TRUE))
#   install.packages('BiocManager')
# 
# BiocManager::install('EnhancedVolcano')
#library(EnhancedVolcano)
#

#
# EnhancedVolcano(DEtrimdata0df,
#                 lab = DEtrimdata0df$protein ,
#                 x = 'log2FC',
#                 y = 'AdjPvalue')

proteingene<- DEtrimdata0df$protein

DEtrimdata0df$protein <- gsub("sp.", "", DEtrimdata0df$protein)
DEtrimdata0df$protein <- gsub("_HUMAN", "", DEtrimdata0df$protein)



VOLCANO<- function(DE,conditionname){
  EnhancedVolcano(DE,
                  lab = DE$protein,
                  x = 'log2FC',
                  y = 'AdjPvalue',
                  title = paste0(conditionname," versus Control"),
                  pCutoff = .05,
                  FCcutoff = 2,
                  pointSize = 5,
                  labSize = 6,
                  ylim = c(0, 2)
                  #xlim = c(-1, 1)
  )
}
emf(file = "biotinC2SL_TMM_noimputVplot.emf")
VOLCANO(DEtrimdata0df,"Biotin-C2 SL & TMM LogFC >2, adj.pval <.05")
dev.off()





##################################################

#####################################
##DE proteins/genes biotin vs Control3 
#####################################

condition<- "hekwt_c"
trimdata0_CON <- biotinvsc3_sl_tmm[ ,grepl(condition , colnames(biotinvsc3_sl_tmm) ) ]
dim(trimdata0_CON)

condition<- "biotin"
trimdata0_biotin <- biotinvsc3_sl_tmm[ ,grepl(condition , colnames(biotinvsc3_sl_tmm) ) ]
dim(trimdata0_biotin)

trimdata0convsbiotin<- cbind.data.frame(trimdata0_CON, trimdata0_biotin)



# CRvsC<- cbind.data.frame(UMHET3MCR, UMHET3MC)
# head (CRvsC)
# dim(CRvsC)


# dim(UMHET3MCR)
# dim(UMHET3MC)
# group = factor(rep(c("CR","C"),c(9,10)))

group_trimdata0 = factor(rep(c("CON","biotin"),c(3,6)))


# res = rowttests(as.matrix(CRvsC),factor(group))
# res$adjP = p.adjust(res$p.value,"BH")
# head(res)



######Assuming there is difference in variance in two groups

library(broom)
res_trimdata0 = apply(biotinvsc3_sl_tmm,1,function(i)tidy(t.test(i ~ group_trimdata0)))



res_trimdata0 = do.call(rbind,res_trimdata0)
res_trimdata0$adjP = p.adjust(res_trimdata0$p.value,"BH")




head(res_trimdata0)
#res1[,c("statistic","p.value","adjP")]

#row.names(res1)<- rownames(FCR)
rownames(res_trimdata0)<- rownames(biotinvsc3_sl_tmm)



# head(res1)
# 
# write.csv(res1,"welch2samplettestGHRCOvsCfemalec57bl6.csv", row.names = TRUE, col.names = TRUE)


#####################################
##DE proteins/genes CR vs C VOLVANO PLOT
#####################################
#DE_CRvsC<- cbind.data.frame(rownames(CRvsC),CRlog2FCratio, res$adjP)
DEtrimdata0df<-cbind.data.frame(rownames(biotinvsc3_sl_tmm),LFbiotinvsc3_sl_tmm ,res_trimdata0$adjP)


colnames(DEtrimdata0df)<- c("protein", "log2FC", "AdjPvalue")


write.csv(DEtrimdata0df,"DEbiotinvsc3_sl_tmm.csv")




# 
# if (!requireNamespace('BiocManager', quietly = TRUE))
#   install.packages('BiocManager')
# 
# BiocManager::install('EnhancedVolcano')
#library(EnhancedVolcano)
#

#
# EnhancedVolcano(DEtrimdata0df,
#                 lab = DEtrimdata0df$protein ,
#                 x = 'log2FC',
#                 y = 'AdjPvalue')

proteingene<- DEtrimdata0df$protein

DEtrimdata0df$protein <- gsub("sp.", "", DEtrimdata0df$protein)
DEtrimdata0df$protein <- gsub("_HUMAN", "", DEtrimdata0df$protein)



VOLCANO<- function(DE,conditionname){
  EnhancedVolcano(DE,
                  lab = DE$protein,
                  x = 'log2FC',
                  y = 'AdjPvalue',
                  title = paste0(conditionname," versus Control"),
                  pCutoff = .05,
                  FCcutoff = 1,
                  pointSize = 5,
                  labSize = 6,
                  ylim = c(0, 2)
                  #xlim = c(-1, 1)
  )
}

emf(file = "biotinC3SL_TMM_noimputVplot.emf")
VOLCANO(DEtrimdata0df,"Biotin-C3 SL & TMM LogFC >2, adj.pval <.05")
dev.off()


##################################################
#**************STARTING IMPUTATION***************#

##Do the category wise missing analysis

dim(biotinvsc1)
dim(biotinvsc2)
dim(biotinvsc3)
head (biotinvsc1)
dim(dataset_nocontam)


REMOVEEMPTYROWS<-function(data){
  data[data==0] <- NA
  data<-data[rowSums(is.na(data)) != ncol(data), ]
  #head (data)
  #dim(data)
  return(data)
}


###REMOVING EMPTY ROWS
#proteotypic2<- REMOVEEMPTYROWS(proteotypic)
allcondt<- REMOVEEMPTYROWS(dataset_nocontam)
dim(allcondt)
allcondt[1:5,1:5]

biotinvsc1_noemptyrow<- REMOVEEMPTYROWS(biotinvsc1)
biotinvsc2_noemptyrow<- REMOVEEMPTYROWS(biotinvsc2)
biotinvsc3_noemptyrow<- REMOVEEMPTYROWS(biotinvsc3)
dim(biotinvsc1_noemptyrow)
dim(biotinvsc2_noemptyrow)
dim(biotinvsc3_noemptyrow)
biotinvsc1_noemptyrow[1:5,1:5]



#############COMPLETE DATA MISSINGNESS
##NOW focus on control groups and impute only 20% missing data

condition<- "control_"
control1_4impute <- biotinvsc1_noemptyrow[ ,grepl(condition , colnames(biotinvsc1_noemptyrow) ) ]
dim(control1_4impute)
head (control1_4impute)

boxplot (log(control1_4impute,10), ylim = c(4,9), main = "control1",
         par(cex.lab=3), # is for y-axis
         par(cex.axis=2), # is for x-axis)
         ylab = "abundance in log10"
)

# boxplot (log(biotin_4impute[,1:6],10), ylim = c(4,9), main = "biotin",
#          par(cex.lab=3), # is for y-axis
#          par(cex.axis=2),
#          ylab = "abundance in log10"# is for x-axis)
# )



condition<- "biotin_"
biotinc1_4impute <- biotinvsc1_noemptyrow[ ,grepl(condition , colnames(biotinvsc1_noemptyrow) ) ]
dim(biotinc1_4impute)
head (biotinc1_4impute)

#MISSING DIAGNOSIS biotinc1

biotinc1_4impute$na_count <-rowSums(is.na(biotinc1_4impute)) #counts missing values in a rows
head(biotinc1_4impute[1:5,])
head(biotinc1_4impute$na_count)
countNA <- data.frame(cbind(rownames(biotinc1_4impute), biotinc1_4impute$na_count,
                            round(biotinc1_4impute$na_count *100/(ncol(biotinc1_4impute)-1),1)))#-2 is protein and nacount column,-1 if protein name is row name 
head(countNA)
colnames(countNA)<- c("Protein", "Count", "Percentage")
countNA[countNA==0] <- NA #0s are complete rows
head(countNA$Count)
head(countNA$Percentage)
tail(countNA)
dim(countNA)
countNA2<- na.omit(countNA) #removed complete matrix protein info
nrowcountNA<- nrow(countNA)- nrow(countNA2)
dim(countNA2)# 

write.csv(countNA2, "countNAbiotinc1withMissinginfo.txt", row.names = TRUE)


write.csv(countNA, "countNAbiotinc1.txt", row.names = TRUE)

biotinc1_4impute[1:5,]
biotinc1_4impute$na_count<-NULL
dim(biotinc1_4impute)
dim(countNA)
dim(countNA2)
countNA[is.na(countNA)] <- 0 #make all the complete proteins to 0%missing and retained the protein info

#############COMPLETE DATA MISSINGNESS ENDS####################


TRIMDATA<- function(data1,data2,p) #data1=complete matrix data2=countNA, p = cutoff value to trim data
{
  data2<- data.frame(data2)
  data2$Count<- as.integer(data2$Count)
  data2$Percentage<- as.integer(data2$Percentage)
  #write.csv(data2, "missing.txt", row.names = TRUE)
  ############PLOTTING HISTOGRAM
  hist(data2$Percentage, col= "#b5ede6", border= "#249487", main = "Missing data distribution", labels = TRUE, xlab = "Protein missing data(%)", ylab = "Frequency", cex.lab= 2, cex.axis=2, cex.main=2, cex.sub=2)
  data3<-subset(data2, data2$Percentage<= p)
  hist(data3$Percentage, col= "#b5ede6", border= "#249487", main = "Missing data distribution after removing 20%", labels = TRUE, xlab = "Protein missing data(%)", ylab = "Frequency", cex.lab= 2, cex.axis=2, cex.main=2, cex.sub=2)
  proteins<- data3$Protein
  #data4<- subset(data1, data1[,1] %in% proteins)# first column is not protein names it is stored as rownames
  data4<- subset(data1, rownames(data1) %in% proteins)
  return(data4)
}


#TRANSPIOSE DATA
# 
# TRANSPOSEDATA<-function(data){
#   colproteinsname<-data[,1]
#   rowsamples = names(data[,2:ncol(data)])
#   data<- data[,2:ncol(data)]
#   rownames(data)=c()
#   colnames(data)=c()
#   data
#   tdata<- t(data)
#   tdata <- data.frame(tdata)
#   dim(tdata)
#   tdata
#   colnames(tdata)<-colproteinsname
#   tdata
#   tdata<- cbind(rowsamples, tdata)
#   tdata
#   return(tdata)
# }


trimdata20<-TRIMDATA(biotinc1_4impute, countNA, 20)
head(trimdata20)
dim(trimdata20)
#write.csv(trimdata20,"//harp/proteomics/cmidha/Borrelia/4Helisa/biotin/unnorm.csv", row.names = TRUE, col.names = TRUE)

##################################
biotinc1_4impute_20<- trimdata20 #
##################################


dim(biotinc1_4impute_20) #655

data20missbiotinc1<- subset(biotinvsc1_noemptyrow, rownames(biotinvsc1_noemptyrow) %in% rownames(biotinc1_4impute_20))
data20missbiotinc2<- subset(biotinvsc2_noemptyrow, rownames(biotinvsc2_noemptyrow) %in% rownames(biotinc1_4impute_20))
data20missbiotinc3<- subset(biotinvsc3_noemptyrow, rownames(biotinvsc3_noemptyrow) %in% rownames(biotinc1_4impute_20))


dim(data20missbiotinc1)
dim(data20missbiotinc2)
dim(data20missbiotinc3)

##IMPUTATION
library(limma)
library(pcaMethods)
library(VIM)


# IMPUTE 20%missing biotin

biotin_20<- data20missbiotinc1[1:6]

biotin_20KNN <- kNN(biotin_20, variable = colnames(biotin_20),k=10,  dist_var = colnames(biotin_20), numFun = median)

dim(biotin_20KNN)
head(biotin_20KNN[1:5,1:5])
imputeddata<- biotin_20KNN[,1:(ncol(biotin_20KNN)/2)]
rownames(imputeddata)<-  rownames(biotin_20)                      
head(imputeddata[1:5,])


biotin_20KNNimputed<- imputeddata



# kNN IMPUTE corresponding control1 of 20%missing biotin

control1_20<- data20missbiotinc1[7:12]
control1_20<- REMOVEEMPTYROWS(control1_20)

dim(control1_20) #367 LEFT

control1_20KNN <- kNN(control1_20, variable = colnames(control1_20),k=10,  dist_var = colnames(control1_20), numFun = median)

dim(control1_20KNN)
head(control1_20KNN[1:5,1:5])
imputeddata<- control1_20KNN[,1:(ncol(control1_20KNN)/2)]
rownames(imputeddata)<-  rownames(control1_20)                      
head(imputeddata[1:5,])


CONTROL1_20KNNimputed<- imputeddata

#AVERAGE VALUE IMPUTE IN CONTROLS ONLY 1
library(zoo)
control1_20rowavgimpute <- t(na.aggregate(t(control1_20)))
control1_20rowavgimpute

# kNN IMPUTE corresponding control2 of 20%missing biotin

control2_20<- data20missbiotinc2[7:9]
control2_20<- REMOVEEMPTYROWS(control2_20)

dim(control2_20) #232 LEFT

control2_20KNN <- kNN(control2_20, variable = colnames(control2_20),k=10,  dist_var = colnames(control2_20), numFun = median)

dim(control2_20KNN)
head(control2_20KNN[1:5,1:5])
imputeddata<- control2_20KNN[,1:(ncol(control2_20KNN)/2)]
rownames(imputeddata)<-  rownames(control2_20)                      
head(imputeddata[1:5,])


CONTROL2_20KNNimputed<- imputeddata

#AVERAGE VALUE IMPUTE IN CONTROLS ONLY 1
#library(zoo)
control2_20rowavgimpute <- t(na.aggregate(t(control2_20)))
control2_20rowavgimpute

# kNN IMPUTE corresponding control3 of 20%missing biotin

control3_20<- data20missbiotinc3[7:9]
control3_20<- REMOVEEMPTYROWS(control3_20)

dim(control3_20) #221 LEFT 

control3_20KNN <- kNN(control3_20, variable = colnames(control3_20),k=10,  dist_var = colnames(control3_20), numFun = median)

dim(control3_20KNN)
head(control3_20KNN[1:5,1:5])
imputeddata<- control3_20KNN[,1:(ncol(control3_20KNN)/2)]
rownames(imputeddata)<-  rownames(control3_20)                      
head(imputeddata[1:5,])


CONTROL3_20KNNimputed<- imputeddata

#AVERAGE VALUE IMPUTE IN CONTROLS ONLY 1
#library(zoo)
control3_20rowavgimpute <- t(na.aggregate(t(control3_20)))
control3_20rowavgimpute

####RETERIVING COMPLETE SETS

head(biotin_20KNNimputed)
head(CONTROL1_20KNNimputed)
head(CONTROL2_20KNNimputed)
head(CONTROL3_20KNNimputed)


head(control1_20rowavgimpute)
head(control2_20rowavgimpute)
head(control3_20rowavgimpute)

dim(biotin_20KNNimputed)
dim(CONTROL1_20KNNimputed)
dim(CONTROL2_20KNNimputed)
dim(CONTROL3_20KNNimputed)

dim(control1_20rowavgimpute)
dim(control2_20rowavgimpute)
dim(control3_20rowavgimpute)


##COMPARE AVERAGE IMPUTE AND kNN IMPUTE
C1kNNbyavg<- rowMeans(CONTROL1_20KNNimputed)/rowMeans(control1_20rowavgimpute)
C2kNNbyavg<- rowMeans(CONTROL2_20KNNimputed)/rowMeans(control2_20rowavgimpute)
C3kNNbyavg<- rowMeans(CONTROL3_20KNNimputed)/rowMeans(control3_20rowavgimpute)
C1kNNbyavg<- sort(C1kNNbyavg, decreasing = FALSE)
C2kNNbyavg<- sort(C2kNNbyavg, decreasing = FALSE)
C3kNNbyavg<- sort(C3kNNbyavg, decreasing = FALSE)
plot(C1kNNbyavg)
plot(C2kNNbyavg)
plot(C3kNNbyavg)

hist(C1kNNbyavg)
hist(C2kNNbyavg)
hist(C3kNNbyavg)

plot(C1kNNbyavg, type='p', col= "#E15759", cex= 3, lwd= 2, xlab="Number of protein", 
     ylab='Ratio of kNN / avg impute values', cex.axis=2, cex.lab=2)
points(C2kNNbyavg, type='p', cex= 3, lwd= 2,col = "#F28E2B") ##blueviolet", "", "deepskyblue", "powderblue", "slateblue"
points(C3kNNbyavg, type='p', cex= 3, lwd= 2,col = "#4E79A7")

C1KNNavgimputed20<- cbind.data.frame(rowMeans(CONTROL1_20KNNimputed), rowMeans(control1_20rowavgimpute))
C2KNNavgimputed20<- cbind.data.frame(rowMeans(CONTROL2_20KNNimputed), rowMeans(control2_20rowavgimpute))
C3KNNavgimputed20<- cbind.data.frame(rowMeans(CONTROL3_20KNNimputed), rowMeans(control3_20rowavgimpute))

write.csv(C1KNNavgimputed20 , "C1KNNavgimputed20.csv")
write.csv(C2KNNavgimputed20 , "C2KNNavgimputed20.csv")
write.csv(C3KNNavgimputed20 , "C3KNNavgimputed20.csv")

plot(x = rowMeans(CONTROL1_20KNNimputed) ,y = rowMeans(control1_20rowavgimpute),
     xlab = "kNN",
     ylab = "Average impute",
     col= "#E15759",
     type='p', cex= 2, lwd= 3,
#     xlim = c(2.5,5),
 #    ylim = c(15,30),		 
     main = "kNN vs Average Impute control1"
)
points(x = rowMeans(CONTROL2_20KNNimputed) ,y = rowMeans(control2_20rowavgimpute),
       col = "#F28E2B",
       type='p', cex= 2, lwd= 3)
points(x = rowMeans(CONTROL3_20KNNimputed) ,y = rowMeans(control3_20rowavgimpute),
       col = "#4E79A7",
       type='p', cex= 2, lwd= 3)

#########################################################
######Extracted protein matrix wrt 20%imputed kNN biotin
#########################################################

Imputed20biotinc1 <- merge(CONTROL1_20KNNimputed, biotin_20KNNimputed,
                          by = 'row.names', all = TRUE)

row.names(Imputed20biotinc1)<- Imputed20biotinc1$Row.names
Imputed20biotinc1$Row.names<- NULL
Imputed20biotinc1<- na.omit(Imputed20biotinc1)

dim(Imputed20biotinc1)
head(Imputed20biotinc1)


Imputed20biotinc2 <- merge(CONTROL2_20KNNimputed, biotin_20KNNimputed,
                           by = 'row.names', all = TRUE)

row.names(Imputed20biotinc2)<- Imputed20biotinc2$Row.names
Imputed20biotinc2$Row.names<- NULL
Imputed20biotinc2<- na.omit(Imputed20biotinc2)
dim(Imputed20biotinc2)
head(Imputed20biotinc2)


Imputed20biotinc3 <- merge(CONTROL3_20KNNimputed, biotin_20KNNimputed,
                           by = 'row.names', all = TRUE)

row.names(Imputed20biotinc3)<- Imputed20biotinc3$Row.names
Imputed20biotinc3$Row.names<- NULL
Imputed20biotinc3<- na.omit(Imputed20biotinc3)
dim(Imputed20biotinc3)
head(Imputed20biotinc3)



# FINAL MATRICES after 20% missing imputation in kNN
#NOW apply sample loading normalization and TMM and perform DE analysis
Imputed20biotinc1
Imputed20biotinc2
Imputed20biotinc3


#SL normalizationSL

Imputed20biotinc1_sl <- SL_norm(data.frame(Imputed20biotinc1))
Imputed20biotinc2_sl <- SL_norm(data.frame(Imputed20biotinc2))
Imputed20biotinc3_sl <- SL_norm(data.frame(Imputed20biotinc3))



##RAW vs SL in log 10 comparison
emf(file = "biotinC1_20kNN_SLboxplot.emf")
boxplot(c(log(Imputed20biotinc1,10),log(Imputed20biotinc1_sl,10)) ,
        col = c(rep("red", 6), rep("blue",6),rep("green", 6), rep("orange",6)),
        ylab  ="MS1 abundance in log 10" , 
        par(cex.lab=1.5), # is for y-axis
        par(cex.axis=1), # is for x-axis
        main = "Abundance before SL: red(biotin),blue(c1), after SL: green(biotin_SL),orange(cl1_SL)",
        names = c(rep("biotin", 6), rep("Control1", 6),rep("biotin_SL",6), rep("control1_SL",6)))
dev.off()

emf(file = "biotinC2_20kNN_SLboxplot.emf")
boxplot(c(log(Imputed20biotinc2,10),log(Imputed20biotinc2_sl,10)) ,
        col = c(rep("red", 6), rep("blue",3),rep("green", 6), rep("orange",3)),
        ylab  ="MS1 abundance in log 10" , 
        par(cex.lab=1.5), # is for y-axis
        par(cex.axis=1), # is for x-axis
        main = "Abundance before SL: red(biotin),blue(c2), after SL: green(biotin_SL),orange(cl2_SL)",
        names = c(rep("biotin", 6), rep("Control2", 3),rep("biotin_SL",6), rep("control2_SL",3)))

dev.off()

emf(file = "biotinC3_20kNN_SLboxplot.emf")
boxplot(c(log(Imputed20biotinc3,10),log(Imputed20biotinc3_sl,10)) ,
        col = c(rep("red", 6), rep("blue",3),rep("green", 6), rep("orange",3)),
        ylab  ="MS1 abundance in log 10" , 
        par(cex.lab=1.5), # is for y-axis
        par(cex.axis=1), # is for x-axis
        main = "Abundance before SL: red(biotin),blue(c3), after SL: green(biotin_SL),orange(cl3_SL)",
        names = c(rep("biotin", 6), rep("Control3", 3),rep("biotin_SL",6), rep("control3_SL",3)))
dev.off()


format(round(colSums(Imputed20biotinc1_sl), digits = 0), big.mark = ",") #ALL sum sould be same
format(round(colSums(Imputed20biotinc2_sl), digits = 0), big.mark = ",")
format(round(colSums(Imputed20biotinc3_sl), digits = 0), big.mark = ",")

#TMM NORMLAIZATION

sl_tmm <- calcNormFactors(Imputed20biotinc1_sl)
Imputed20biotinc1_sl_tmm <- sweep(Imputed20biotinc1_sl, 2, sl_tmm, FUN = "/") # data after SL and TMM on original scale

sl_tmm <- calcNormFactors(Imputed20biotinc2_sl)
Imputed20biotinc2_sl_tmm <- sweep(Imputed20biotinc2_sl, 2, sl_tmm, FUN = "/") # data after SL and TMM on original scale

sl_tmm <- calcNormFactors(Imputed20biotinc3_sl)
Imputed20biotinc3_sl_tmm <- sweep(Imputed20biotinc3_sl, 2, sl_tmm, FUN = "/") # data after SL and TMM on original scale


write.csv(Imputed20biotinc1_sl_tmm, "Imputed20kNNbiotinc1_sl_tmm.csv")
write.csv(Imputed20biotinc2_sl_tmm, "Imputed20kNNbiotinc2_sl_tmm.csv")
write.csv(Imputed20biotinc3_sl_tmm, "Imputed20kNNbiotinc3_sl_tmm.csv")

##RAW vs SL & TMM in log 10 comparison

emf(file = "biotinC1_20kNN_SL_TMMboxplot.emf")
boxplot(c(log(Imputed20biotinc1,10),log(Imputed20biotinc1_sl_tmm,10)) ,
        col = c(rep("red", 6), rep("blue",6),rep("green", 6), rep("orange",6)),
        ylab  ="MS1 abundance in log 10" , 
        par(cex.lab=1.5), # is for y-axis
        par(cex.axis=1), # is for x-axis
        main = "20%kNN red(biotin),blue(c1), SL & TMM: green(biotin),orange(c1)",
        names = c(rep("biotin", 6), rep("Control1", 6),rep("biotin_SL_TMM",6), rep("control1_SL_TMM",6)))

dev.off()

emf(file = "biotinC2_20kNN_SL_TMMboxplot.emf")
boxplot(c(log(Imputed20biotinc2,10),log(Imputed20biotinc2_sl_tmm,10)) ,
        col = c(rep("red", 6), rep("blue",3),rep("green", 6), rep("orange",3)),
        ylab  ="MS1 abundance in log 10" , 
        par(cex.lab=1.5), # is for y-axis
        par(cex.axis=1), # is for x-axis
        main = "20%kNN red(biotin),blue(c2), SL & TMM: green(biotin),orange(c2)",
        names = c(rep("biotin", 6), rep("Control2", 3),rep("biotin_SL_TMM",6), rep("control2_SL_TMM",3)))

dev.off()

emf(file = "biotinC3_20kNN_SL_TMMboxplot.emf")
boxplot(c(log(Imputed20biotinc3,10),log(Imputed20biotinc3_sl_tmm,10)) ,
        col = c(rep("red", 6), rep("blue",3),rep("green", 6), rep("orange",3)),
        ylab  ="MS1 abundance in log 10" , 
        par(cex.lab=1.5), # is for y-axis
        par(cex.axis=1), # is for x-axis
        main = "20%kNN red(biotin),blue(c3), SL & TMM: green(biotin),orange(c3)",
        names = c(rep("biotin", 6), rep("Control3", 3),rep("biotin_SL_TMM",6), rep("control3_SL_TMM",3)))
dev.off()


dim(Imputed20biotinc1_sl_tmm)
dim(Imputed20biotinc2_sl_tmm)
dim(Imputed20biotinc3_sl_tmm)

Imputed20biotinc1_sl_tmm_biotin<- rowMeans(Imputed20biotinc1_sl_tmm[,1:6], na.rm = TRUE)
Imputed20biotinc1_sl_tmm_c1<- rowMeans(Imputed20biotinc1_sl_tmm[,7:12], na.rm = TRUE)

Imputed20biotinc2_sl_tmm_biotin<- rowMeans(Imputed20biotinc2_sl_tmm[,1:6], na.rm = TRUE)
Imputed20biotinc2_sl_tmm_c2<- rowMeans(Imputed20biotinc2_sl_tmm[,7:9], na.rm = TRUE)

Imputed20biotinc3_sl_tmm_biotin<- rowMeans(Imputed20biotinc3_sl_tmm[,1:6], na.rm = TRUE)
Imputed20biotinc3_sl_tmm_c3<- rowMeans(Imputed20biotinc3_sl_tmm[,7:9], na.rm = TRUE)

LFImputed20biotinc1_sl_tmm = log(Imputed20biotinc1_sl_tmm_biotin/Imputed20biotinc1_sl_tmm_c1, 2)
LFImputed20biotinc2_sl_tmm = log(Imputed20biotinc2_sl_tmm_biotin/Imputed20biotinc2_sl_tmm_c2, 2)
LFImputed20biotinc3_sl_tmm = log(Imputed20biotinc3_sl_tmm_biotin/Imputed20biotinc3_sl_tmm_c3, 2)





# log(10/2,2)
# log(10,2)-log(2,2)
#####################################
##DE proteins/genes biotin vs Control1 
#####################################

condition<- "control"
trimdata0_CON <- Imputed20biotinc1_sl_tmm[ ,grepl(condition , colnames(Imputed20biotinc1_sl_tmm) ) ]
dim(trimdata0_CON)

condition<- "biotin"
trimdata0_biotin <- Imputed20biotinc1_sl_tmm[ ,grepl(condition , colnames(Imputed20biotinc1_sl_tmm) ) ]
dim(trimdata0_biotin)

trimdata0convsbiotin<- cbind.data.frame(trimdata0_CON, trimdata0_biotin)



# CRvsC<- cbind.data.frame(UMHET3MCR, UMHET3MC)
# head (CRvsC)
# dim(CRvsC)


# dim(UMHET3MCR)
# dim(UMHET3MC)
# group = factor(rep(c("CR","C"),c(9,10)))

group_trimdata0 = factor(rep(c("CON","biotin"),c(6,6)))


# res = rowttests(as.matrix(CRvsC),factor(group))
# res$adjP = p.adjust(res$p.value,"BH")
# head(res)



######Assuming there is difference in variance in two groups

library(broom)
res_trimdata0 = apply(trimdata0convsbiotin,1,function(i)tidy(t.test(i ~ group_trimdata0)))



res_trimdata0 = do.call(rbind,res_trimdata0)
res_trimdata0$adjP = p.adjust(res_trimdata0$p.value,"BH")




head(res_trimdata0)
#res1[,c("statistic","p.value","adjP")]

#row.names(res1)<- rownames(FCR)
row.names(res_trimdata0)<- row.names(trimdata0convsbiotin)



# head(res1)
# 
# write.csv(res1,"welch2samplettestGHRCOvsCfemalec57bl6.csv", row.names = TRUE, col.names = TRUE)


#####################################
##DE proteins/genes CR vs C VOLVANO PLOT
#####################################
#DE_CRvsC<- cbind.data.frame(rownames(CRvsC),CRlog2FCratio, res$adjP)
DEtrimdata0df<-cbind.data.frame(rownames(trimdata0convsbiotin),LFImputed20biotinc1_sl_tmm ,res_trimdata0$adjP)


colnames(DEtrimdata0df)<- c("protein", "log2FC", "AdjPvalue")


write.csv(DEtrimdata0df,"DEkNNimputed20biotinvsc1_sl_tmm.csv")
-log(0.05,10)



# 
# if (!requireNamespace('BiocManager', quietly = TRUE))
#   install.packages('BiocManager')
# 
# BiocManager::install('EnhancedVolcano')
library(EnhancedVolcano)
#

#
# EnhancedVolcano(DEtrimdata0df,
#                 lab = DEtrimdata0df$protein ,
#                 x = 'log2FC',
#                 y = 'AdjPvalue')

proteingene<- DEtrimdata0df$protein

DEtrimdata0df$protein <- gsub("sp.", "", DEtrimdata0df$protein)
DEtrimdata0df$protein <- gsub("_HUMAN", "", DEtrimdata0df$protein)



VOLCANO<- function(DE,conditionname){
  EnhancedVolcano(DE,
                  lab = DE$protein,
                  x = 'log2FC',
                  y = 'AdjPvalue',
                  title = paste0(conditionname," versus Control"),
                  pCutoff = .05,
                  FCcutoff = 2,
                  pointSize = 5,
                  labSize = 6,
                  ylim = c(0, 7)
                  #xlim = c(-1, 1)
  )
}

emf(file = "biotinC1_20kNN_SL_TMMVplot.emf")
VOLCANO(DEtrimdata0df,"biotin-C1 20%kNN SL & TMM LogFC >2 and adj.pval <.05")
dev.off()

write.csv (DEtrimdata0df,"DEbiotinC120kNNimputeSL_TMM.csv")


##################################################

#####################################
##DE proteins/genes biotin vs Control2 
#####################################

condition<- "hekwt_b"
trimdata0_CON <- Imputed20biotinc2_sl_tmm[ ,grepl(condition , colnames(Imputed20biotinc2_sl_tmm) ) ]
dim(trimdata0_CON)

condition<- "biotin"
trimdata0_biotin <- Imputed20biotinc2_sl_tmm[ ,grepl(condition , colnames(Imputed20biotinc2_sl_tmm) ) ]
dim(trimdata0_biotin)

trimdata0convsbiotin<- cbind.data.frame(trimdata0_CON, trimdata0_biotin)



# CRvsC<- cbind.data.frame(UMHET3MCR, UMHET3MC)
# head (CRvsC)
# dim(CRvsC)


# dim(UMHET3MCR)
# dim(UMHET3MC)
# group = factor(rep(c("CR","C"),c(9,10)))

group_trimdata0 = factor(rep(c("CON","biotin"),c(3,6)))


# res = rowttests(as.matrix(CRvsC),factor(group))
# res$adjP = p.adjust(res$p.value,"BH")
# head(res)



######Assuming there is difference in variance in two groups

library(broom)
res_trimdata0 = apply(trimdata0convsbiotin,1,function(i)tidy(t.test(i ~ group_trimdata0)))



res_trimdata0 = do.call(rbind,res_trimdata0)
res_trimdata0$adjP = p.adjust(res_trimdata0$p.value,"BH")




head(res_trimdata0)
#res1[,c("statistic","p.value","adjP")]

#row.names(res1)<- rownames(FCR)
row.names(res_trimdata0)<- row.names(trimdata0convsbiotin)



# head(res1)
# 
# write.csv(res1,"welch2samplettestGHRCOvsCfemalec57bl6.csv", row.names = TRUE, col.names = TRUE)


#####################################
##DE proteins/genes CR vs C VOLVANO PLOT
#####################################
#DE_CRvsC<- cbind.data.frame(rownames(CRvsC),CRlog2FCratio, res$adjP)
DEtrimdata0df<-cbind.data.frame(rownames(trimdata0convsbiotin),LFImputed20biotinc2_sl_tmm ,res_trimdata0$adjP)


colnames(DEtrimdata0df)<- c("protein", "log2FC", "AdjPvalue")


write.csv(DEtrimdata0df,"DEkNNimputed20biotinvsc2_sl_tmm.csv")
-log(0.05,10)



# 
# if (!requireNamespace('BiocManager', quietly = TRUE))
#   install.packages('BiocManager')
# 
# BiocManager::install('EnhancedVolcano')
library(EnhancedVolcano)
#

#
# EnhancedVolcano(DEtrimdata0df,
#                 lab = DEtrimdata0df$protein ,
#                 x = 'log2FC',
#                 y = 'AdjPvalue')

proteingene<- DEtrimdata0df$protein

DEtrimdata0df$protein <- gsub("sp.", "", DEtrimdata0df$protein)
DEtrimdata0df$protein <- gsub("_HUMAN", "", DEtrimdata0df$protein)



VOLCANO<- function(DE,conditionname){
  EnhancedVolcano(DE,
                  lab = DE$protein,
                  x = 'log2FC',
                  y = 'AdjPvalue',
                  title = paste0(conditionname," versus Control"),
                  pCutoff = .05,
                  FCcutoff = 2,
                  pointSize = 5,
                  labSize = 6,
                  ylim = c(0, 7)
                  #xlim = c(-1, 1)
  )
}

emf(file = "biotinC2_20kNN_SL_TMMVplot.emf")
VOLCANO(DEtrimdata0df,"biotin-C2 20%kNN SL & TMM LogFC >2 and adj.pval <.05")
dev.off()


write.csv (DEtrimdata0df,"DEbiotinC220kNNimputeSL_TMM.csv")


##################################################

#####################################
##DE proteins/genes biotin vs Control3 
#####################################

condition<- "hekwt_c"
trimdata0_CON <- Imputed20biotinc3_sl_tmm[ ,grepl(condition , colnames(Imputed20biotinc3_sl_tmm) ) ]
dim(trimdata0_CON)

condition<- "biotin"
trimdata0_biotin <- Imputed20biotinc3_sl_tmm[ ,grepl(condition , colnames(Imputed20biotinc3_sl_tmm) ) ]
dim(trimdata0_biotin)

trimdata0convsbiotin<- cbind.data.frame(trimdata0_CON, trimdata0_biotin)



# CRvsC<- cbind.data.frame(UMHET3MCR, UMHET3MC)
# head (CRvsC)
# dim(CRvsC)


# dim(UMHET3MCR)
# dim(UMHET3MC)
# group = factor(rep(c("CR","C"),c(9,10)))

group_trimdata0 = factor(rep(c("CON","biotin"),c(3,6)))


# res = rowttests(as.matrix(CRvsC),factor(group))
# res$adjP = p.adjust(res$p.value,"BH")
# head(res)



######Assuming there is difference in variance in two groups

library(broom)
res_trimdata0 = apply(trimdata0convsbiotin,1,function(i)tidy(t.test(i ~ group_trimdata0)))



res_trimdata0 = do.call(rbind,res_trimdata0)
res_trimdata0$adjP = p.adjust(res_trimdata0$p.value,"BH")




head(res_trimdata0)
#res1[,c("statistic","p.value","adjP")]

#row.names(res1)<- rownames(FCR)
row.names(res_trimdata0)<- row.names(trimdata0convsbiotin)



# head(res1)
# 
# write.csv(res1,"welch2samplettestGHRCOvsCfemalec57bl6.csv", row.names = TRUE, col.names = TRUE)


#####################################
##DE proteins/genes CR vs C VOLVANO PLOT
#####################################
#DE_CRvsC<- cbind.data.frame(rownames(CRvsC),CRlog2FCratio, res$adjP)
DEtrimdata0df<-cbind.data.frame(rownames(trimdata0convsbiotin),LFImputed20biotinc3_sl_tmm ,res_trimdata0$adjP)


colnames(DEtrimdata0df)<- c("protein", "log2FC", "AdjPvalue")


write.csv(DEtrimdata0df,"DEkNNimputed20biotinvsc3_sl_tmm.csv")
-log(0.05,10)



# 
# if (!requireNamespace('BiocManager', quietly = TRUE))
#   install.packages('BiocManager')
# 
# BiocManager::install('EnhancedVolcano')
library(EnhancedVolcano)
#

#
# EnhancedVolcano(DEtrimdata0df,
#                 lab = DEtrimdata0df$protein ,
#                 x = 'log2FC',
#                 y = 'AdjPvalue')

proteingene<- DEtrimdata0df$protein

DEtrimdata0df$protein <- gsub("sp.", "", DEtrimdata0df$protein)
DEtrimdata0df$protein <- gsub("_HUMAN", "", DEtrimdata0df$protein)



VOLCANO<- function(DE,conditionname){
  EnhancedVolcano(DE,
                  lab = DE$protein,
                  x = 'log2FC',
                  y = 'AdjPvalue',
                  title = paste0(conditionname," versus Control"),
                  pCutoff = .05,
                  FCcutoff = 2,
                  pointSize = 5,
                  labSize = 6,
                  ylim = c(0, 7)
                  #xlim = c(-1, 1)
  )
}

emf(file = "biotinC3_20kNN_SL_TMMVplot.emf")
VOLCANO(DEtrimdata0df,"biotin-C3 20%kNN SL & TMM LogFC >2 and adj.pval <.05")
dev.off()


write.csv (DEtrimdata0df,"DEbiotinC320kNNimputeSL_TMM.csv")









#########################################################
######Extracted protein matrix wrt 20% imputed Avg biotin
#########################################################

Imputed20biotinc1 <- merge(control1_20rowavgimpute, biotin_20KNNimputed,
                           by = 'row.names', all = TRUE)

row.names(Imputed20biotinc1)<- Imputed20biotinc1$Row.names
Imputed20biotinc1$Row.names<- NULL
Imputed20biotinc1<- na.omit(Imputed20biotinc1)
dim(Imputed20biotinc1)
head(Imputed20biotinc1)


Imputed20biotinc2 <- merge(control2_20rowavgimpute, biotin_20KNNimputed,
                           by = 'row.names', all = TRUE)

row.names(Imputed20biotinc2)<- Imputed20biotinc2$Row.names
Imputed20biotinc2$Row.names<- NULL
Imputed20biotinc2<- na.omit(Imputed20biotinc2)
dim(Imputed20biotinc2)
head(Imputed20biotinc2)


Imputed20biotinc3 <- merge(control3_20rowavgimpute, biotin_20KNNimputed,
                           by = 'row.names', all = TRUE)

row.names(Imputed20biotinc3)<- Imputed20biotinc3$Row.names
Imputed20biotinc3$Row.names<- NULL
Imputed20biotinc3<- na.omit(Imputed20biotinc3)
dim(Imputed20biotinc3)
head(Imputed20biotinc3)



# FINAL MATRICES after 20% missing imputation in kNN
#NOW apply sample loading normalization and TMM and perform DE analysis
Imputed20biotinc1
Imputed20biotinc2
Imputed20biotinc3


#SL normalizationSL

Imputed20biotinc1_sl <- SL_norm(data.frame(Imputed20biotinc1))
Imputed20biotinc2_sl <- SL_norm(data.frame(Imputed20biotinc2))
Imputed20biotinc3_sl <- SL_norm(data.frame(Imputed20biotinc3))



##RAW vs SL in log 10 comparison
boxplot(c(log(Imputed20biotinc1,10),log(Imputed20biotinc1_sl,10)) ,
        col = c(rep("red", 6), rep("blue",6),rep("green", 6), rep("orange",6)),
        ylab  ="MS1 abundance in log 10" , 
        par(cex.lab=1.5), # is for y-axis
        par(cex.axis=1), # is for x-axis
        main = "Abundance before SL: red(biotin),blue(c1), after SL: green(biotin_SL),orange(cl1_SL)",
        names = c(rep("biotin", 6), rep("Control1", 6),rep("biotin_SL",6), rep("control1_SL",6)))

boxplot(c(log(Imputed20biotinc2,10),log(Imputed20biotinc2_sl,10)) ,
        col = c(rep("red", 6), rep("blue",3),rep("green", 6), rep("orange",3)),
        ylab  ="MS1 abundance in log 10" , 
        par(cex.lab=1.5), # is for y-axis
        par(cex.axis=1), # is for x-axis
        main = "Abundance before SL: red(biotin),blue(c2), after SL: green(biotin_SL),orange(cl2_SL)",
        names = c(rep("biotin", 6), rep("Control2", 3),rep("biotin_SL",6), rep("control2_SL",3)))

boxplot(c(log(Imputed20biotinc3,10),log(Imputed20biotinc3_sl,10)) ,
        col = c(rep("red", 6), rep("blue",3),rep("green", 6), rep("orange",3)),
        ylab  ="MS1 abundance in log 10" , 
        par(cex.lab=1.5), # is for y-axis
        par(cex.axis=1), # is for x-axis
        main = "Abundance before SL: red(biotin),blue(c3), after SL: green(biotin_SL),orange(cl3_SL)",
        names = c(rep("biotin", 6), rep("Control3", 3),rep("biotin_SL",6), rep("control3_SL",3)))


format(round(colSums(Imputed20biotinc1_sl), digits = 0), big.mark = ",") #ALL sum sould be same
format(round(colSums(Imputed20biotinc2_sl), digits = 0), big.mark = ",")
format(round(colSums(Imputed20biotinc3_sl), digits = 0), big.mark = ",")

#TMM NORMLAIZATION

sl_tmm <- calcNormFactors(Imputed20biotinc1_sl)
Imputed20biotinc1_sl_tmm <- sweep(Imputed20biotinc1_sl, 2, sl_tmm, FUN = "/") # data after SL and TMM on original scale

sl_tmm <- calcNormFactors(Imputed20biotinc2_sl)
Imputed20biotinc2_sl_tmm <- sweep(Imputed20biotinc2_sl, 2, sl_tmm, FUN = "/") # data after SL and TMM on original scale

sl_tmm <- calcNormFactors(Imputed20biotinc3_sl)
Imputed20biotinc3_sl_tmm <- sweep(Imputed20biotinc3_sl, 2, sl_tmm, FUN = "/") # data after SL and TMM on original scale

write.csv(Imputed20biotinc1_sl_tmm, "Imputed20Avgbiotinc1_sl_tmm.csv")
write.csv(Imputed20biotinc2_sl_tmm, "Imputed20Avgbiotinc2_sl_tmm.csv")
write.csv(Imputed20biotinc3_sl_tmm, "Imputed20Avgbiotinc3_sl_tmm.csv")

##RAW vs SL & TMM in log 10 comparison
emf(file= "biotinC1_20Avgimute_SL_TMMboxplot.emf")
boxplot(c(log(Imputed20biotinc1,10),log(Imputed20biotinc1_sl_tmm,10)) ,
        col = c(rep("red", 6), rep("blue",6),rep("green", 6), rep("orange",6)),
        ylab  ="MS1 abundance in log 10" , 
        par(cex.lab=1.5), # is for y-axis
        par(cex.axis=1), # is for x-axis
        main = "20% Avg impute: red(biotin),blue(c1), SL & TMM: green(biotin),orange(c1)",
        names = c(rep("biotin", 6), rep("Control1", 6),rep("biotin_SL_TMM",6), rep("control1_SL_TMM",6)))
dev.off()

emf(file= "biotinC2_20Avgimute_SL_TMMboxplot.emf")
boxplot(c(log(Imputed20biotinc2,10),log(Imputed20biotinc2_sl_tmm,10)) ,
        col = c(rep("red", 6), rep("blue",3),rep("green", 6), rep("orange",3)),
        ylab  ="MS1 abundance in log 10" , 
        par(cex.lab=1.5), # is for y-axis
        par(cex.axis=1), # is for x-axis
        main = "20% Avg impute: red(biotin),blue(c2), SL & TMM: green(biotin),orange(c2)",
        names = c(rep("biotin", 6), rep("Control2", 3),rep("biotin_SL_TMM",6), rep("control2_SL_TMM",3)))
dev.off()

emf(file= "biotinC3_20Avgimute_SL_TMMboxplot.emf")
boxplot(c(log(Imputed20biotinc3,10),log(Imputed20biotinc3_sl_tmm,10)) ,
        col = c(rep("red", 6), rep("blue",3),rep("green", 6), rep("orange",3)),
        ylab  ="MS1 abundance in log 10" , 
        par(cex.lab=1.5), # is for y-axis
        par(cex.axis=1), # is for x-axis
        main = "20% Avg impute: red(biotin),blue(c3), SL & TMM: green(biotin),orange(c3)",
        names = c(rep("biotin", 6), rep("Control3", 3),rep("biotin_SL_TMM",6), rep("control3_SL_TMM",3)))
dev.off()


dim(Imputed20biotinc1_sl_tmm)
dim(Imputed20biotinc2_sl_tmm)
dim(Imputed20biotinc3_sl_tmm)

Imputed20biotinc1_sl_tmm_biotin<- rowMeans(Imputed20biotinc1_sl_tmm[,1:6], na.rm = TRUE)
Imputed20biotinc1_sl_tmm_c1<- rowMeans(Imputed20biotinc1_sl_tmm[,7:12], na.rm = TRUE)

Imputed20biotinc2_sl_tmm_biotin<- rowMeans(Imputed20biotinc2_sl_tmm[,1:6], na.rm = TRUE)
Imputed20biotinc2_sl_tmm_c2<- rowMeans(Imputed20biotinc2_sl_tmm[,7:9], na.rm = TRUE)

Imputed20biotinc3_sl_tmm_biotin<- rowMeans(Imputed20biotinc3_sl_tmm[,1:6], na.rm = TRUE)
Imputed20biotinc3_sl_tmm_c3<- rowMeans(Imputed20biotinc3_sl_tmm[,7:9], na.rm = TRUE)

LFImputed20biotinc1_sl_tmm = log(Imputed20biotinc1_sl_tmm_biotin/Imputed20biotinc1_sl_tmm_c1, 2)
LFImputed20biotinc2_sl_tmm = log(Imputed20biotinc2_sl_tmm_biotin/Imputed20biotinc2_sl_tmm_c2, 2)
LFImputed20biotinc3_sl_tmm = log(Imputed20biotinc3_sl_tmm_biotin/Imputed20biotinc3_sl_tmm_c3, 2)





# log(10/2,2)
# log(10,2)-log(2,2)
#####################################
##DE proteins/genes biotin vs Control1 
#####################################

condition<- "control"
trimdata0_CON <- Imputed20biotinc1_sl_tmm[ ,grepl(condition , colnames(Imputed20biotinc1_sl_tmm) ) ]
dim(trimdata0_CON)

condition<- "biotin"
trimdata0_biotin <- Imputed20biotinc1_sl_tmm[ ,grepl(condition , colnames(Imputed20biotinc1_sl_tmm) ) ]
dim(trimdata0_biotin)

trimdata0convsbiotin<- cbind.data.frame(trimdata0_CON, trimdata0_biotin)



# CRvsC<- cbind.data.frame(UMHET3MCR, UMHET3MC)
# head (CRvsC)
# dim(CRvsC)


# dim(UMHET3MCR)
# dim(UMHET3MC)
# group = factor(rep(c("CR","C"),c(9,10)))

group_trimdata0 = factor(rep(c("CON","biotin"),c(6,6)))


# res = rowttests(as.matrix(CRvsC),factor(group))
# res$adjP = p.adjust(res$p.value,"BH")
# head(res)



######Assuming there is difference in variance in two groups

library(broom)
res_trimdata0 = apply(trimdata0convsbiotin,1,function(i)tidy(t.test(i ~ group_trimdata0)))



res_trimdata0 = do.call(rbind,res_trimdata0)
res_trimdata0$adjP = p.adjust(res_trimdata0$p.value,"BH")




head(res_trimdata0)
#res1[,c("statistic","p.value","adjP")]

#row.names(res1)<- rownames(FCR)
row.names(res_trimdata0)<- row.names(trimdata0convsbiotin)



# head(res1)
# 
# write.csv(res1,"welch2samplettestGHRCOvsCfemalec57bl6.csv", row.names = TRUE, col.names = TRUE)


#####################################
##DE proteins/genes CR vs C VOLVANO PLOT
#####################################
#DE_CRvsC<- cbind.data.frame(rownames(CRvsC),CRlog2FCratio, res$adjP)
DEtrimdata0df<-cbind.data.frame(rownames(trimdata0convsbiotin),LFImputed20biotinc1_sl_tmm ,res_trimdata0$adjP)


colnames(DEtrimdata0df)<- c("protein", "log2FC", "AdjPvalue")


write.csv(DEtrimdata0df,"DEAVGimputed20biotinvsc1_sl_tmm.csv")
-log(0.05,10)



# 
# if (!requireNamespace('BiocManager', quietly = TRUE))
#   install.packages('BiocManager')
# 
# BiocManager::install('EnhancedVolcano')
library(EnhancedVolcano)
#

#
# EnhancedVolcano(DEtrimdata0df,
#                 lab = DEtrimdata0df$protein ,
#                 x = 'log2FC',
#                 y = 'AdjPvalue')

proteingene<- DEtrimdata0df$protein

DEtrimdata0df$protein <- gsub("sp.", "", DEtrimdata0df$protein)
DEtrimdata0df$protein <- gsub("_HUMAN", "", DEtrimdata0df$protein)



VOLCANO<- function(DE,conditionname){
  EnhancedVolcano(DE,
                  lab = DE$protein,
                  x = 'log2FC',
                  y = 'AdjPvalue',
                  title = paste0(conditionname," versus Control"),
                  pCutoff = .05,
                  FCcutoff = 2,
                  pointSize = 5,
                  labSize = 6,
                  ylim = c(0, 7)
                  #xlim = c(-1, 1)
  )
}
emf(file= "biotinC1_20Avgimute_SL_TMMVplot.emf")

VOLCANO(DEtrimdata0df,"biotin-C1 Avg imputed 20% SL & TMM LogFC >2, adj.pval <.05")
dev.off()

write.csv (DEtrimdata0df,"DEbiotinC120AvgimputeSL_TMM.csv")


##################################################

#####################################
##DE proteins/genes biotin vs Control2 
#####################################


condition<- "hekwt_b"
trimdata0_CON <- Imputed20biotinc2_sl_tmm[ ,grepl(condition , colnames(Imputed20biotinc2_sl_tmm) ) ]
dim(trimdata0_CON)

condition<- "biotin"
trimdata0_biotin <- Imputed20biotinc2_sl_tmm[ ,grepl(condition , colnames(Imputed20biotinc2_sl_tmm) ) ]
dim(trimdata0_biotin)

trimdata0convsbiotin<- cbind.data.frame(trimdata0_CON, trimdata0_biotin)



# CRvsC<- cbind.data.frame(UMHET3MCR, UMHET3MC)
# head (CRvsC)
# dim(CRvsC)


# dim(UMHET3MCR)
# dim(UMHET3MC)
# group = factor(rep(c("CR","C"),c(9,10)))

group_trimdata0 = factor(rep(c("CON","biotin"),c(3,6)))


# res = rowttests(as.matrix(CRvsC),factor(group))
# res$adjP = p.adjust(res$p.value,"BH")
# head(res)



######Assuming there is difference in variance in two groups

library(broom)
res_trimdata0 = apply(trimdata0convsbiotin,1,function(i)tidy(t.test(i ~ group_trimdata0)))



res_trimdata0 = do.call(rbind,res_trimdata0)
res_trimdata0$adjP = p.adjust(res_trimdata0$p.value,"BH")




head(res_trimdata0)
#res1[,c("statistic","p.value","adjP")]

#row.names(res1)<- rownames(FCR)
row.names(res_trimdata0)<- row.names(trimdata0convsbiotin)



# head(res1)
# 
# write.csv(res1,"welch2samplettestGHRCOvsCfemalec57bl6.csv", row.names = TRUE, col.names = TRUE)


#####################################
##DE proteins/genes CR vs C VOLVANO PLOT
#####################################
#DE_CRvsC<- cbind.data.frame(rownames(CRvsC),CRlog2FCratio, res$adjP)
DEtrimdata0df<-cbind.data.frame(rownames(trimdata0convsbiotin),LFImputed20biotinc2_sl_tmm ,res_trimdata0$adjP)


colnames(DEtrimdata0df)<- c("protein", "log2FC", "AdjPvalue")


write.csv(DEtrimdata0df,"DEAvgimputed20biotinvsc2_sl_tmm.csv")
-log(0.05,10)



# 
# if (!requireNamespace('BiocManager', quietly = TRUE))
#   install.packages('BiocManager')
# 
# BiocManager::install('EnhancedVolcano')
library(EnhancedVolcano)
#

#
# EnhancedVolcano(DEtrimdata0df,
#                 lab = DEtrimdata0df$protein ,
#                 x = 'log2FC',
#                 y = 'AdjPvalue')

proteingene<- DEtrimdata0df$protein

DEtrimdata0df$protein <- gsub("sp.", "", DEtrimdata0df$protein)
DEtrimdata0df$protein <- gsub("_HUMAN", "", DEtrimdata0df$protein)




VOLCANO<- function(DE,conditionname){
  EnhancedVolcano(DE,
                  lab = DE$protein,
                  x = 'log2FC',
                  y = 'AdjPvalue',
                  title = paste0(conditionname," versus Control"),
                  pCutoff = .05,
                  FCcutoff = 2,
                  pointSize = 5,
                  labSize = 6,
                  ylim = c(0, 7)
                  #xlim = c(-1, 1)
  )
}
emf(file= "biotinC2_20Avgimute_SL_TMMVplot.emf")

VOLCANO(DEtrimdata0df,"biotin-C2 Avg imputed 20% SL & TMM LogFC >2, adj.pval <.05")
dev.off()


write.csv (DEtrimdata0df,"DEbiotinC220AvgimputeSL_TMM.csv")


##################################################

#####################################
##DE proteins/genes biotin vs Control3 
#####################################

condition<- "hekwt_c"
trimdata0_CON <- Imputed20biotinc3_sl_tmm[ ,grepl(condition , colnames(Imputed20biotinc3_sl_tmm) ) ]
dim(trimdata0_CON)

condition<- "biotin"
trimdata0_biotin <- Imputed20biotinc3_sl_tmm[ ,grepl(condition , colnames(Imputed20biotinc3_sl_tmm) ) ]
dim(trimdata0_biotin)

trimdata0convsbiotin<- cbind.data.frame(trimdata0_CON, trimdata0_biotin)



# CRvsC<- cbind.data.frame(UMHET3MCR, UMHET3MC)
# head (CRvsC)
# dim(CRvsC)


# dim(UMHET3MCR)
# dim(UMHET3MC)
# group = factor(rep(c("CR","C"),c(9,10)))

group_trimdata0 = factor(rep(c("CON","biotin"),c(3,6)))


# res = rowttests(as.matrix(CRvsC),factor(group))
# res$adjP = p.adjust(res$p.value,"BH")
# head(res)



######Assuming there is difference in variance in two groups

library(broom)
res_trimdata0 = apply(trimdata0convsbiotin,1,function(i)tidy(t.test(i ~ group_trimdata0)))



res_trimdata0 = do.call(rbind,res_trimdata0)
res_trimdata0$adjP = p.adjust(res_trimdata0$p.value,"BH")




head(res_trimdata0)
#res1[,c("statistic","p.value","adjP")]

#row.names(res1)<- rownames(FCR)
row.names(res_trimdata0)<- row.names(trimdata0convsbiotin)



# head(res1)
# 
# write.csv(res1,"welch2samplettestGHRCOvsCfemalec57bl6.csv", row.names = TRUE, col.names = TRUE)


#####################################
##DE proteins/genes CR vs C VOLVANO PLOT
#####################################
#DE_CRvsC<- cbind.data.frame(rownames(CRvsC),CRlog2FCratio, res$adjP)
DEtrimdata0df<-cbind.data.frame(rownames(trimdata0convsbiotin),LFImputed20biotinc3_sl_tmm ,res_trimdata0$adjP)


colnames(DEtrimdata0df)<- c("protein", "log2FC", "AdjPvalue")


write.csv(DEtrimdata0df,"DEAvgimputed20biotinvsc3_sl_tmm.csv")
-log(0.05,10)



# 
# if (!requireNamespace('BiocManager', quietly = TRUE))
#   install.packages('BiocManager')
# 
# BiocManager::install('EnhancedVolcano')
library(EnhancedVolcano)
#

#
# EnhancedVolcano(DEtrimdata0df,
#                 lab = DEtrimdata0df$protein ,
#                 x = 'log2FC',
#                 y = 'AdjPvalue')

proteingene<- DEtrimdata0df$protein

DEtrimdata0df$protein <- gsub("sp.", "", DEtrimdata0df$protein)
DEtrimdata0df$protein <- gsub("_HUMAN", "", DEtrimdata0df$protein)




VOLCANO<- function(DE,conditionname){
  EnhancedVolcano(DE,
                  lab = DE$protein,
                  x = 'log2FC',
                  y = 'AdjPvalue',
                  title = paste0(conditionname," versus Control"),
                  pCutoff = .05,
                  FCcutoff = 2,
                  pointSize = 5,
                  labSize = 6,
                  ylim = c(0, 7)
                  #xlim = c(-1, 1)
  )
}
emf(file= "biotinC3_20Avgimute_SL_TMMVplot.emf")

VOLCANO(DEtrimdata0df,"biotin-C3 Avg imputed 20% SL & TMM LogFC >2, adj.pval <.05")
dev.off()

write.csv (DEtrimdata0df,"DEbiotinC320AvgimputeSL_TMM.csv")


###############EXTRACT exclusively expressed Biotin samples
EMPTYROWS<-function(data){
  data[data==0] <- NA
  data<-data[rowSums(is.na(data)) == ncol(data), ]
  #head (data)
  #dim(data)
  return(data)
}


dim(biotinvsc1_noemptyrow)
head(biotinvsc1_noemptyrow)
c1<- biotinvsc1_noemptyrow[7:12]
c1emptyrows<- EMPTYROWS(c1)
dim(c1emptyrows)

head(c1emptyrows)
Exclusivebiotinc1<- subset(biotinvsc1_noemptyrow[,1:6], rownames(biotinvsc1_noemptyrow) %in% rownames(c1emptyrows))
head(Exclusivebiotinc1)

write.csv(Exclusivebiotinc1, "Exclusivebiotinc1.csv")


dim(biotinvsc2_noemptyrow)
head(biotinvsc2_noemptyrow)
c2<- biotinvsc2_noemptyrow[7:9]
c2emptyrows<- EMPTYROWS(c2)
dim(c2emptyrows)

head(c2emptyrows)
Exclusivebiotinc2<- subset(biotinvsc2_noemptyrow[,1:6], rownames(biotinvsc2_noemptyrow) %in% rownames(c2emptyrows))
head(Exclusivebiotinc2)
write.csv(Exclusivebiotinc2, "Exclusivebiotinc2.csv")


dim(biotinvsc3_noemptyrow)
head(biotinvsc3_noemptyrow)
c3<- biotinvsc3_noemptyrow[7:9]
c3emptyrows<- EMPTYROWS(c3)
dim(c3emptyrows)

head(c3emptyrows)
Exclusivebiotinc3<- subset(biotinvsc3_noemptyrow[,1:6], rownames(biotinvsc3_noemptyrow) %in% rownames(c3emptyrows))
head(Exclusivebiotinc3)
write.csv(Exclusivebiotinc3, "Exclusivebiotinc3.csv")

###############EXTRACT non normalized data matrics and append to normalized dataset without imputation
dim(biotinvsc1_sl_tmm)
dim(biotinvsc1_noemptyrow)

head(biotinvsc1_sl_tmm)
head(biotinvsc1_noemptyrow)


Remaining_nonnormbiotinC1<- subset(biotinvsc1_noemptyrow, !(rownames(biotinvsc1_noemptyrow) %in% rownames(biotinvsc1_sl_tmm)))
head(Remaining_nonnormbiotinC1)
derivedSL_TMM_unnormBiotinC1<- rbind.data.frame(biotinvsc1_sl_tmm, Remaining_nonnormbiotinC1)
dim(derivedSL_TMM_unnormBiotinC1)


Remaining_nonnormbiotinC2<- subset(biotinvsc2_noemptyrow, !(rownames(biotinvsc2_noemptyrow) %in% rownames(biotinvsc2_sl_tmm)))
head(Remaining_nonnormbiotinC2)
derivedSL_TMM_unnormBiotinC2<- rbind.data.frame(biotinvsc2_sl_tmm, Remaining_nonnormbiotinC2)
dim(derivedSL_TMM_unnormBiotinC2)

Remaining_nonnormbiotinC3<- subset(biotinvsc3_noemptyrow, !(rownames(biotinvsc3_noemptyrow) %in% rownames(biotinvsc3_sl_tmm)))
head(Remaining_nonnormbiotinC3)
derivedSL_TMM_unnormBiotinC3<- rbind.data.frame(biotinvsc3_sl_tmm, Remaining_nonnormbiotinC3)
dim(derivedSL_TMM_unnormBiotinC3)


write.csv(derivedSL_TMM_unnormBiotinC1, "SL_TMM_nonormBiotinC1.csv")
write.csv(derivedSL_TMM_unnormBiotinC2, "SL_TMM_nonormBiotinC2.csv")
write.csv(derivedSL_TMM_unnormBiotinC3, "SL_TMM_nonormBiotinC3.csv")
