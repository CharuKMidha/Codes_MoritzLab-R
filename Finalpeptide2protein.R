########## finaL Peptide to protein -SUM SAME PEPTIDE WITH DIFF CHARGES, THEN TOP3 aggregation as average FOR PROTEIN ###########################
'/AUTHOR CHARU KAPIL MIDHA 11 May 2023
Script takes the library xpress generated matrics from zhi in peptide form 
<protein><peptide.charge><intensities>
- same peptide with different charge state is summed
- top 3 peptides are extracted and Averaged to quantify proteins
- protein matrix generated /'



rm(list = ls())
data_start <- read.csv("inputAll_Intensity2.csv")

#data_start <- read.csv("testraw.txt", sep = '\t')

head(data_start)


data_start[1:8,1:10]
data_start$Fragment<- NULL
dim(data_start)

data_start<- data_start[!grepl("CONTAM.*", data_start$Protein),]
dim(data_start)


##### GET PROTEOTYPIC
#proteotypic<- data_start[!grepl(";", data_start$protein),]
##### Remove DECOY from PROTEOTYPIC
# head(proteotypic)
# proteotypic[1:5,1:5]
# proteotypic<- proteotypic[!grepl("DECOY", proteotypic$protein),]
# proteotypic<- proteotypic[grepl("_HUMAN", proteotypic$protein),]
# dim(proteotypic)
# proteotypic$protein

pepprotid<- (paste(data_start$Peptide, data_start$Protein, sep = "#"))
#ncol(proteotypic)
#pepprot<- cbind(pepprotid, proteotypic[,4:ncol(proteotypic)] )
pepprot<- cbind(pepprotid, data_start )
#pepprot[ ,-1]
dim(pepprot)
colnames(pepprot)
#rownames(pepprot)<- pepprot$pepprotid
head(pepprot)
dim(pepprot)
pepprot[1:5,1:5]

pepprot$pepprotid<- gsub("..#", "#", pepprot$pepprotid)
length(rep(pepprot$pepprotid, (ncol(pepprot)-1)))
length(rep(names(pepprot[,-1]), nrow(pepprot)))
length(unlist(pepprot[ ,-1]))

#CONVERT MATRIX TO THREE ROWS
pepprot$Protein<- NULL
pepprot$Peptide<- NULL
result <- data.frame(protpepid = rep(pepprot$pepprotid, (ncol(pepprot)-1)),
                     name = rep(names(pepprot[,-1]), each = nrow(pepprot)), 
                     value = unlist(pepprot[ ,-1]), row.names = NULL)


result[1:5,]
pepprotsample<- (paste(result$protpepid, result$name, sep = "@"))

result$protpepid <- NULL
result$name <- NULL
#result
result2<- cbind (pepprotsample, result)

head(result2)
library(dplyr)
dim(result2)
#summary(result_nozero)
#result_nozero<- result2[-grep(0.000, result2$value), ]
#result_zero<- result2[grep(0.000, result2$value), ]
#head(result_nozero)

###REMOVED 0VALUES AND DECOYS

data.frame(result2)
result2[result2==0] <- NA
result_nozero<- na.omit(result2)


#result_nozero<- result2 %>% filter(!if_any(starts_with("value"), ~ . == NA))
dim(result_nozero)
#result_nozero<- result_nozero %>% filter(!if_any(starts_with("pepprotsample"), ~ . == "DECOY"))

#dim(result_nozero)
head(result_nozero)

#values<- c(pepprot[2:ncol(pepprot)])
require(data.table)
#dt <- data.table(pepprot)
dt <- data.table(result_nozero)
dt[1:5]
#typeof(dt$pepprotsample)
#as.numeric(dt$value)

#dt[dt == 0] <- NA

#########SUM OF same PEPTIDES with different charge state
#dt.mean <- dt[, lapply(.SD, mean), by = pepprotsample]
dim(dt)
dt.sum<- aggregate( value ~ pepprotsample, dt, sum)

dim(dt.sum)
head(dt.sum)

#WORKING Eg of keeping upto TOP 3 values

# data <- data.frame(group = c(rep(letters[1:3], each = 5), "d"),    # Create example data
#                    value = 1:16)
# data   
# 
# data_new1 <- data[order(data$value, decreasing = TRUE), ]
# data_new1 <- Reduce(rbind,                                 # Top N highest values by group
#                     by(data_new1,
#                        data_new1["group"],
#                        head,
#                        n = 3))




library(tidyr)
dtpepprotsample<- dt.sum %>%
  as.data.frame() %>%
  separate(pepprotsample, c("pep", "protsample"), "#")

dim(dtpepprotsample)
head(dtpepprotsample)
tail(dtpepprotsample)
dtpepprotsample$pep<- NULL
#dtpepprotsample<- data.table(dtpepprotsample)
tail(dt[order(dt$value, decreasing = TRUE), ])
data_new1 <- dtpepprotsample[order(dtpepprotsample$value, decreasing = TRUE), ]
tail(data_new1)


data_new1 <- Reduce(rbind,                                 # Top N highest values by group
                    by(data_new1,
                       data_new1["protsample"],
                       head,
                       n = 3))



dim(data_new1)
head(data_new1)


dt2.mean <- aggregate( value ~ protsample, data_new1, mean)
head(dt2.mean)

dtprotsample<- dt2.mean %>%
  as.data.frame() %>%
  separate(protsample, c("prot", "sample"), "@")

head(dtprotsample)

dtprotsample<- dtprotsample[, c("sample", "prot", "value")]

#library(tidyr)
head(dtprotsample)
datauniq <- unique(dtprotsample)
dim(dtprotsample)
dim(datauniq)

finaldata = datauniq %>%
  spread(key = prot,
         value = value)

head(finaldata)
dim(finaldata)

finaldata<- t(finaldata)
dim(finaldata)
rownames(finaldata)
finaldata[1:5,1:5]
colnames(finaldata)<-finaldata[1,];finaldata<-finaldata[2:nrow(finaldata),]
finaldata
protein<- rownames(finaldata)

finaldata<- cbind(protein, finaldata)
finaldata[1:5,1:5]

typeof(finaldata)
row.names(finaldata)<- NULL


finaldata[1:5,1:5]
length(colnames(finaldata))
#as.character(matrixdata)
dim(finaldata)
write.csv(finaldata, file = "proteinMEAN2304x19Osm28.csv", row.names = FALSE)

#Tinitialdata$rowsamples<- str_replace(Tinitialdata$rowsamples, "_200_","#200_")
#########SUMMATION OF PEPTIDES TO AGGREGATE TO PROTEINS
# dt.sum <- dt[, lapply(.SD, sum), by = pepprotsample]
# head(dt.sum)
# library(tidyr)
# dtpepprotsample<- dt.sum %>%
#   as.data.frame() %>%
#   separate(pepprotsample, c("pep", "protsample"), "#")
# 
# head(dtpepprotsample)
# dtpepprotsample$pep<- NULL
# dtpepprotsample<- data.table(dtpepprotsample)
# 
# dt2.sum <- dtpepprotsample[, lapply(.SD, sum), by = protsample]
# head(dt2.sum)
# 
# dtprotsample<- dt2.sum %>%
#   as.data.frame() %>%
#   separate(protsample, c("prot", "sample"), "@")
# 
# head(dtprotsample)
# 
# dtprotsample<- dtprotsample[, c("sample", "prot", "value")]
# 
# 
# 
# 
# 
# #library(tidyr)
# head(dtprotsample)
# datauniq <- unique(dtprotsample)
# dim(dtprotsample)
# dim(datauniq)
# 
# finaldata = datauniq %>%
#   spread(key = prot,
#          value = value)
# 
# head(finaldata)
# dim(finaldata)
# 
# finaldata<- t(finaldata)
# dim(finaldata)
# rownames(finaldata)
# finaldata[1:5,1:5]
# colnames(finaldata)<-finaldata[1,];finaldata<-finaldata[2:nrow(finaldata),]
# finaldata
# protein<- rownames(finaldata)
# 
# finaldata<- cbind(protein, finaldata)
# finaldata[1:5,1:5]
# 
# typeof(finaldata)
# row.names(finaldata)<- NULL
# 
# 
# finaldata[1:5,1:5]
# length(colnames(finaldata))
# #as.character(matrixdata)
# dim(finaldata)
# write.csv(finaldata, file = "proteinSum770_5fdr.csv", row.names = FALSE)
