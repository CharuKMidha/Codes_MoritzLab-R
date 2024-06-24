#Plate sepatation
rm(list = ls())

#TRANSPIOSE DATA take care of peptide or proteins

TRANSPOSEDATA<-function(data){
  colproteinsname<-data[,2] 
  rowsamples = names(data[,2:ncol(data)])
  data<- data[,2:ncol(data)]
  rownames(data)=c()
  colnames(data)=c()
  data
  tdata<- t(data)
  tdata <- data.frame(tdata)
  dim(tdata)
  tdata
  colnames(tdata)<-colproteinsname
  tdata
  tdata<- cbind(rowsamples, tdata)
  tdata
  return(tdata)
}

#protein#1, removed contams and Debrujin decoy

plate123<- read.csv("plate3.csv", header = TRUE, sep = ',')
dim(plate123)
plate123[1:5,]
##Bogdan data has these column names
# plate123$Abundance..126 <- NULL
# plate123$Modifications<- NULL
# plate123$File.ID<- NULL
# plate123$Abundance..127N <- NULL #empty channel
# plate123$Abundance..127C<- NULL #carrier channel 
##Bogdan data has these column names ENDS

#DS output file after libra
plate123$libra1 <- NULL #empty channel
plate123$libra2 <- NULL #carrier channel

plate123<- (plate123[plate123$isolation_interference < .7,]) #Isolation Interference less than 70%
plate123$isolation_interference <- NULL
#plate123<- head(plate123[plate123[,6:16] < 1])
dim(plate123)
plate123[1:5,]

plate123<- (plate123[plate123$num_prots == 1,]) #proteotypic proteins
plate123$num_prots <- NULL


#non-phosphopeptides

plate123 <- plate123[grep("-", plate123$localized_mods_STY.79.9663), ]
plate123$localized_mods_STY.79.9663  <- NULL

plate123<- plate123[- grep ("DEBRUIJN", plate123$protein), ]
plate123<- plate123[- grep ("CONTAM", plate123$protein), ]
dim(plate123)
plate123[1:5,]
library (stringr)
plate123$spectrum<- str_replace(plate123$spectrum,"190510X_JvsU_targeted15K_1x1_","")
plate123$spectrum<- str_replace(plate123$spectrum,"\\..*",".")
plate123[1:5,]
tail(plate123)
# For Bogdan's file 
#Added 1 to all
# min(plate123[,4:11], na.rm = TRUE ) 
# 
# allplatesplus1<- plate123[,4:11] +1 
# allplatesplus1[1:5,]
# #substitute NA to 1
# allplatesplus1[is.na(allplatesplus1)] <- 1
# 
# dim(allplatesplus1)
#  For Bogdan's file 


#allplatesplus1$rowsum <-rowSums(allplatesplus1[,1:10])
#allplatesplus1$carrierfilter <- allplatesplus1$Abundance..127C *100/ allplatesplus1$rowsum
#allplatesplus1[1:5,]
#allplatesplus1$rowsum <- NULL
#allplates<- cbind.data.frame(plate123[,1:3],allplatesplus1)
#allplates[1:5,]
#dim(allplates)


#1 Remove contams and decoys if you have any
#allplates<- allplates[- grep("CONTAM_sp", allplates$Protein.Accessions), ]
#dim(allplates)
#allplates<- allplates[- grep("DEBRUIJN", allplates$Protein.Accessions), ]
#dim(allplates)

#allplates[1:5,]

#allplatescarrierlessthan70<- (allplates[allplates$carrierfilter< 70,]) #less than 70%
#dim(allplatescarrierlessthan70)

#EXtract plates 
#plate1 <- allplatescarrierlessthan70[grep("1x1_Plate3", allplatescarrierlessthan70$Spectrum.File), ]
#head(plate1)
#plate2 <- allplatescarrierlessthan70[grep("2x2_Plate2", allplatescarrierlessthan70$Spectrum.File), ]
#plate3 <- allplatescarrierlessthan70[grep("4x4_Plate1", allplatescarrierlessthan70$Spectrum.File), ]





#1 Remove contams and decoys
#datanocontam<- plate1[- grep("CONTAM_sp", plate1$protein), ]
#datanocontamdecoy<- datanocontam[- grep("DEBRUIJN", datanocontam$protein), ]
#dim(datanocontam)
#dim(datanocontamdecoy)
#head(datanocontamdecoy)



#EXtract plates 
# plate3 <- allplates[grep("1x1_Plate3", allplates$Spectrum.File), ]
# head(plate3)
# plate2 <- allplates[grep("2x2_Plate2", allplates$Spectrum.File), ]
# plate1 <- allplates[grep("4x4_Plate1", allplates$Spectrum.File), ]
# 
# dim(plate1)
# dim(plate2)
# dim(plate3)
#2 Extract nonphospho peptides
#nonphosphopeptides <- datanocontamdecoy[grep("-", datanocontamdecoy$norm_info_gain_STY.79.9663), ]
#dim(nonphosphopeptides)
#phosphopeptides <- datanocontamdecoy[- grep("-", datanocontamdecoy$norm_info_gain_STY.79.9663), ]
#nonphosphopeptides$protein = NULL
#nonphosphopeptides$Empty127N = NULL
#nonphosphopeptides$Carrier127C = NULL
#nonphosphopeptides$norm_info_gain_STY.79.9663 = NULL
#dim(nonphosphopeptides)

plate1[1:5,]
plate2[1:5,]
plate3[1:5,]


#NON-PHOSPHO peptides
#nonphosphopep1 <- plate1[grep("-", plate1$norm_info_gain_STY.79.9663), ]
#dim(nonphosphopep1)

#####PLATE1 STARTS
npplate1<- cbind.data.frame(plate1$Spectrum.File,
                            plate1$Annotated.Sequence,
                            plate1$Abundance..128N,
                            plate1$Abundance..128C,#have to corrected
                            plate1$Abundance..129N,
                            plate1$Abundance..129C,
                            plate1$Abundance..130N,
                            plate1$Abundance..130C,
                            #plate1$Abundance..131N #dont include,
                            plate1$Abundance..131C,
                            plate1$Protein.Accessions
                
)

npplate1[1:5,1:10]

dim(npplate1)
tail(npplate1)

colnames(npplate1)<- c("spectrum","peptide",
                       "J128N","J128C",
                       "J129N","J129C",
                       "U130N","U130C",
                       "U131C", "protein")



#Extract the groups in each plate
npplate1[1:5,]

head(npplate1$spectrum)
library (stringr)
npplate1$spectrum<- str_replace(npplate1$spectrum,"190510X_JvsU_targeted15K_4x4_","")
npplate1$spectrum<- str_replace(npplate1$spectrum,"raw","")
npplate1[1:5,]

plate1_1 <- npplate1[grep("Plate1_1.", npplate1$spectrum), ]
plate1_1 <- plate1_1[- grep("Plate1_10.", plate1_1$spectrum), ]
plate1_1 <- plate1_1[- grep("Plate1_11.", plate1_1$spectrum), ]
plate1_1 <- plate1_1[- grep("Plate1_12.", plate1_1$spectrum), ]

plate1_2 <- npplate1[grep("Plate1_2.", npplate1$spectrum), ]
plate1_3 <- npplate1[grep("Plate1_3.", npplate1$spectrum), ]
plate1_4 <- npplate1[grep("Plate1_4.", npplate1$spectrum), ]
plate1_5 <- npplate1[grep("Plate1_5.", npplate1$spectrum), ]
plate1_6 <- npplate1[grep("Plate1_6.", npplate1$spectrum), ]
plate1_7 <- npplate1[grep("Plate1_7.", npplate1$spectrum), ]
plate1_8 <- npplate1[grep("Plate1_8.", npplate1$spectrum), ]
plate1_9 <- npplate1[grep("Plate1_9.", npplate1$spectrum), ]
plate1_10 <- npplate1[grep("Plate1_10.", npplate1$spectrum), ]
plate1_11 <- npplate1[grep("Plate1_11.", npplate1$spectrum), ]
plate1_12 <- npplate1[grep("Plate1_12.", npplate1$spectrum), ]

dim(plate1_1)
dim(plate1_2)
dim(plate1_3)
dim(plate1_4)
dim(plate1_5)
dim(plate1_6)
dim(plate1_7)
dim(plate1_8)
dim(plate1_8)
dim(plate1_9)
dim(plate1_10)
dim(plate1_11)
dim(plate1_12)

head(plate1_1)

plate1_1$spectrum1<- paste(plate1_1$spectrum,"J128N",sep = "" )
plate1_1$spectrum2<- paste(plate1_1$spectrum,"J128C",sep = "" )
plate1_1$spectrum3<- paste(plate1_1$spectrum,"J129N",sep = "" )
plate1_1$spectrum4<- paste(plate1_1$spectrum,"J129C",sep = "" )
plate1_1$spectrum5<- paste(plate1_1$spectrum,"U130N",sep = "" )
plate1_1$spectrum6<- paste(plate1_1$spectrum,"U130C",sep = "" )
plate1_1$spectrum7<- paste(plate1_1$spectrum,"U131C",sep = "" )

plate1_2$spectrum1<- paste(plate1_2$spectrum,"J128N",sep = "" )
plate1_2$spectrum2<- paste(plate1_2$spectrum,"J128C",sep = "" )
plate1_2$spectrum3<- paste(plate1_2$spectrum,"J129N",sep = "" )
plate1_2$spectrum4<- paste(plate1_2$spectrum,"J129C",sep = "" )
plate1_2$spectrum5<- paste(plate1_2$spectrum,"U130N",sep = "" )
plate1_2$spectrum6<- paste(plate1_2$spectrum,"U130C",sep = "" )
plate1_2$spectrum7<- paste(plate1_2$spectrum,"U131C",sep = "" )

plate1_3$spectrum1<- paste(plate1_3$spectrum,"J128N",sep = "" )
plate1_3$spectrum2<- paste(plate1_3$spectrum,"J128C",sep = "" )
plate1_3$spectrum3<- paste(plate1_3$spectrum,"J129N",sep = "" )
plate1_3$spectrum4<- paste(plate1_3$spectrum,"J129C",sep = "" )
plate1_3$spectrum5<- paste(plate1_3$spectrum,"U130N",sep = "" )
plate1_3$spectrum6<- paste(plate1_3$spectrum,"U130C",sep = "" )
plate1_3$spectrum7<- paste(plate1_3$spectrum,"U131C",sep = "" )


plate1_4$spectrum1<- paste(plate1_4$spectrum,"J128N",sep = "" )
plate1_4$spectrum2<- paste(plate1_4$spectrum,"J128C",sep = "" )
plate1_4$spectrum3<- paste(plate1_4$spectrum,"J129N",sep = "" )
plate1_4$spectrum4<- paste(plate1_4$spectrum,"J129C",sep = "" )
plate1_4$spectrum5<- paste(plate1_4$spectrum,"U130N",sep = "" )
plate1_4$spectrum6<- paste(plate1_4$spectrum,"U130C",sep = "" )
plate1_4$spectrum7<- paste(plate1_4$spectrum,"U131C",sep = "" )

plate1_5$spectrum1<- paste(plate1_5$spectrum,"J128N",sep = "" )
plate1_5$spectrum2<- paste(plate1_5$spectrum,"J128C",sep = "" )
plate1_5$spectrum3<- paste(plate1_5$spectrum,"J129N",sep = "" )
plate1_5$spectrum4<- paste(plate1_5$spectrum,"J129C",sep = "" )
plate1_5$spectrum5<- paste(plate1_5$spectrum,"U130N",sep = "" )
plate1_5$spectrum6<- paste(plate1_5$spectrum,"U130C",sep = "" )
plate1_5$spectrum7<- paste(plate1_5$spectrum,"U131C",sep = "" )

plate1_6$spectrum1<- paste(plate1_6$spectrum,"J128N",sep = "" )
plate1_6$spectrum2<- paste(plate1_6$spectrum,"J128C",sep = "" )
plate1_6$spectrum3<- paste(plate1_6$spectrum,"J129N",sep = "" )
plate1_6$spectrum4<- paste(plate1_6$spectrum,"J129C",sep = "" )
plate1_6$spectrum5<- paste(plate1_6$spectrum,"U130N",sep = "" )
plate1_6$spectrum6<- paste(plate1_6$spectrum,"U130C",sep = "" )
plate1_6$spectrum7<- paste(plate1_6$spectrum,"U131C",sep = "" )


plate1_7$spectrum1<- paste(plate1_7$spectrum,"J128N",sep = "" )
plate1_7$spectrum2<- paste(plate1_7$spectrum,"J128C",sep = "" )
plate1_7$spectrum3<- paste(plate1_7$spectrum,"J129N",sep = "" )
plate1_7$spectrum4<- paste(plate1_7$spectrum,"J129C",sep = "" )
plate1_7$spectrum5<- paste(plate1_7$spectrum,"U130N",sep = "" )
plate1_7$spectrum6<- paste(plate1_7$spectrum,"U130C",sep = "" )
plate1_7$spectrum7<- paste(plate1_7$spectrum,"U131C",sep = "" )

plate1_8$spectrum1<- paste(plate1_8$spectrum,"J128N",sep = "" )
plate1_8$spectrum2<- paste(plate1_8$spectrum,"J128C",sep = "" )
plate1_8$spectrum3<- paste(plate1_8$spectrum,"J129N",sep = "" )
plate1_8$spectrum4<- paste(plate1_8$spectrum,"J129C",sep = "" )
plate1_8$spectrum5<- paste(plate1_8$spectrum,"U130N",sep = "" )
plate1_8$spectrum6<- paste(plate1_8$spectrum,"U130C",sep = "" )
plate1_8$spectrum7<- paste(plate1_8$spectrum,"U131C",sep = "" )


plate1_9$spectrum1<- paste(plate1_9$spectrum,"J128N",sep = "" )
plate1_9$spectrum2<- paste(plate1_9$spectrum,"J128C",sep = "" )
plate1_9$spectrum3<- paste(plate1_9$spectrum,"J129N",sep = "" )
plate1_9$spectrum4<- paste(plate1_9$spectrum,"J129C",sep = "" )
plate1_9$spectrum5<- paste(plate1_9$spectrum,"U130N",sep = "" )
plate1_9$spectrum6<- paste(plate1_9$spectrum,"U130C",sep = "" )
plate1_9$spectrum7<- paste(plate1_9$spectrum,"U131C",sep = "" )

plate1_10$spectrum1<- paste(plate1_10$spectrum,"J128N",sep = "" )
plate1_10$spectrum2<- paste(plate1_10$spectrum,"J128C",sep = "" )
plate1_10$spectrum3<- paste(plate1_10$spectrum,"J129N",sep = "" )
plate1_10$spectrum4<- paste(plate1_10$spectrum,"J129C",sep = "" )
plate1_10$spectrum5<- paste(plate1_10$spectrum,"U130N",sep = "" )
plate1_10$spectrum6<- paste(plate1_10$spectrum,"U130C",sep = "" )
plate1_10$spectrum7<- paste(plate1_10$spectrum,"U131C",sep = "" )

plate1_11$spectrum1<- paste(plate1_11$spectrum,"J128N",sep = "" )
plate1_11$spectrum2<- paste(plate1_11$spectrum,"J128C",sep = "" )
plate1_11$spectrum3<- paste(plate1_11$spectrum,"J129N",sep = "" )
plate1_11$spectrum4<- paste(plate1_11$spectrum,"J129C",sep = "" )
plate1_11$spectrum5<- paste(plate1_11$spectrum,"U130N",sep = "" )
plate1_11$spectrum6<- paste(plate1_11$spectrum,"U130C",sep = "" )
plate1_11$spectrum7<- paste(plate1_11$spectrum,"U131C",sep = "" )

plate1_12$spectrum1<- paste(plate1_12$spectrum,"J128N",sep = "" )
plate1_12$spectrum2<- paste(plate1_12$spectrum,"J128C",sep = "" )
plate1_12$spectrum3<- paste(plate1_12$spectrum,"J129N",sep = "" )
plate1_12$spectrum4<- paste(plate1_12$spectrum,"J129C",sep = "" )
plate1_12$spectrum5<- paste(plate1_12$spectrum,"U130N",sep = "" )
plate1_12$spectrum6<- paste(plate1_12$spectrum,"U130C",sep = "" )
plate1_12$spectrum7<- paste(plate1_12$spectrum,"U131C",sep = "" )

allspectrum<- c(
  plate1_1$spectrum1,plate1_1$spectrum2,
  plate1_1$spectrum3,plate1_1$spectrum4,
  plate1_1$spectrum5,plate1_1$spectrum6,
  plate1_1$spectrum7,
  plate1_2$spectrum1,plate1_2$spectrum2,
  plate1_2$spectrum3,plate1_2$spectrum4,
  plate1_2$spectrum5,plate1_2$spectrum6,
  plate1_2$spectrum7,
  plate1_3$spectrum1,plate1_3$spectrum2,
  plate1_3$spectrum3,plate1_3$spectrum4,
  plate1_3$spectrum5,plate1_3$spectrum6,
  plate1_3$spectrum7,
  plate1_4$spectrum1,plate1_4$spectrum2,
  plate1_4$spectrum3,plate1_4$spectrum4,
  plate1_4$spectrum5,plate1_4$spectrum6,
  plate1_4$spectrum7,
  plate1_5$spectrum1,plate1_5$spectrum2,
  plate1_5$spectrum3,plate1_5$spectrum4,
  plate1_5$spectrum5,plate1_5$spectrum6,
  plate1_5$spectrum7,
  plate1_6$spectrum1,plate1_6$spectrum2,
  plate1_6$spectrum3,plate1_6$spectrum4,
  plate1_6$spectrum5,plate1_6$spectrum6,
  plate1_6$spectrum7,
  plate1_7$spectrum1,plate1_7$spectrum2,
  plate1_7$spectrum3,plate1_7$spectrum4,
  plate1_7$spectrum5,plate1_7$spectrum6,
  plate1_7$spectrum7,
  plate1_8$spectrum1,plate1_8$spectrum2,
  plate1_8$spectrum3,plate1_8$spectrum4,
  plate1_8$spectrum5,plate1_8$spectrum6,
  plate1_8$spectrum7,
  plate1_9$spectrum1,plate1_9$spectrum2,
  plate1_9$spectrum3,plate1_9$spectrum4,
  plate1_9$spectrum5,plate1_9$spectrum6,
  plate1_9$spectrum7,
  plate1_10$spectrum1,plate1_10$spectrum2,
  plate1_10$spectrum3,plate1_10$spectrum4,
  plate1_10$spectrum5,plate1_10$spectrum6,
  plate1_10$spectrum7,
  plate1_11$spectrum1,plate1_11$spectrum2,
  plate1_11$spectrum3,plate1_11$spectrum4,
  plate1_11$spectrum5,plate1_11$spectrum6,
  plate1_11$spectrum7,
  plate1_12$spectrum1,plate1_12$spectrum2,
  plate1_12$spectrum3,plate1_12$spectrum4,
  plate1_12$spectrum5,plate1_12$spectrum6,
  plate1_12$spectrum7
)

length(allspectrum)



peptide<- c(plate1_1[,2], plate1_1[,2],
            plate1_1[,2],plate1_1[,2],
            plate1_1[,2],plate1_1[,2],
            plate1_1[,2],
            plate1_2[,2], plate1_2[,2],
            plate1_2[,2],plate1_2[,2],
            plate1_2[,2],plate1_2[,2],
            plate1_2[,2],
            plate1_3[,2], plate1_3[,2],
            plate1_3[,2],plate1_3[,2],
            plate1_3[,2],plate1_3[,2],
            plate1_3[,2],
            plate1_4[,2], plate1_4[,2],
            plate1_4[,2],plate1_4[,2],
            plate1_4[,2],plate1_4[,2],
            plate1_4[,2],
            plate1_5[,2], plate1_5[,2],
            plate1_5[,2],plate1_5[,2],
            plate1_5[,2],plate1_5[,2],
            plate1_5[,2],
            plate1_6[,2], plate1_6[,2],
            plate1_6[,2],plate1_6[,2],
            plate1_6[,2],plate1_6[,2],
            plate1_6[,2],
            plate1_7[,2], plate1_7[,2],
            plate1_7[,2],plate1_7[,2],
            plate1_7[,2],plate1_7[,2],
            plate1_7[,2],
            plate1_8[,2], plate1_8[,2],
            plate1_8[,2],plate1_8[,2],
            plate1_8[,2],plate1_8[,2],
            plate1_8[,2],
            plate1_9[,2], plate1_9[,2],
            plate1_9[,2],plate1_9[,2],
            plate1_9[,2],plate1_9[,2],
            plate1_9[,2],
            plate1_10[,2], plate1_10[,2],
            plate1_10[,2],plate1_10[,2],
            plate1_10[,2],plate1_10[,2],
            plate1_10[,2],
            plate1_11[,2], plate1_11[,2],
            plate1_11[,2],plate1_11[,2],
            plate1_11[,2],plate1_11[,2],
            plate1_11[,2],
            plate1_12[,2], plate1_12[,2],
            plate1_12[,2],plate1_12[,2],
            plate1_12[,2],plate1_12[,2],
            plate1_12[,2]
)

length(peptide)

plate1values<- c(plate1_1$J128N, plate1_1$J128C,
                 plate1_1$J129N,plate1_1$J129C,
                 plate1_1$U130N, plate1_1$U130C,
                 plate1_1$U131C,
                 plate1_2$J128N, plate1_2$J128C,
                 plate1_2$J129N,plate1_2$J129C,
                 plate1_2$U130N, plate1_2$U130C,
                 plate1_2$U131C,
                 plate1_3$J128N, plate1_3$J128C,
                 plate1_3$J129N,plate1_3$J129C,
                 plate1_3$U130N, plate1_3$U130C,
                 plate1_3$U131C,
                 plate1_4$J128N, plate1_4$J128C,
                 plate1_4$J129N,plate1_4$J129C,
                 plate1_4$U130N, plate1_4$U130C,
                 plate1_4$U131C,
                 plate1_5$J128N, plate1_5$J128C,
                 plate1_5$J129N,plate1_5$J129C,
                 plate1_5$U130N, plate1_5$U130C,
                 plate1_5$U131C,
                 plate1_6$J128N, plate1_6$J128C,
                 plate1_6$J129N,plate1_6$J129C,
                 plate1_6$U130N, plate1_6$U130C,
                 plate1_6$U131C,
                 plate1_7$J128N, plate1_7$J128C,
                 plate1_7$J129N,plate1_7$J129C,
                 plate1_7$U130N, plate1_7$U130C,
                 plate1_7$U131C,
                 plate1_8$J128N, plate1_8$J128C,
                 plate1_8$J129N,plate1_8$J129C,
                 plate1_8$U130N, plate1_8$U130C,
                 plate1_8$U131C,
                 plate1_9$J128N, plate1_9$J128C,
                 plate1_9$J129N,plate1_9$J129C,
                 plate1_9$U130N, plate1_9$U130C,
                 plate1_9$U131C,
                 plate1_10$J128N, plate1_10$J128C,
                 plate1_10$J129N,plate1_10$J129C,
                 plate1_10$U130N, plate1_10$U130C,
                 plate1_10$U131C,
                 plate1_11$J128N, plate1_11$J128C,
                 plate1_11$J129N,plate1_11$J129C,
                 plate1_11$U130N, plate1_11$U130C,
                 plate1_11$U131C,
                 plate1_12$J128N, plate1_12$J128C,
                 plate1_12$J129N,plate1_12$J129C,
                 plate1_12$U130N, plate1_12$U130C,
                 plate1_12$U131C
                 
)

protein<- c(plate1_1[,10], plate1_1[,10],
            plate1_1[,10],plate1_1[,10],
            plate1_1[,10],plate1_1[,10],
            plate1_1[,10],
            plate1_2[,10], plate1_2[,10],
            plate1_2[,10],plate1_2[,10],
            plate1_2[,10],plate1_2[,10],
            plate1_2[,10],
            plate1_3[,10], plate1_3[,10],
            plate1_3[,10],plate1_3[,10],
            plate1_3[,10],plate1_3[,10],
            plate1_3[,10],
            plate1_4[,10], plate1_4[,10],
            plate1_4[,10],plate1_4[,10],
            plate1_4[,10],plate1_4[,10],
            plate1_4[,10],
            plate1_5[,10], plate1_5[,10],
            plate1_5[,10],plate1_5[,10],
            plate1_5[,10],plate1_5[,10],
            plate1_5[,10],
            plate1_6[,10], plate1_6[,10],
            plate1_6[,10],plate1_6[,10],
            plate1_6[,10],plate1_6[,10],
            plate1_6[,10],
            plate1_7[,10], plate1_7[,10],
            plate1_7[,10],plate1_7[,10],
            plate1_7[,10],plate1_7[,10],
            plate1_7[,10],
            plate1_8[,10], plate1_8[,10],
            plate1_8[,10],plate1_8[,10],
            plate1_8[,10],plate1_8[,10],
            plate1_8[,10],
            plate1_9[,10], plate1_9[,10],
            plate1_9[,10],plate1_9[,10],
            plate1_9[,10],plate1_9[,10],
            plate1_9[,10],
            plate1_10[,10], plate1_10[,10],
            plate1_10[,10],plate1_10[,10],
            plate1_10[,10],plate1_10[,10],
            plate1_10[,10],
            plate1_11[,10], plate1_11[,10],
            plate1_11[,10],plate1_11[,10],
            plate1_11[,10],plate1_11[,10],
            plate1_11[,10],
            plate1_12[,10], plate1_12[,10],
            plate1_12[,10],plate1_12[,10],
            plate1_12[,10],plate1_12[,10],
            plate1_12[,10]
)

pepprot<- str_c(peptide,'#', protein)

plate1all<- cbind.data.frame(allspectrum,pepprot,plate1values)

tail(plate1all)
head (plate1all)
dim(plate1all)
plate1all[50:60,]

library(dplyr)
#plate1all<- plate1all[order(plate1all$protein),]


dim(plate1all)
plate1allU<- unique(plate1all)
dim(plate1allU)

plate1allnona<- (plate1allU[plate1allU$plate1values >1 ,]) #no value with 1 which I added
plate1allnona[1:5,]
plate1allnona$plate1values <- plate1allnona$plate1values - 1 


write.csv(plate1allU,"plate1column.csv", row.names = FALSE)
write.csv(plate1allnona,"plate1columnnona.csv", row.names = FALSE)

####################PLATE 1 ENDS


#####PLATE2 STARTS
rm(list = ls())
npplate2<- cbind.data.frame(plate2$Spectrum.File,
                            plate2$Annotated.Sequence,
                            plate2$Abundance..128N,
                            plate2$Abundance..128C,#have to corrected
                            plate2$Abundance..129N,
                            plate2$Abundance..129C,
                            plate2$Abundance..130N,
                            plate2$Abundance..130C,
                            #plate2$Abundance..131N #dont include,
                            plate2$Abundance..131C,
                            plate2$Protein.Accessions
                            
)

npplate2[1:5,1:10]

dim(npplate2)
tail(npplate2)

colnames(npplate2)<- c("spectrum","peptide",
                       "J128N","J128C",
                       "U129N","U129C",
                       "J130N","J130C",
                       "U131C", "protein")



#Extract the groups in each plate
npplate2[1:5,]

head(npplate2$spectrum)
library (stringr)
npplate2$spectrum<- str_replace(npplate2$spectrum,"190513X_JvsU_targeted15K_2x2_","")
npplate2$spectrum<- str_replace(npplate2$spectrum,"raw","")
npplate2[1:5,]

plate2_1 <- npplate2[grep("Plate2_1.", npplate2$spectrum), ]
plate2_1 <- plate2_1[- grep("Plate2_10.", plate2_1$spectrum), ]
plate2_1 <- plate2_1[- grep("Plate2_11.", plate2_1$spectrum), ]
plate2_1 <- plate2_1[- grep("Plate2_12.", plate2_1$spectrum), ]

plate2_2 <- npplate2[grep("Plate2_2.", npplate2$spectrum), ]
plate2_3 <- npplate2[grep("Plate2_3.", npplate2$spectrum), ]
plate2_4 <- npplate2[grep("Plate2_4.", npplate2$spectrum), ]
plate2_5 <- npplate2[grep("Plate2_5.", npplate2$spectrum), ]
plate2_6 <- npplate2[grep("Plate2_6.", npplate2$spectrum), ]
plate2_7 <- npplate2[grep("Plate2_7.", npplate2$spectrum), ]
plate2_8 <- npplate2[grep("Plate2_8.", npplate2$spectrum), ]
plate2_9 <- npplate2[grep("Plate2_9.", npplate2$spectrum), ]
plate2_10 <- npplate2[grep("Plate2_10.", npplate2$spectrum), ]
plate2_11 <- npplate2[grep("Plate2_11.", npplate2$spectrum), ]
plate2_12 <- npplate2[grep("Plate2_12.", npplate2$spectrum), ]

dim(plate2_1)
dim(plate2_2)
dim(plate2_3)
dim(plate2_4)
dim(plate2_5)
dim(plate2_6)
dim(plate2_7)
dim(plate2_8)
dim(plate2_8)
dim(plate2_9)
dim(plate2_10)
dim(plate2_11)
dim(plate2_12)

head(plate2_1)

plate2_1$spectrum1<- paste(plate2_1$spectrum,"J128N",sep = "" )
plate2_1$spectrum2<- paste(plate2_1$spectrum,"J128C",sep = "" )
plate2_1$spectrum3<- paste(plate2_1$spectrum,"U129N",sep = "" )
plate2_1$spectrum4<- paste(plate2_1$spectrum,"U129C",sep = "" )
plate2_1$spectrum5<- paste(plate2_1$spectrum,"J130N",sep = "" )
plate2_1$spectrum6<- paste(plate2_1$spectrum,"J130C",sep = "" )
plate2_1$spectrum7<- paste(plate2_1$spectrum,"U131C",sep = "" )

plate2_2$spectrum1<- paste(plate2_2$spectrum,"J128N",sep = "" )
plate2_2$spectrum2<- paste(plate2_2$spectrum,"J128C",sep = "" )
plate2_2$spectrum3<- paste(plate2_2$spectrum,"U129N",sep = "" )
plate2_2$spectrum4<- paste(plate2_2$spectrum,"U129C",sep = "" )
plate2_2$spectrum5<- paste(plate2_2$spectrum,"J130N",sep = "" )
plate2_2$spectrum6<- paste(plate2_2$spectrum,"J30C",sep = "" )
plate2_2$spectrum7<- paste(plate2_2$spectrum,"U131C",sep = "" )

plate2_3$spectrum1<- paste(plate2_3$spectrum,"J128N",sep = "" )
plate2_3$spectrum2<- paste(plate2_3$spectrum,"J128C",sep = "" )
plate2_3$spectrum3<- paste(plate2_3$spectrum,"U129N",sep = "" )
plate2_3$spectrum4<- paste(plate2_3$spectrum,"U129C",sep = "" )
plate2_3$spectrum5<- paste(plate2_3$spectrum,"J130N",sep = "" )
plate2_3$spectrum6<- paste(plate2_3$spectrum,"J130C",sep = "" )
plate2_3$spectrum7<- paste(plate2_3$spectrum,"U131C",sep = "" )


plate2_4$spectrum1<- paste(plate2_4$spectrum,"J128N",sep = "" )
plate2_4$spectrum2<- paste(plate2_4$spectrum,"J128C",sep = "" )
plate2_4$spectrum3<- paste(plate2_4$spectrum,"U129N",sep = "" )
plate2_4$spectrum4<- paste(plate2_4$spectrum,"U129C",sep = "" )
plate2_4$spectrum5<- paste(plate2_4$spectrum,"J130N",sep = "" )
plate2_4$spectrum6<- paste(plate2_4$spectrum,"J130C",sep = "" )
plate2_4$spectrum7<- paste(plate2_4$spectrum,"U131C",sep = "" )

plate2_5$spectrum1<- paste(plate2_5$spectrum,"J128N",sep = "" )
plate2_5$spectrum2<- paste(plate2_5$spectrum,"J128C",sep = "" )
plate2_5$spectrum3<- paste(plate2_5$spectrum,"U129N",sep = "" )
plate2_5$spectrum4<- paste(plate2_5$spectrum,"U129C",sep = "" )
plate2_5$spectrum5<- paste(plate2_5$spectrum,"J130N",sep = "" )
plate2_5$spectrum6<- paste(plate2_5$spectrum,"J130C",sep = "" )
plate2_5$spectrum7<- paste(plate2_5$spectrum,"U131C",sep = "" )

plate2_6$spectrum1<- paste(plate2_6$spectrum,"J128N",sep = "" )
plate2_6$spectrum2<- paste(plate2_6$spectrum,"J128C",sep = "" )
plate2_6$spectrum3<- paste(plate2_6$spectrum,"U129N",sep = "" )
plate2_6$spectrum4<- paste(plate2_6$spectrum,"U129C",sep = "" )
plate2_6$spectrum5<- paste(plate2_6$spectrum,"J130N",sep = "" )
plate2_6$spectrum6<- paste(plate2_6$spectrum,"J130C",sep = "" )
plate2_6$spectrum7<- paste(plate2_6$spectrum,"U131C",sep = "" )


plate2_7$spectrum1<- paste(plate2_7$spectrum,"J128N",sep = "" )
plate2_7$spectrum2<- paste(plate2_7$spectrum,"J128C",sep = "" )
plate2_7$spectrum3<- paste(plate2_7$spectrum,"U129N",sep = "" )
plate2_7$spectrum4<- paste(plate2_7$spectrum,"U129C",sep = "" )
plate2_7$spectrum5<- paste(plate2_7$spectrum,"J130N",sep = "" )
plate2_7$spectrum6<- paste(plate2_7$spectrum,"J130C",sep = "" )
plate2_7$spectrum7<- paste(plate2_7$spectrum,"U131C",sep = "" )

plate2_8$spectrum1<- paste(plate2_8$spectrum,"J128N",sep = "" )
plate2_8$spectrum2<- paste(plate2_8$spectrum,"J128C",sep = "" )
plate2_8$spectrum3<- paste(plate2_8$spectrum,"U129N",sep = "" )
plate2_8$spectrum4<- paste(plate2_8$spectrum,"U129C",sep = "" )
plate2_8$spectrum5<- paste(plate2_8$spectrum,"J130N",sep = "" )
plate2_8$spectrum6<- paste(plate2_8$spectrum,"J130C",sep = "" )
plate2_8$spectrum7<- paste(plate2_8$spectrum,"U131C",sep = "" )


plate2_9$spectrum1<- paste(plate2_9$spectrum,"J128N",sep = "" )
plate2_9$spectrum2<- paste(plate2_9$spectrum,"J128C",sep = "" )
plate2_9$spectrum3<- paste(plate2_9$spectrum,"U129N",sep = "" )
plate2_9$spectrum4<- paste(plate2_9$spectrum,"U129C",sep = "" )
plate2_9$spectrum5<- paste(plate2_9$spectrum,"J130N",sep = "" )
plate2_9$spectrum6<- paste(plate2_9$spectrum,"J130C",sep = "" )
plate2_9$spectrum7<- paste(plate2_9$spectrum,"U131C",sep = "" )

plate2_10$spectrum1<- paste(plate2_10$spectrum,"J128N",sep = "" )
plate2_10$spectrum2<- paste(plate2_10$spectrum,"J128C",sep = "" )
plate2_10$spectrum3<- paste(plate2_10$spectrum,"U129N",sep = "" )
plate2_10$spectrum4<- paste(plate2_10$spectrum,"U129C",sep = "" )
plate2_10$spectrum5<- paste(plate2_10$spectrum,"J130N",sep = "" )
plate2_10$spectrum6<- paste(plate2_10$spectrum,"J130C",sep = "" )
plate2_10$spectrum7<- paste(plate2_10$spectrum,"U131C",sep = "" )

plate2_11$spectrum1<- paste(plate2_11$spectrum,"J128N",sep = "" )
plate2_11$spectrum2<- paste(plate2_11$spectrum,"J128C",sep = "" )
plate2_11$spectrum3<- paste(plate2_11$spectrum,"U129N",sep = "" )
plate2_11$spectrum4<- paste(plate2_11$spectrum,"U129C",sep = "" )
plate2_11$spectrum5<- paste(plate2_11$spectrum,"J130N",sep = "" )
plate2_11$spectrum6<- paste(plate2_11$spectrum,"J130C",sep = "" )
plate2_11$spectrum7<- paste(plate2_11$spectrum,"U131C",sep = "" )

plate2_12$spectrum1<- paste(plate2_12$spectrum,"J128N",sep = "" )
plate2_12$spectrum2<- paste(plate2_12$spectrum,"J128C",sep = "" )
plate2_12$spectrum3<- paste(plate2_12$spectrum,"U129N",sep = "" )
plate2_12$spectrum4<- paste(plate2_12$spectrum,"U129C",sep = "" )
plate2_12$spectrum5<- paste(plate2_12$spectrum,"J130N",sep = "" )
plate2_12$spectrum6<- paste(plate2_12$spectrum,"J130C",sep = "" )
plate2_12$spectrum7<- paste(plate2_12$spectrum,"U131C",sep = "" )

allspectrum<- c(
  plate2_1$spectrum1,plate2_1$spectrum2,
  plate2_1$spectrum3,plate2_1$spectrum4,
  plate2_1$spectrum5,plate2_1$spectrum6,
  plate2_1$spectrum7,
  plate2_2$spectrum1,plate2_2$spectrum2,
  plate2_2$spectrum3,plate2_2$spectrum4,
  plate2_2$spectrum5,plate2_2$spectrum6,
  plate2_2$spectrum7,
  plate2_3$spectrum1,plate2_3$spectrum2,
  plate2_3$spectrum3,plate2_3$spectrum4,
  plate2_3$spectrum5,plate2_3$spectrum6,
  plate2_3$spectrum7,
  plate2_4$spectrum1,plate2_4$spectrum2,
  plate2_4$spectrum3,plate2_4$spectrum4,
  plate2_4$spectrum5,plate2_4$spectrum6,
  plate2_4$spectrum7,
  plate2_5$spectrum1,plate2_5$spectrum2,
  plate2_5$spectrum3,plate2_5$spectrum4,
  plate2_5$spectrum5,plate2_5$spectrum6,
  plate2_5$spectrum7,
  plate2_6$spectrum1,plate2_6$spectrum2,
  plate2_6$spectrum3,plate2_6$spectrum4,
  plate2_6$spectrum5,plate2_6$spectrum6,
  plate2_6$spectrum7,
  plate2_7$spectrum1,plate2_7$spectrum2,
  plate2_7$spectrum3,plate2_7$spectrum4,
  plate2_7$spectrum5,plate2_7$spectrum6,
  plate2_7$spectrum7,
  plate2_8$spectrum1,plate2_8$spectrum2,
  plate2_8$spectrum3,plate2_8$spectrum4,
  plate2_8$spectrum5,plate2_8$spectrum6,
  plate2_8$spectrum7,
  plate2_9$spectrum1,plate2_9$spectrum2,
  plate2_9$spectrum3,plate2_9$spectrum4,
  plate2_9$spectrum5,plate2_9$spectrum6,
  plate2_9$spectrum7,
  plate2_10$spectrum1,plate2_10$spectrum2,
  plate2_10$spectrum3,plate2_10$spectrum4,
  plate2_10$spectrum5,plate2_10$spectrum6,
  plate2_10$spectrum7,
  plate2_11$spectrum1,plate2_11$spectrum2,
  plate2_11$spectrum3,plate2_11$spectrum4,
  plate2_11$spectrum5,plate2_11$spectrum6,
  plate2_11$spectrum7,
  plate2_12$spectrum1,plate2_12$spectrum2,
  plate2_12$spectrum3,plate2_12$spectrum4,
  plate2_12$spectrum5,plate2_12$spectrum6,
  plate2_12$spectrum7
)

length(allspectrum)



peptide<- c(plate2_1[,2], plate2_1[,2],
            plate2_1[,2],plate2_1[,2],
            plate2_1[,2],plate2_1[,2],
            plate2_1[,2],
            plate2_2[,2], plate2_2[,2],
            plate2_2[,2],plate2_2[,2],
            plate2_2[,2],plate2_2[,2],
            plate2_2[,2],
            plate2_3[,2], plate2_3[,2],
            plate2_3[,2],plate2_3[,2],
            plate2_3[,2],plate2_3[,2],
            plate2_3[,2],
            plate2_4[,2], plate2_4[,2],
            plate2_4[,2],plate2_4[,2],
            plate2_4[,2],plate2_4[,2],
            plate2_4[,2],
            plate2_5[,2], plate2_5[,2],
            plate2_5[,2],plate2_5[,2],
            plate2_5[,2],plate2_5[,2],
            plate2_5[,2],
            plate2_6[,2], plate2_6[,2],
            plate2_6[,2],plate2_6[,2],
            plate2_6[,2],plate2_6[,2],
            plate2_6[,2],
            plate2_7[,2], plate2_7[,2],
            plate2_7[,2],plate2_7[,2],
            plate2_7[,2],plate2_7[,2],
            plate2_7[,2],
            plate2_8[,2], plate2_8[,2],
            plate2_8[,2],plate2_8[,2],
            plate2_8[,2],plate2_8[,2],
            plate2_8[,2],
            plate2_9[,2], plate2_9[,2],
            plate2_9[,2],plate2_9[,2],
            plate2_9[,2],plate2_9[,2],
            plate2_9[,2],
            plate2_10[,2], plate2_10[,2],
            plate2_10[,2],plate2_10[,2],
            plate2_10[,2],plate2_10[,2],
            plate2_10[,2],
            plate2_11[,2], plate2_11[,2],
            plate2_11[,2],plate2_11[,2],
            plate2_11[,2],plate2_11[,2],
            plate2_11[,2],
            plate2_12[,2], plate2_12[,2],
            plate2_12[,2],plate2_12[,2],
            plate2_12[,2],plate2_12[,2],
            plate2_12[,2]
)

length(peptide)

plate2values<- c(plate2_1$J128N, plate2_1$J128C,
                 plate2_1$U129N,plate2_1$U129C,
                 plate2_1$J130N, plate2_1$J130C,
                 plate2_1$U131C,
                 plate2_2$J128N, plate2_2$J128C,
                 plate2_2$U129N,plate2_2$U129C,
                 plate2_2$J130N, plate2_2$J130C,
                 plate2_2$U131C,
                 plate2_3$J128N, plate2_3$J128C,
                 plate2_3$U129N,plate2_3$U129C,
                 plate2_3$J130N, plate2_3$J130C,
                 plate2_3$U131C,
                 plate2_4$J128N, plate2_4$J128C,
                 plate2_4$U129N,plate2_4$U129C,
                 plate2_4$J130N, plate2_4$J130C,
                 plate2_4$U131C,
                 plate2_5$J128N, plate2_5$J128C,
                 plate2_5$U129N,plate2_5$U129C,
                 plate2_5$J130N, plate2_5$J130C,
                 plate2_5$U131C,
                 plate2_6$J128N, plate2_6$J128C,
                 plate2_6$U129N,plate2_6$U129C,
                 plate2_6$J130N, plate2_6$J130C,
                 plate2_6$U131C,
                 plate2_7$J128N, plate2_7$J128C,
                 plate2_7$U129N,plate2_7$U129C,
                 plate2_7$J130N, plate2_7$J130C,
                 plate2_7$U131C,
                 plate2_8$J128N, plate2_8$J128C,
                 plate2_8$U129N,plate2_8$U129C,
                 plate2_8$J130N, plate2_8$J130C,
                 plate2_8$U131C,
                 plate2_9$J128N, plate2_9$J128C,
                 plate2_9$U129N,plate2_9$U129C,
                 plate2_9$J130N, plate2_9$J130C,
                 plate2_9$U131C,
                 plate2_10$J128N, plate2_10$J128C,
                 plate2_10$U129N,plate2_10$U129C,
                 plate2_10$J130N, plate2_10$J130C,
                 plate2_10$U131C,
                 plate2_11$J128N, plate2_11$J128C,
                 plate2_11$U129N,plate2_11$U129C,
                 plate2_11$J130N, plate2_11$J130C,
                 plate2_11$U131C,
                 plate2_12$J128N, plate2_12$J128C,
                 plate2_12$U129N,plate2_12$U129C,
                 plate2_12$J130N, plate2_12$J130C,
                 plate2_12$U131C
                 
)

protein<- c(plate2_1[,10], plate2_1[,10],
            plate2_1[,10],plate2_1[,10],
            plate2_1[,10],plate2_1[,10],
            plate2_1[,10],
            plate2_2[,10], plate2_2[,10],
            plate2_2[,10],plate2_2[,10],
            plate2_2[,10],plate2_2[,10],
            plate2_2[,10],
            plate2_3[,10], plate2_3[,10],
            plate2_3[,10],plate2_3[,10],
            plate2_3[,10],plate2_3[,10],
            plate2_3[,10],
            plate2_4[,10], plate2_4[,10],
            plate2_4[,10],plate2_4[,10],
            plate2_4[,10],plate2_4[,10],
            plate2_4[,10],
            plate2_5[,10], plate2_5[,10],
            plate2_5[,10],plate2_5[,10],
            plate2_5[,10],plate2_5[,10],
            plate2_5[,10],
            plate2_6[,10], plate2_6[,10],
            plate2_6[,10],plate2_6[,10],
            plate2_6[,10],plate2_6[,10],
            plate2_6[,10],
            plate2_7[,10], plate2_7[,10],
            plate2_7[,10],plate2_7[,10],
            plate2_7[,10],plate2_7[,10],
            plate2_7[,10],
            plate2_8[,10], plate2_8[,10],
            plate2_8[,10],plate2_8[,10],
            plate2_8[,10],plate2_8[,10],
            plate2_8[,10],
            plate2_9[,10], plate2_9[,10],
            plate2_9[,10],plate2_9[,10],
            plate2_9[,10],plate2_9[,10],
            plate2_9[,10],
            plate2_10[,10], plate2_10[,10],
            plate2_10[,10],plate2_10[,10],
            plate2_10[,10],plate2_10[,10],
            plate2_10[,10],
            plate2_11[,10], plate2_11[,10],
            plate2_11[,10],plate2_11[,10],
            plate2_11[,10],plate2_11[,10],
            plate2_11[,10],
            plate2_12[,10], plate2_12[,10],
            plate2_12[,10],plate2_12[,10],
            plate2_12[,10],plate2_12[,10],
            plate2_12[,10]
)

pepprot<- str_c(peptide,'#', protein)

plate2all<- cbind.data.frame(allspectrum,pepprot,plate2values)

tail(plate2all)
head (plate2all)
dim(plate2all)
plate2all[50:60,]

library(dplyr)
#plate2all<- plate2all[order(plate2all$protein),]


dim(plate2all)
plate2allU<- unique(plate2all)
dim(plate2allU)

plate2allnona<- (plate2allU[plate2allU$plate2values >1 ,]) #no value with 1 which I added
plate2allnona[1:5,]
plate2allnona$plate2values <- plate2allnona$plate2values - 1 


write.csv(plate2allU,"plate2column.csv", row.names = FALSE)
write.csv(plate2allnona,"plate2columnnona.csv", row.names = FALSE)

################PLATE 2 ENDS HERE

#####PLATE 3 STARTS
#rm(list = ls())
#Removes row with <1 value in DS libra output
plate123<- filter(plate123, plate123$libra3 >=1 ) #value more than 1 
plate123<- filter(plate123, plate123$libra4 >=1 ) #value more than 1 
plate123<- filter(plate123, plate123$libra5 >=1 ) #value more than 1 
plate123<- filter(plate123, plate123$libra6 >=1 ) #value more than 1 
plate123<- filter(plate123, plate123$libra7 >=1 ) #value more than 1 
plate123<- filter(plate123, plate123$libra8 >=1 ) #value more than 1 
plate123<- filter(plate123, plate123$libra9 >=1 ) #value more than 1 
plate123<- filter(plate123, plate123$libra10 >=1 ) #value more than 1 

dim(plate123)
npplate3<- cbind.data.frame(plate123$spectrum,
                            plate123$peptide,
                            plate123$libra3,
                            plate123$libra4,#have to corrected
                            plate123$libra5,
                            plate123$libra6,
                            plate123$libra7,
                            plate123$libra8,
                            #plate123$libra9 #dont include,
                            plate123$libra10,
                            plate123$protein
                            
)

npplate3[1:5,1:10]

dim(npplate3)
tail(npplate3)

colnames(npplate3)<- c("spectrum","peptide",
                       "J128N","U128C",
                       "J129N","U129C",
                       "J130N","U130C",
                       "U131C", "protein")



#Extract the groups in each plate
npplate3[1:5,]

head(npplate3$spectrum)

npplate3[1:5,]

plate3_1 <- npplate3[grep("Plate3_1.", npplate3$spectrum), ]
plate3_1 <- plate3_1[- grep("Plate3_10.", plate3_1$spectrum), ]
plate3_1 <- plate3_1[- grep("Plate3_11.", plate3_1$spectrum), ]
plate3_1 <- plate3_1[- grep("Plate3_12.", plate3_1$spectrum), ]

plate3_2 <- npplate3[grep("Plate3_2.", npplate3$spectrum), ]
plate3_3 <- npplate3[grep("Plate3_3.", npplate3$spectrum), ]
plate3_4 <- npplate3[grep("Plate3_4.", npplate3$spectrum), ]
plate3_5 <- npplate3[grep("Plate3_5.", npplate3$spectrum), ]
plate3_6 <- npplate3[grep("Plate3_6.", npplate3$spectrum), ]
plate3_7 <- npplate3[grep("Plate3_7.", npplate3$spectrum), ]
plate3_8 <- npplate3[grep("Plate3_8.", npplate3$spectrum), ]
plate3_9 <- npplate3[grep("Plate3_9.", npplate3$spectrum), ]
plate3_10 <- npplate3[grep("Plate3_10.", npplate3$spectrum), ]
plate3_11 <- npplate3[grep("Plate3_11.", npplate3$spectrum), ]
plate3_12 <- npplate3[grep("Plate3_12.", npplate3$spectrum), ]

dim(plate3_1)
dim(plate3_2)
dim(plate3_3)
dim(plate3_4)
dim(plate3_5)
dim(plate3_6)
dim(plate3_7)
dim(plate3_8)
dim(plate3_8)
dim(plate3_9)
dim(plate3_10)
dim(plate3_11)
dim(plate3_12)

head(plate3_1)

plate3_1$spectrum1<- paste(plate3_1$spectrum,"J128N",sep = "" )
plate3_1$spectrum2<- paste(plate3_1$spectrum,"U128C",sep = "" )
plate3_1$spectrum3<- paste(plate3_1$spectrum,"J129N",sep = "" )
plate3_1$spectrum4<- paste(plate3_1$spectrum,"U129C",sep = "" )
plate3_1$spectrum5<- paste(plate3_1$spectrum,"J130N",sep = "" )
plate3_1$spectrum6<- paste(plate3_1$spectrum,"U130C",sep = "" )
plate3_1$spectrum7<- paste(plate3_1$spectrum,"U131C",sep = "" )

plate3_2$spectrum1<- paste(plate3_2$spectrum,"J128N",sep = "" )
plate3_2$spectrum2<- paste(plate3_2$spectrum,"U128C",sep = "" )
plate3_2$spectrum3<- paste(plate3_2$spectrum,"J129N",sep = "" )
plate3_2$spectrum4<- paste(plate3_2$spectrum,"U129C",sep = "" )
plate3_2$spectrum5<- paste(plate3_2$spectrum,"J130N",sep = "" )
plate3_2$spectrum6<- paste(plate3_2$spectrum,"U130C",sep = "" )
plate3_2$spectrum7<- paste(plate3_2$spectrum,"U131C",sep = "" )

plate3_3$spectrum1<- paste(plate3_3$spectrum,"J128N",sep = "" )
plate3_3$spectrum2<- paste(plate3_3$spectrum,"U128C",sep = "" )
plate3_3$spectrum3<- paste(plate3_3$spectrum,"J129N",sep = "" )
plate3_3$spectrum4<- paste(plate3_3$spectrum,"U129C",sep = "" )
plate3_3$spectrum5<- paste(plate3_3$spectrum,"J130N",sep = "" )
plate3_3$spectrum6<- paste(plate3_3$spectrum,"U130C",sep = "" )
plate3_3$spectrum7<- paste(plate3_3$spectrum,"U131C",sep = "" )


plate3_4$spectrum1<- paste(plate3_4$spectrum,"J128N",sep = "" )
plate3_4$spectrum2<- paste(plate3_4$spectrum,"U128C",sep = "" )
plate3_4$spectrum3<- paste(plate3_4$spectrum,"J129N",sep = "" )
plate3_4$spectrum4<- paste(plate3_4$spectrum,"U129C",sep = "" )
plate3_4$spectrum5<- paste(plate3_4$spectrum,"J130N",sep = "" )
plate3_4$spectrum6<- paste(plate3_4$spectrum,"U130C",sep = "" )
plate3_4$spectrum7<- paste(plate3_4$spectrum,"U131C",sep = "" )

plate3_5$spectrum1<- paste(plate3_5$spectrum,"J128N",sep = "" )
plate3_5$spectrum2<- paste(plate3_5$spectrum,"U128C",sep = "" )
plate3_5$spectrum3<- paste(plate3_5$spectrum,"J129N",sep = "" )
plate3_5$spectrum4<- paste(plate3_5$spectrum,"U129C",sep = "" )
plate3_5$spectrum5<- paste(plate3_5$spectrum,"J130N",sep = "" )
plate3_5$spectrum6<- paste(plate3_5$spectrum,"U130C",sep = "" )
plate3_5$spectrum7<- paste(plate3_5$spectrum,"U131C",sep = "" )

plate3_6$spectrum1<- paste(plate3_6$spectrum,"J128N",sep = "" )
plate3_6$spectrum2<- paste(plate3_6$spectrum,"U128C",sep = "" )
plate3_6$spectrum3<- paste(plate3_6$spectrum,"J129N",sep = "" )
plate3_6$spectrum4<- paste(plate3_6$spectrum,"U129C",sep = "" )
plate3_6$spectrum5<- paste(plate3_6$spectrum,"J130N",sep = "" )
plate3_6$spectrum6<- paste(plate3_6$spectrum,"U130C",sep = "" )
plate3_6$spectrum7<- paste(plate3_6$spectrum,"U131C",sep = "" )


plate3_7$spectrum1<- paste(plate3_7$spectrum,"J128N",sep = "" )
plate3_7$spectrum2<- paste(plate3_7$spectrum,"U128C",sep = "" )
plate3_7$spectrum3<- paste(plate3_7$spectrum,"J129N",sep = "" )
plate3_7$spectrum4<- paste(plate3_7$spectrum,"U129C",sep = "" )
plate3_7$spectrum5<- paste(plate3_7$spectrum,"J130N",sep = "" )
plate3_7$spectrum6<- paste(plate3_7$spectrum,"U130C",sep = "" )
plate3_7$spectrum7<- paste(plate3_7$spectrum,"U131C",sep = "" )

plate3_8$spectrum1<- paste(plate3_8$spectrum,"J128N",sep = "" )
plate3_8$spectrum2<- paste(plate3_8$spectrum,"U128C",sep = "" )
plate3_8$spectrum3<- paste(plate3_8$spectrum,"J129N",sep = "" )
plate3_8$spectrum4<- paste(plate3_8$spectrum,"U129C",sep = "" )
plate3_8$spectrum5<- paste(plate3_8$spectrum,"J130N",sep = "" )
plate3_8$spectrum6<- paste(plate3_8$spectrum,"U130C",sep = "" )
plate3_8$spectrum7<- paste(plate3_8$spectrum,"U131C",sep = "" )


plate3_9$spectrum1<- paste(plate3_9$spectrum,"J128N",sep = "" )
plate3_9$spectrum2<- paste(plate3_9$spectrum,"U128C",sep = "" )
plate3_9$spectrum3<- paste(plate3_9$spectrum,"J129N",sep = "" )
plate3_9$spectrum4<- paste(plate3_9$spectrum,"U129C",sep = "" )
plate3_9$spectrum5<- paste(plate3_9$spectrum,"J130N",sep = "" )
plate3_9$spectrum6<- paste(plate3_9$spectrum,"U130C",sep = "" )
plate3_9$spectrum7<- paste(plate3_9$spectrum,"U131C",sep = "" )

plate3_10$spectrum1<- paste(plate3_10$spectrum,"J128N",sep = "" )
plate3_10$spectrum2<- paste(plate3_10$spectrum,"U128C",sep = "" )
plate3_10$spectrum3<- paste(plate3_10$spectrum,"J129N",sep = "" )
plate3_10$spectrum4<- paste(plate3_10$spectrum,"U129C",sep = "" )
plate3_10$spectrum5<- paste(plate3_10$spectrum,"J130N",sep = "" )
plate3_10$spectrum6<- paste(plate3_10$spectrum,"U130C",sep = "" )
plate3_10$spectrum7<- paste(plate3_10$spectrum,"U131C",sep = "" )

plate3_11$spectrum1<- paste(plate3_11$spectrum,"J128N",sep = "" )
plate3_11$spectrum2<- paste(plate3_11$spectrum,"U128C",sep = "" )
plate3_11$spectrum3<- paste(plate3_11$spectrum,"J129N",sep = "" )
plate3_11$spectrum4<- paste(plate3_11$spectrum,"U129C",sep = "" )
plate3_11$spectrum5<- paste(plate3_11$spectrum,"J130N",sep = "" )
plate3_11$spectrum6<- paste(plate3_11$spectrum,"U130C",sep = "" )
plate3_11$spectrum7<- paste(plate3_11$spectrum,"U131C",sep = "" )

plate3_12$spectrum1<- paste(plate3_12$spectrum,"J128N",sep = "" )
plate3_12$spectrum2<- paste(plate3_12$spectrum,"U128C",sep = "" )
plate3_12$spectrum3<- paste(plate3_12$spectrum,"J129N",sep = "" )
plate3_12$spectrum4<- paste(plate3_12$spectrum,"U129C",sep = "" )
plate3_12$spectrum5<- paste(plate3_12$spectrum,"J130N",sep = "" )
plate3_12$spectrum6<- paste(plate3_12$spectrum,"U130C",sep = "" )
plate3_12$spectrum7<- paste(plate3_12$spectrum,"U131C",sep = "" )

allspectrum<- c(
  plate3_1$spectrum1,plate3_1$spectrum2,
  plate3_1$spectrum3,plate3_1$spectrum4,
  plate3_1$spectrum5,plate3_1$spectrum6,
  plate3_1$spectrum7,
  plate3_2$spectrum1,plate3_2$spectrum2,
  plate3_2$spectrum3,plate3_2$spectrum4,
  plate3_2$spectrum5,plate3_2$spectrum6,
  plate3_2$spectrum7,
  plate3_3$spectrum1,plate3_3$spectrum2,
  plate3_3$spectrum3,plate3_3$spectrum4,
  plate3_3$spectrum5,plate3_3$spectrum6,
  plate3_3$spectrum7,
  plate3_4$spectrum1,plate3_4$spectrum2,
  plate3_4$spectrum3,plate3_4$spectrum4,
  plate3_4$spectrum5,plate3_4$spectrum6,
  plate3_4$spectrum7,
  plate3_5$spectrum1,plate3_5$spectrum2,
  plate3_5$spectrum3,plate3_5$spectrum4,
  plate3_5$spectrum5,plate3_5$spectrum6,
  plate3_5$spectrum7,
  plate3_6$spectrum1,plate3_6$spectrum2,
  plate3_6$spectrum3,plate3_6$spectrum4,
  plate3_6$spectrum5,plate3_6$spectrum6,
  plate3_6$spectrum7,
  plate3_7$spectrum1,plate3_7$spectrum2,
  plate3_7$spectrum3,plate3_7$spectrum4,
  plate3_7$spectrum5,plate3_7$spectrum6,
  plate3_7$spectrum7,
  plate3_8$spectrum1,plate3_8$spectrum2,
  plate3_8$spectrum3,plate3_8$spectrum4,
  plate3_8$spectrum5,plate3_8$spectrum6,
  plate3_8$spectrum7,
  plate3_9$spectrum1,plate3_9$spectrum2,
  plate3_9$spectrum3,plate3_9$spectrum4,
  plate3_9$spectrum5,plate3_9$spectrum6,
  plate3_9$spectrum7,
  plate3_10$spectrum1,plate3_10$spectrum2,
  plate3_10$spectrum3,plate3_10$spectrum4,
  plate3_10$spectrum5,plate3_10$spectrum6,
  plate3_10$spectrum7,
  plate3_11$spectrum1,plate3_11$spectrum2,
  plate3_11$spectrum3,plate3_11$spectrum4,
  plate3_11$spectrum5,plate3_11$spectrum6,
  plate3_11$spectrum7,
  plate3_12$spectrum1,plate3_12$spectrum2,
  plate3_12$spectrum3,plate3_12$spectrum4,
  plate3_12$spectrum5,plate3_12$spectrum6,
  plate3_12$spectrum7
)

length(allspectrum)



peptide<- c(plate3_1[,2], plate3_1[,2],
            plate3_1[,2],plate3_1[,2],
            plate3_1[,2],plate3_1[,2],
            plate3_1[,2],
            plate3_2[,2], plate3_2[,2],
            plate3_2[,2],plate3_2[,2],
            plate3_2[,2],plate3_2[,2],
            plate3_2[,2],
            plate3_3[,2], plate3_3[,2],
            plate3_3[,2],plate3_3[,2],
            plate3_3[,2],plate3_3[,2],
            plate3_3[,2],
            plate3_4[,2], plate3_4[,2],
            plate3_4[,2],plate3_4[,2],
            plate3_4[,2],plate3_4[,2],
            plate3_4[,2],
            plate3_5[,2], plate3_5[,2],
            plate3_5[,2],plate3_5[,2],
            plate3_5[,2],plate3_5[,2],
            plate3_5[,2],
            plate3_6[,2], plate3_6[,2],
            plate3_6[,2],plate3_6[,2],
            plate3_6[,2],plate3_6[,2],
            plate3_6[,2],
            plate3_7[,2], plate3_7[,2],
            plate3_7[,2],plate3_7[,2],
            plate3_7[,2],plate3_7[,2],
            plate3_7[,2],
            plate3_8[,2], plate3_8[,2],
            plate3_8[,2],plate3_8[,2],
            plate3_8[,2],plate3_8[,2],
            plate3_8[,2],
            plate3_9[,2], plate3_9[,2],
            plate3_9[,2],plate3_9[,2],
            plate3_9[,2],plate3_9[,2],
            plate3_9[,2],
            plate3_10[,2], plate3_10[,2],
            plate3_10[,2],plate3_10[,2],
            plate3_10[,2],plate3_10[,2],
            plate3_10[,2],
            plate3_11[,2], plate3_11[,2],
            plate3_11[,2],plate3_11[,2],
            plate3_11[,2],plate3_11[,2],
            plate3_11[,2],
            plate3_12[,2], plate3_12[,2],
            plate3_12[,2],plate3_12[,2],
            plate3_12[,2],plate3_12[,2],
            plate3_12[,2]
)

length(peptide)

plate3values<- c(plate3_1$J128N, plate3_1$U128C,
                 plate3_1$J129N,plate3_1$U129C,
                 plate3_1$J130N, plate3_1$U130C,
                 plate3_1$U131C,
                 plate3_2$J128N, plate3_2$U128C,
                 plate3_2$J129N,plate3_2$U129C,
                 plate3_2$J130N, plate3_2$U130C,
                 plate3_2$U131C,
                 plate3_3$J128N, plate3_3$U128C,
                 plate3_3$J129N,plate3_3$U129C,
                 plate3_3$J130N, plate3_3$U130C,
                 plate3_3$U131C,
                 plate3_4$J128N, plate3_4$U128C,
                 plate3_4$J129N,plate3_4$U129C,
                 plate3_4$J130N, plate3_4$U130C,
                 plate3_4$U131C,
                 plate3_5$J128N, plate3_5$U128C,
                 plate3_5$J129N,plate3_5$U129C,
                 plate3_5$J130N, plate3_5$U130C,
                 plate3_5$U131C,
                 plate3_6$J128N, plate3_6$U128C,
                 plate3_6$J129N,plate3_6$U129C,
                 plate3_6$J130N, plate3_6$U130C,
                 plate3_6$U131C,
                 plate3_7$J128N, plate3_7$U128C,
                 plate3_7$J129N,plate3_7$U129C,
                 plate3_7$J130N, plate3_7$U130C,
                 plate3_7$U131C,
                 plate3_8$J128N, plate3_8$U128C,
                 plate3_8$J129N,plate3_8$U129C,
                 plate3_8$J130N, plate3_8$U130C,
                 plate3_8$U131C,
                 plate3_9$J128N, plate3_9$U128C,
                 plate3_9$J129N,plate3_9$U129C,
                 plate3_9$J130N, plate3_9$U130C,
                 plate3_9$U131C,
                 plate3_10$J128N, plate3_10$U128C,
                 plate3_10$J129N,plate3_10$U129C,
                 plate3_10$J130N, plate3_10$U130C,
                 plate3_10$U131C,
                 plate3_11$J128N, plate3_11$U128C,
                 plate3_11$J129N,plate3_11$U129C,
                 plate3_11$J130N, plate3_11$U130C,
                 plate3_11$U131C,
                 plate3_12$J128N, plate3_12$U128C,
                 plate3_12$J129N,plate3_12$U129C,
                 plate3_12$J130N, plate3_12$U130C,
                 plate3_12$U131C
                 
)

protein<- c(plate3_1[,10], plate3_1[,10],
            plate3_1[,10],plate3_1[,10],
            plate3_1[,10],plate3_1[,10],
            plate3_1[,10],
            plate3_2[,10], plate3_2[,10],
            plate3_2[,10],plate3_2[,10],
            plate3_2[,10],plate3_2[,10],
            plate3_2[,10],
            plate3_3[,10], plate3_3[,10],
            plate3_3[,10],plate3_3[,10],
            plate3_3[,10],plate3_3[,10],
            plate3_3[,10],
            plate3_4[,10], plate3_4[,10],
            plate3_4[,10],plate3_4[,10],
            plate3_4[,10],plate3_4[,10],
            plate3_4[,10],
            plate3_5[,10], plate3_5[,10],
            plate3_5[,10],plate3_5[,10],
            plate3_5[,10],plate3_5[,10],
            plate3_5[,10],
            plate3_6[,10], plate3_6[,10],
            plate3_6[,10],plate3_6[,10],
            plate3_6[,10],plate3_6[,10],
            plate3_6[,10],
            plate3_7[,10], plate3_7[,10],
            plate3_7[,10],plate3_7[,10],
            plate3_7[,10],plate3_7[,10],
            plate3_7[,10],
            plate3_8[,10], plate3_8[,10],
            plate3_8[,10],plate3_8[,10],
            plate3_8[,10],plate3_8[,10],
            plate3_8[,10],
            plate3_9[,10], plate3_9[,10],
            plate3_9[,10],plate3_9[,10],
            plate3_9[,10],plate3_9[,10],
            plate3_9[,10],
            plate3_10[,10], plate3_10[,10],
            plate3_10[,10],plate3_10[,10],
            plate3_10[,10],plate3_10[,10],
            plate3_10[,10],
            plate3_11[,10], plate3_11[,10],
            plate3_11[,10],plate3_11[,10],
            plate3_11[,10],plate3_11[,10],
            plate3_11[,10],
            plate3_12[,10], plate3_12[,10],
            plate3_12[,10],plate3_12[,10],
            plate3_12[,10],plate3_12[,10],
            plate3_12[,10]
)

pepprot<- str_c(peptide,'#', protein)

plate3all<- cbind.data.frame(allspectrum,pepprot,plate3values)

tail(plate3all)
head (plate3all)
dim(plate3all)
plate3all[50:60,]

library(dplyr)
#plate3all<- plate3all[order(plate3all$protein),]


dim(plate3all)
plate3allU<- unique(plate3all)
dim(plate3allU)

# plate3allnona<- (plate3allU[plate3allU$plate3values >1 ,]) #no value with 1 which I added
# plate3allnona[1:5,]
# plate3allnona$plate3values <- plate3allnona$plate3values - 1 


write.csv(plate3allU,"plate3column.csv", row.names = FALSE)
#write.csv(plate3allnona,"plate3columnnona.csv", row.names = FALSE)

################PLATE 3 ENDS HERE

#1)Open the above file and  remove the first V1... row
#2) run for other two plates combine them then 
#3) run perl script from David to remove the duplicates from above file


####################COMBINE ALL 123 plates ###############################

#-------------------------------------------------------------
#preprocessing matrix PLATE1

rm(list = ls())
####plate1columnU.csv.avgU.sum  plate2columnU.csv.avgU.sum plate3column.csv.avgU.sum
pretable<- read.csv("plate1columnU.csv.avgU.sum", header = TRUE, sep = "\t")
#pretable$pepprot <- str_replace(pretable$pepprot,"^.*#","")
#length(unique(pretable$pepprot))

#write.csv(pretable,"plate1columnprot.csv", row.names = FALSE)
install.packages("scales") 
library(tidyr)
head(pretable)
datauniq <- unique(pretable)
dim(pretable)
dim(datauniq)
#colnames(pretable) <- c("allspectrum", "pepprot", "values")
#head(pretable)
finaldata = datauniq %>%
  spread(key = pepprot,
         value = plate1values)

dim(finaldata)
head(t(finaldata))
write.csv(t(finaldata), file = "plate1pepprotreducematrix.csv", col.names = FALSE)
#Open the above file to remove V1 row and remove the . to _

#-------------------------------------------------------------
#preprocessing matrix PLATE2

rm(list = ls())
####plate1columnU.csv.avgU.sum  plate2columnU.csv.avgU.sum plate3column.csv.avgU.sum
pretable<- read.csv("plate2columnU.csv.avgU.sum", header = TRUE, sep = "\t")
#pretable$pepprot <- str_replace(pretable$pepprot,"^.*#","")
#length(unique(pretable$pepprot))

#write.csv(pretable,"plate1columnprot.csv", row.names = FALSE)
#install.packages("scales") 
library(tidyr)
head(pretable)
datauniq <- unique(pretable)
dim(pretable)
dim(datauniq)
#colnames(pretable) <- c("allspectrum", "pepprot", "values")
#head(pretable)
finaldata = datauniq %>%
  spread(key = pepprot,
         value = plate2values)

dim(finaldata)
head(t(finaldata))
write.csv(t(finaldata), file = "plate2pepprotreducematrix.csv", col.names = FALSE)
#Open the above file to remove V1 row and remove the . to _

#-------------------------------------------------------------
#preprocessing matrix PLATE3

rm(list = ls())
####plate1columnU.csv.avgU.sum  plate2columnU.csv.avgU.sum plate3column.csv.avgU.sum
pretable<- read.csv("plate3column.csv.gavg.sum", header = TRUE, sep = "\t")
#pretable$pepprot <- str_replace(pretable$pepprot,"^.*#","")
#length(unique(pretable$pepprot))

#write.csv(pretable,"plate1columnprot.csv", row.names = FALSE)
#install.packages("scales") 
library(tidyr)
head(pretable)
datauniq <- unique(pretable)
dim(pretable)
dim(datauniq)
#colnames(pretable) <- c("allspectrum", "pepprot", "values")
#head(pretable)
finaldata = datauniq %>%
  spread(key = pepprot,
         value = plate3values)

dim(finaldata)
head(t(finaldata))
write.csv(t(finaldata), file = "plate3pepprotreducematrix.csv", col.names = FALSE)
#Open the above file to remove V1 row and remove the . to _

#####################################     PCA ########################
rm(list = ls())
#TRANSPIOSE DATA

TRANSPOSEDATA<-function(data){
  colproteinsname<-data[,1]
  rowsamples = names(data[,2:ncol(data)])
  data<- data[,2:ncol(data)]
  rownames(data)=c()
  colnames(data)=c()
  data
  tdata<- t(data)
  tdata <- data.frame(tdata)
  dim(tdata)
  tdata
  colnames(tdata)<-colproteinsname
  tdata
  tdata<- cbind(rowsamples, tdata)
  tdata
  return(tdata)
}

library(ggplot2)
library(grid)
library(gridExtra)
library(tidyverse)
theme_set(theme_bw(16))
#dim(nonphosphopeptides)
#initialdata<- read.csv("allplatematrix.csv", header = TRUE, sep = ',')
initialdata<- read.csv("plate3pepprotreducematrix.csv", header = TRUE, sep = ',')

initialdata[1:5,]
dim(initialdata)

min(initialdata[,2:85], na.rm = TRUE)
#the matrix is like proteins as rows and conditions in the columns
#REMOVE 12 columns
# initialdata$Plate1_1_U131C= NULL
# initialdata$Plate1_12_J128C= NULL
# initialdata$Plate1_12_U130N= NULL
# initialdata$Plate1_5_J128C= NULL
# initialdata$Plate1_5_J129C= NULL
# initialdata$Plate2_1_U129C= NULL
# initialdata$Plate2_11_J128N= NULL
# initialdata$Plate2_12_J128C= NULL
# initialdata$Plate3_1_U128C= NULL
# initialdata$Plate3_1_U129C= NULL
# initialdata$Plate3_11_U129C= NULL
# initialdata$Plate3_4_J129N= NULL

# 
# 
# #Plate1 extract
# Plate1matrix <- cbind.data.frame(allspectrum= initialdata$allspectrum, initialdata[, grep(pattern="^Plate1_", colnames(initialdata))])
# dim(Plate1matrix)
# Plate1matrix<- na.omit(Plate1matrix)
# dim(Plate1matrix)
# 
# head(Plate1matrix)
# 
# #Plate2 extract
# Plate2matrix <- cbind.data.frame(allspectrum= initialdata$allspectrum, initialdata[, grep(pattern="^Plate2_", colnames(initialdata))])
# dim(Plate2matrix)
# Plate2matrix<- na.omit(Plate2matrix)
# dim(Plate2matrix)
# head(Plate2matrix)
# 
# #Plate3 extract
# Plate3matrix <- cbind.data.frame(allspectrum= initialdata$allspectrum, initialdata[, grep(pattern="^Plate3_", colnames(initialdata))])
# Plate3matrix<- na.omit(Plate3matrix)
# dim(Plate3matrix)
# head(Plate3matrix)

#allplates
head (initialdata)
dim(initialdata)
initialdata <- na.omit(initialdata)
dim (initialdata)
head(initialdata)

#Transpose the the matrix to make proteins as columns and conditions as Rows
#Tinitialdata<- TRANSPOSEDATA(initialdata2)


Tinitialdata<- TRANSPOSEDATA(initialdata)
#Tinitialdata<- TRANSPOSEDATA(Plate1matrix)
#Tinitialdata<- TRANSPOSEDATA(Plate2matrix)
#Tinitialdata<- TRANSPOSEDATA(Plate3matrix)



Tinitialdata[1:5,1:5]
dim(Tinitialdata)


#get the name of the conditions as row name

tags<- substr(sapply( strsplit(as.character(Tinitialdata$rowsamples), "_"), "[[", 3 ),1,1)
length(tags)
#row.names(Tinitialdata) <- paste(Tinitialdata$rowsamples, row.names(Tinitialdata), sep="_") 
row.names(Tinitialdata) <- paste(Tinitialdata$rowsamples, c(tags), sep="_") 

row.names(Tinitialdata)
#Remove the rowsamples column as we save this information in the rownames
Tinitialdata$rowsamples <- NULL
#class(Tinitialdata)
#class(Tinitialdata$`1`)

#write.csv(Tinitialdata, "Tdataout.csv")
sum(!is.finite(scale(Tinitialdata))) 
which(apply(Tinitialdata, 2, var)==0)

Tdatazerovar<- Tinitialdata[ , which(apply(Tinitialdata, 2, var) != 0)]
dim(Tdatazerovar)
Tdatazerovar[1:5,]


#we compute the PCA
Tinitialdata_pca<- prcomp(Tdatazerovar, center = TRUE, scale. = TRUE)
#summary(Tinitialdata_pca)
dim(Tinitialdata_pca)
Tinitialdata_pca$sdev
#Compute Variance Explained by Each PC
var_explained_df <- data.frame(PC= paste0("PC",1:length(tags)),
                               var_explained=round(((Tinitialdata_pca$sdev)^2/sum((Tinitialdata_pca$sdev)^2))*100, 1))

head(  var_explained_df)
#write.csv(var_explained_df, "pcascreeinfoPlate3.csv")
#Scree plot with line plot in R
var_explained_df %>%
  ggplot(aes(x=PC,y=var_explained, group=1))+
  geom_point(size=4)+
  geom_line()+
  labs(title="Scree plot: PCA on scaled data")

#png("screeplotnonphospho2x2.png")
#Screeplot with bar plot in R
var_explained_df %>%
  ggplot(aes(x=PC,y=var_explained))+
  geom_col()+
  labs(title="Scree plot: PCA on scaled data phospho")
#dev.off()


#Ploting PCA1 and PCA2

#plot(Tinitialdata_pca$x[,1], Tinitialdata_pca$x[,2])

##extract x from pca of data frame
Tinitialdata_out <- as.data.frame(Tinitialdata_pca$x)

#Add a new column after removing the replicate information
#Tinitialdata_out$Stages <- sapply( strsplit(as.character(row.names(Tinitialdata_order)), "_"), "[[", 1 )
Tinitialdata_out$plate<- sapply( strsplit(as.character(row.names(Tinitialdata)), "_"), "[[", 1 )
Tinitialdata_out$group<- sapply( strsplit(as.character(row.names(Tinitialdata)), "_"), "[[", 2 )
Tinitialdata_out$names<- sapply( strsplit(as.character(row.names(Tinitialdata)), "_"), "[[", 3 )
Tinitialdata_out$ju<- sapply( strsplit(as.character(row.names(Tinitialdata)), "_"), "[[", 4 )


Tinitialdata_out$plate
Tinitialdata_out$group
Tinitialdata_out$names
Tinitialdata_out$ju
#Plotting the pca based on color of stages or conditions
#png("nonphospho2x2.png", width = 1090, height = 675)
#install.packages("ggforce")

# ggplot(Tinitialdata_out,aes(x=PC1,y=PC2,color=ju) )+
#   geom_point( size=5,shape = 21, stroke = 2.5)+
#   geom_label( label = Tinitialdata_out$names,
#               nudge_x = 5, nudge_y = 6, 
#               check_overlap = T) +
#   geom_mark_ellipse(aes(color = plate,
#                         label=ju),
#                     expand = unit(0.5,"mm"),
#                     label.buffer = unit(-5, 'mm'))+
#   xlab(paste("PC1 ",round(var_explained_df$var_explained[1],0),"%", sep=""))+
#   ylab(paste("PC2 ", round(var_explained_df$var_explained[2],0),"%", sep=""))+
#   theme_bw() +
#   theme(text=element_text(size=30))+
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         panel.background = element_blank(), axis.line = element_line(colour = "black", size=2))
# 
# #dev.off()
#install.packages('ggforce')
library(ggforce)
'/
Tinitialdata_out %>%
  ggplot(aes(x = PC1,
             y = PC2))+
  geom_mark_ellipse(aes(fill = plate,
                        label = group),
                    expand = unit(0.5,"mm"),
                    label.buffer = unit(-5, 'mm'))+
  geom_point(aes(color = ju),labels = names )+
  theme(legend.position = "topright")


Tinitialdata_out %>%
  ggplot(aes(x = PC1,
             y = PC2))+
  geom_polygon(aes(fill = plate,
                   label = group),
               expand = unit(0.5,"mm"),
               label.buffer = unit(-5, 'mm'))+
  geom_point(aes(color = ju))+
  theme(legend.position = "topright")

install.packages("deldir")
Tinitialdata_out %>%
  ggplot(aes(x = PC1,
             y = PC2))+
  geom_delaunay_tile(alpha = 0.3) + 
  geom_delaunay_segment2(aes(colour = group, group = -1), size = 2,
                         lineend = 'round')
/'
#TRY THIS ONE

mytheme = list(
  theme_classic()+
    theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
          legend.title = element_blank(),legend.position="bottom", strip.text = element_text(face="bold", size=15),
          axis.text=element_text(face="bold", size=15),axis.title = element_text(face="bold", size=15),plot.title = element_text(face = "bold", hjust = 1,size=13))
)

Tinitialdata_out %>%
  ggplot(aes(x = PC1,
             y = PC2))+
  
  #geom_voronoi_tile(aes(fill = group, group = -1L), 
   #                 max.radius = 0.2)+ #colour = 'black')
  #geom_mark_ellipse(aes(fill = plate,
   #                     label = plate),
    #                expand = unit(0.5,"mm"),
     #               label.buffer = unit(-5, 'mm'))+
#  geom_point(aes(shape = factor(ju)),size = 8)+
  #geom_point(aes(shape = factor(plate)),size = 8)+
  geom_point(aes(colour = factor(ju)), size = 4)+
 # geom_point(aes(colour = factor(group)), size = 4)+
  ggtitle("Plate3(1x1)")+
  mytheme+
  geom_label( label = Tinitialdata_out$names,
              nudge_x = .2, nudge_y = 2, 
              check_overlap = T) +
  xlab(paste("PC1 ",round(var_explained_df$var_explained[1],0),"%", sep=""))+
  ylab(paste("PC2 ", round(var_explained_df$var_explained[2],0),"%", sep=""))

)
#ggsave("plate1.png")

#PLATE1 PCA


