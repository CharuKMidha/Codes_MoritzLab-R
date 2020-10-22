#R.FileName	PG.ProteinGroups	PG.Quantity
#UI_6HRS_SWATH1	[ RT-Cal protein ]	29985.99609
#UI_6HRS_SWATH1	A0AVT1	320.7033691
#UI_18HRS_SWATH1	A0FGR8	354.5825806
#UI_18HRS_SWATH1	A0MZ66	101.2458038
#UI_30HRS_SWATH1	A2RUS2	59.19276428
#UI_30HRS_SWATH1	A5A3E0	12765.27637
#UI_42HRS_SWATH1	A5YKK6	547.6705322
#UI_42HRS_SWATH1	A5YKK6	547.6705322

data = read.csv("/K565/K562_Good_SN_RT2-3_Unique.txt", header = FALSE, sep ="\t")    
data

library(tidyr)
finaldata = data %>%
  spread(key = V2,
         value = V3)


write.csv(finaldata, file = "PHL_Good_SN_RT2-3_Unique_analysis_matrix.csv")


