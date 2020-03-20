#FileName	Peptide	Quantity
#K562inhouse_1_Good_SN_peptide	GTFIIDPGGVIR.2	188416.1406
#K562inhouse_1_Good_SN_peptide	GTFIIDPAAVIR.2	164120.8594
#K562inhouse_1_Good_SN_peptide	GTFIIDPAAVIR.3	3264.374512
#K562inhouse_1_Good_SN_peptide	YILAGVENSK.2	126287.7969
#K562inhouse_1_Good_SN_peptide	TPVITGAPYEYR.2	67901.125
#K562inhouse_1_Good_SN_peptide	TPVITGAPYEYR.3	2730.161133
#K562inhouse_1_Good_SN_peptide	DGLDAASYYAPVR.2	68274.65625
#K562inhouse_1_Good_SN_peptide	ADVTPADFSEWSK.2	55172.82031
#K562inhouse_1_Good_SN_peptide	LGGNEQVTR.2	43690.59375


data = read.csv("Allgradientprotein.txt", header = TRUE, sep = '\t')    



library(tidyr)
finaldata = data %>%
  spread(key = Protein,
         value = Quantity)


write.csv(t(finaldata), file = "Allgradientprotein.csv")


