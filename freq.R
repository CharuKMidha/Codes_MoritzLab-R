library(data.table)
#data_table <- read.table('proteotypicProtein_peptide.txt', sep = '\t', header = TRUE, fill = TRUE)
data_table <- read.table('fragment_precursordistributionQ1Q3.txt', sep = '\t', header = TRUE, fill = TRUE)

head (data_table)

#data_summary <- summary(data_table)
#tail(data_summary)

freq_list <- table(data_table$Q1)

freq_list

write.csv(freq_list, file="fragment_precursordistributionQ1Q3_freq.txt", sep = '\t')
