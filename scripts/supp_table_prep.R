## Env. Bioinformatics, Fall 2021
## Author: JBW
## Date created: 11/18/2021
## Date updated: 11/19/2021

setwd("Desktop/Macbook/WHOI/Classes/FALL_21/final-project")

######
# SUBSET OUT THE 47 SEQS WE ANALYZED

Supp_tab_1 <- read.csv("supp_table_1.csv", na.strings = c("n.a.",""))
  # saved relevant spreadsheet tab as supp_table_1.csv

Supp_tab_1_single <- subset.data.frame(Supp_tab_1, Supp_tab_1$Colonization_Status == "single")
  # remove any with colonization status of "multiple" or "incomplete"

Supp_tab_1_47 <- subset.data.frame(Supp_tab_1_single, Supp_tab_1_single$total_reads < 80000000)
  # remove the 2 samples with 10x the total reads

colnames(Supp_tab_1_47)[1] = "Sample_ID"
  # fix column name bc it was annoying me

#write.csv(Supp_tab_1_47, file = "supp_table_subset.csv",row.names = F)
  # write new csv file -- only need to do once


######
# RE-ORDER TO MATCH VCF FILE

poptab <- read.csv("supp_table_subset.csv")
  # just to be safe, load in the new (subsetted) csv file

seq_order = read.csv("vcf_seq_order.txt",sep='\n',header=F)
seq_order$index = seq(1:47)
  # load in the seq list from the vcf file header & create numbered index column

poptab$index = NA   #create index column for sample metadata table
for (i in 1:47){
  for (k in 1:47){
    if(poptab$NCBI_SRA_Accesion[i]==seq_order$V1[k]){ #if the SRA seq ID's match...
      poptab$index[i] = seq_order$index[k] #...then assign the vcf index number
    }
  }
}
poptab_ordered = poptab[order(poptab$index),] #re-order by the vcf index number

#write.csv(poptab_ordered, file = "supp_table_ordered.csv",row.names = F)
  # write the new csv -- only need to do once


