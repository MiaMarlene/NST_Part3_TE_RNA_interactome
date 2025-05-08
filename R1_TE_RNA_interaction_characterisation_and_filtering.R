#Attempting to normalise the chimera reads
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("GenomicRanges")
#BiocManager::install("R4RNA")

library(rnaCrosslinkOO)
library(pheatmap)
library(ggplot2)
library(dplyr)
library(biomaRt)
library("circlize")
library(ggrepel)


setwd("C:/Users/miabe/OneDrive/Desktop/Part III Systems/Project")

#reading files in



blood_int_df = read.csv("Blood_int_df.csv", row.names = 1)
beagle_int_df = read.csv("beagle_int_df.csv", row.names = 1)
tS18_int_df = read.csv("tS18_int_df.csv", row.names =1)
total_interaction_df =read.csv("total_interaction_df.csv")





############################
##### Importing the data####
############################

#Chimeric reads input
sampleTabler1 = c("000Scripts/fullData/T1_S5.assembled.tstk_comp_Dmel_hyb_hybrids.hyb", "s", "1", "s1")
sampleTabler2 = c("000Scripts/fullData/T1C_S6.assembled.tstk_comp_Dmel_hyb_hybrids.hyb", "c", "1", "c1")
sampleTabler3 = c("000Scripts/fullData/T2_S7.assembled.tstk_comp_Dmel_hyb_hybrids.hyb", "c", "2", "c2")
sampleTabler4 = c("000Scripts/fullData/T2_S7.assembled.tstk_comp_Dmel_hyb_hybrids.hyb", "s", "2", "s2")
sampleTabler5 = c("000Scripts/fullData/Tr3_S5.assembled.tstk_comp_Dmel_hyb_hybrids.hyb", "s", "3", "s3")
sampleTabler6 = c("000Scripts/fullData/Tr3-C_S6.assembled.tstk_comp_Dmel_hyb_hybrids.hyb", "c", "3", "c3")

# Make the sample table
sampleTable2 = rbind.data.frame(sampleTabler1, sampleTabler2,
                                sampleTabler3, sampleTabler4,
                                sampleTabler5,sampleTabler6)
#sampleTable2 = rbind.data.frame(sampleTabler1, sampleTabler4,
                                sampleTabler5,sampleTabler6)
# # Add the column names
colnames(sampleTable2) = c("file", "group", "sample", "sampleName")
#
RNA_types = c("rRNA", "miRNA","protein_coding","ncRNA","pseudogene","snoRNA","snRNA","tRNA")


TE_names = c("AY180916_AY180916_AY180916_AY180916","FBgn0001207_FBgn0001207_TEHMSbeagle_HMSbeagle","FBgn0005384_FBgn0005384_3S18_3S18")
#Make cds objects and feature counts
cdsBlood = rnaCrosslinkDataSet(rnas = "AY180916_AY180916_AY180916_AY180916",
                               rnaSize = 0,
                               sampleTable =     sampleTable2)

cdsHMSbeagle = rnaCrosslinkDataSet(rnas = "FBgn0001207_FBgn0001207_TEHMSbeagle_HMSbeagle",
                                   rnaSize = 0,
                                   sampleTable =     sampleTable2)

cds3S18 = rnaCrosslinkDataSet(rnas = "FBgn0005384_FBgn0005384_3S18_3S18",
                              rnaSize = 0,
                              sampleTable =     sampleTable2)



###############################################
#### Creating combined interaction dataframes ####
###############################################

# blood_int_lst = c(cdsBlood@InputFiles$AY180916_AY180916_AY180916_AY180916$host$s1$V10,
#               cdsBlood@InputFiles$AY180916_AY180916_AY180916_AY180916$host$s2$V10,
#               cdsBlood@InputFiles$AY180916_AY180916_AY180916_AY180916$host$s3$V10)
# blood_int_df = as.data.frame(table(blood_int_lst))
# colnames(blood_int_df) = c("RowNames", "blood_int")
# 
# beagle_int_lst = c(cdsHMSbeagle@InputFiles$FBgn0001207_FBgn0001207_TEHMSbeagle_HMSbeagle$host$s1$V10,
#                   cdsHMSbeagle@InputFiles$FBgn0001207_FBgn0001207_TEHMSbeagle_HMSbeagle$host$s2$V10,
#                   cdsHMSbeagle@InputFiles$FBgn0001207_FBgn0001207_TEHMSbeagle_HMSbeagle$host$s3$V10)
# beagle_int_df = as.data.frame(table(beagle_int_lst))
# colnames(beagle_int_df) = c("RowNames", "beagle_int")

# tS18_in_lst = c(cds3S18@InputFiles$FBgn0005384_FBgn0005384_3S18_3S18$host$s1$V10,
#                 cds3S18@InputFiles$FBgn0005384_FBgn0005384_3S18_3S18$host$s2$V10,
#                 cds3S18@InputFiles$FBgn0005384_FBgn0005384_3S18_3S18$host$s3$V10)
# tS18_int_df = as.data.frame(table(tS18_in_lst))
# colnames(tS18_int_df) = c("RowNames", "tS18_int")
# 


# 
# write.csv(blood_int_df, "Blood_int_df.csv")
write.csv(beagle_int_df, "beagle_int_df.csv")
# write.csv(tS18_int_df, "tS18_int_df.csv")


#Interaction data frames
blood_int_df = read.csv("Blood_int_df.csv", row.names = 1)
beagle_int_df = read.csv("beagle_int_df.csv", row.names = 1)
tS18_int_df = read.csv("tS18_int_df.csv", row.names =1)

# 
# 
#Total interaction data frames combined

total_ints1 = cdsBlood@InputFiles$all$all$s1
total_ints2 = cdsBlood@InputFiles$all$all$s2
total_ints3 = cdsBlood@InputFiles$all$all$s3

#Total inter-RNA interactions
inter_ints1_indices = which(total_ints1$V4 != total_ints1$V10)
inter_ints2_indices = which(total_ints2$V4 != total_ints2$V10)
inter_ints3_indices = which(total_ints3$V4 != total_ints3$V10)

total_inter_full_df = rbind(total_ints1[inter_ints1_indices,c(4,10)],
                            total_ints2[inter_ints2_indices,c(4,10)],
                            total_ints3[inter_ints3_indices,c(4,10)])

total_inter_table = as.data.frame(table(c(total_inter_full_df$V4,total_inter_full_df$V10)))
colnames(total_inter_table) = c("RowNames", "total_inter_int")

#inter-RNA interactions without TEs
total_inter_full_df


cat = total_inter_full_df$V10

total_inter_full_df = total_inter_full_df %>%
  mutate(
    category2 = case_when(
      grepl("miRNA", cat) ~ "miRNA",
      grepl("protein_coding", cat) ~ "protein_coding",
      grepl("3S18", cat)~ "3S18",
      grepl("AY180916", cat)~ "Blood",
      grepl("ncRNA", cat)~ "ncRNA",
      grepl("pseudogene", cat)~ "pseudogene",
      grepl("rRNA", cat)~ "rRNA",
      grepl("snoRNA", cat)~ "snoRNA",
      grepl("snRNA", cat)~ "snRNA",
      grepl("tRNA", cat)~ "tRNA",
      grepl("HMSbeagle", cat)~ "HMSbeagle",
      TRUE ~ "Other"
    ))


total_inter_categories = data.frame(interactor =c(total_inter_full_df$category1, total_inter_full_df$category2), interacting_with = c(total_inter_full_df$category2, total_inter_full_df$category1))

total_inter_categories_rRNA = total_inter_categories[total_inter_categories$interactor == "rRNA",]
total_inter_categories_snRNA = total_inter_categories[total_inter_categories$interactor == "snRNA",]
total_inter_categories_miRNA = total_inter_categories[total_inter_categories$interactor == "miRNA",]
total_inter_categories_protein_coding = total_inter_categories[total_inter_categories$interactor == "protein_coding",]
total_inter_categories_pseudogene = total_inter_categories[total_inter_categories$interactor == "pseudogene",]
total_inter_categories_tRNA = total_inter_categories[total_inter_categories$interactor == "tRNA",]
total_inter_categories_Blood = total_inter_categories[total_inter_categories$interactor == "Blood",]
total_inter_categories_HMSbeagle = total_inter_categories[total_inter_categories$interactor == "HMSbeagle",]
total_inter_categories_tS18 = total_inter_categories[total_inter_categories$interactor == "3S18",]
total_inter_categories_ncRNA = total_inter_categories[total_inter_categories$interactor == "ncRNA",]
total_inter_categories_snoRNA = total_inter_categories[total_inter_categories$interactor == "snoRNA",]


total_int_df_rRNA = data.frame(table(total_inter_categories_rRNA))
total_int_df_snRNA = data.frame(table(total_inter_categories_snRNA))
total_int_df_miRNA = data.frame(table(total_inter_categories_miRNA))
total_int_df_protein_coding = data.frame(table(total_inter_categories_protein_coding))
total_int_df_pseudogene = data.frame(table(total_inter_categories_pseudogene))
total_int_df_tRNA = data.frame(table(total_inter_categories_tRNA))
total_int_df_Blood = data.frame(table(total_inter_categories_Blood))
total_int_df_HMSbeagle = data.frame(table(total_inter_categories_HMSbeagle))
total_int_df_tS18 = data.frame(table(total_inter_categories_tS18))
total_int_df_ncRNA = data.frame(table(total_inter_categories_ncRNA))
total_int_df_snoRNA = data.frame(table(total_inter_categories_snoRNA))


total_int_cat_table = left_join(total_int_df_rRNA[2:3], total_int_df_snRNA[,2:3], by="interacting_with")
total_int_cat_table = left_join(total_int_cat_table, total_int_df_miRNA[,2:3], by="interacting_with")
total_int_cat_table = left_join(total_int_cat_table, total_int_df_protein_coding[,2:3], by="interacting_with")
total_int_cat_table = left_join(total_int_cat_table, total_int_df_pseudogene[,2:3], by="interacting_with")
total_int_cat_table = left_join(total_int_cat_table, total_int_df_tRNA[,2:3], by="interacting_with")
total_int_cat_table = left_join(total_int_cat_table, total_int_df_Blood[,2:3], by="interacting_with")
total_int_cat_table = left_join(total_int_cat_table, total_int_df_HMSbeagle[,2:3], by="interacting_with")
total_int_cat_table = left_join(total_int_cat_table, total_int_df_tS18[,2:3], by="interacting_with")
total_int_cat_table = left_join(total_int_cat_table, total_int_df_ncRNA[,2:3], by="interacting_with")
total_int_cat_table = left_join(total_int_cat_table, total_int_df_snoRNA[,2:3], by="interacting_with")
rownames(total_int_cat_table) = total_int_cat_table$interacting_with
total_int_cat_table = total_int_cat_table[2:12]
colnames(total_int_cat_table) = c("rRNA", "snRNA", "miRNA", "protein_coding", "pseudogene", "tRNA", "Blood", "HMSbeagle", "tS18", "ncRNA", "snoRNA")
total_int_cat_table[is.na(total_int_cat_table)] = 0

unique(total_inter_categories$interactor)

write.table(total_int_cat_table, "total_inter_categories.txt", row.names = FALSE, col.names = F)


excl_TE_df = total_inter_full_df[total_inter_full_df$V4 != "AY180916_AY180916_AY180916_AY180916" &
                                   total_inter_full_df$V10 != "AY180916_AY180916_AY180916_AY180916" &
                                   total_inter_full_df$V4 != "FBgn0001207_FBgn0001207_TEHMSbeagle_HMSbeagle" &
                                   total_inter_full_df$V10 != "FBgn0001207_FBgn0001207_TEHMSbeagle_HMSbeagle" &
                                   total_inter_full_df$V4 != "FBgn0005384_FBgn0005384_3S18_3S18" &
                                   total_inter_full_df$V10 != "FBgn0005384_FBgn0005384_3S18_3S18",]
head(excl_TE_df)
incl_TE_df = total_inter_full_df[total_inter_full_df$V4 == "AY180916_AY180916_AY180916_AY180916" |
                                   total_inter_full_df$V10 == "AY180916_AY180916_AY180916_AY180916" |
                                   total_inter_full_df$V4 == "FBgn0001207_FBgn0001207_TEHMSbeagle_HMSbeagle" |
                                   total_inter_full_df$V10 == "FBgn0001207_FBgn0001207_TEHMSbeagle_HMSbeagle" |
                                   total_inter_full_df$V4 == "FBgn0005384_FBgn0005384_3S18_3S18" |
                                   total_inter_full_df$V10 == "FBgn0005384_FBgn0005384_3S18_3S18",]
interTE_df_indices = which(incl_TE_df$V4 =="AY180916_AY180916_AY180916_AY180916" & incl_TE_df$V10 %in% c("FBgn0001207_FBgn0001207_TEHMSbeagle_HMSbeagle","FBgn0005384_FBgn0005384_3S18_3S18") |
                             incl_TE_df$V4 =="FBgn0001207_FBgn0001207_TEHMSbeagle_HMSbeagle" & incl_TE_df$V10 %in% c("AY180916_AY180916_AY180916_AY180916","FBgn0005384_FBgn0005384_3S18_3S18")|
                           incl_TE_df$V4 =="FBgn0005384_FBgn0005384_3S18_3S18" & incl_TE_df$V10 %in% c("FBgn0001207_FBgn0001207_TEHMSbeagle_HMSbeagle","AY180916_AY180916_AY180916_AY180916"))

TE_andothers_df = incl_TE_df[-interTE_df_indices,]
TEself_df = incl_TE_df[interTE_df_indices,]


inter_excl_TEs_lst = c(excl_TE_df$V4, excl_TE_df$V10)
exclTE_inter_table = as.data.frame(table(inter_excl_TEs_lst))
colnames(exclTE_inter_table) = c("RowNames", "exclTE_inter_int")

inter_incl_TEs_lst = c(TE_andothers_df$V4, TE_andothers_df$V10)
inclTE_inter_table = as.data.frame(table(inter_incl_TEs_lst))
colnames(inclTE_inter_table) = c("RowNames", "TE_and_others_int")

TEself_lst = c(TEself_df$V4, TEself_df$V10)
TEself_table = as.data.frame(table(TEself_lst))


U3_interaction_df = total_inter_full_df[c(grep("FBgn0065046", total_inter_full_df$V4),
                                                grep("FBgn0065047", total_inter_full_df$V4),
                                                grep("FBgn0065048", total_inter_full_df$V4),
                                                grep("FBgn0065046", total_inter_full_df$V10),
                                                grep("FBgn0065047", total_inter_full_df$V10),
                                                grep("FBgn0065048", total_inter_full_df$V10)),c("V4", "V10")]

U3_int_table = as.data.frame(table(c(U3_interaction_df$V4, U3_interaction_df$V10)))


colnames(U3_int_table) = c("RowNames", "U3_inter_interactions")


#snoRNA

snoRNA_interaction_df = total_inter_full_df[c(grep("snoRNA", total_inter_full_df$V4),
                                          grep("snoRNA", total_inter_full_df$V10)),c("V4", "V10")]

snoRNA_int_table = as.data.frame(table(c(snoRNA_interaction_df$V4, snoRNA_interaction_df$V10)))

colnames(snoRNA_int_table) = c("RowNames", "snoRNA_int")


# #Total intra-RNA interactions
inter_intra1_indices = which(total_ints1$V4 == total_ints1$V10)
inter_intra2_indices = which(total_ints1$V4 == total_ints1$V10)
inter_intra3_indices = which(total_ints1$V4 == total_ints1$V10)

total_intra_lst = c(total_ints1$V4[inter_intra1_indices],total_ints1$V10[inter_intra1_indices],
                    total_ints2$V4[inter_intra2_indices],total_ints2$V10[inter_intra2_indices],
                    total_ints3$V4[inter_intra3_indices],total_ints3$V10[inter_intra3_indices])

total_intra_table = as.data.frame(table(total_intra_lst))
colnames(total_intra_table) = c("RowNames", "total_intra_int")


#Combining all types of interactions into a dataframe

total_interaction_df = full_join(total_inter_table, total_intra_table, by="RowNames")
total_interaction_df = full_join(total_interaction_df, total_inter_intra_df, by="RowNames")
total_interaction_df = full_join(total_interaction_df, exclTE_inter_table, by="RowNames")
total_interaction_df = full_join(total_interaction_df, inclTE_inter_table, by="RowNames")
total_interaction_df = full_join(total_interaction_df, TEself_table, by="RowNames")
total_interaction_df = full_join(total_interaction_df, U3_int_table, by="RowNames")
total_interaction_df = full_join(total_interaction_df, tRNA_int_table, by="RowNames")
total_interaction_df = full_join(total_interaction_df, RNA28S_int_table, by="RowNames")
total_interaction_df = full_join(total_interaction_df, snoRNA_int_table, by="RowNames")

#  total_interaction_df[is.na(total_interaction_df)] =0
#
# cat = total_interaction_df$RowNames
#
# total_interaction_df = total_interaction_df %>%
#   mutate(
#     category = case_when(
#       grepl("miRNA", cat) ~ "miRNA",
#       grepl("protein_coding", cat) ~ "protein_coding",
#       grepl("3S18", cat)~ "3S18",
#       grepl("AY180916", cat)~ "Blood",
#       grepl("ncRNA", cat)~ "ncRNA",
#       grepl("pseudogene", cat)~ "pseudogene",
#       grepl("rRNA", cat)~ "rRNA",
#       grepl("snoRNA", cat)~ "snoRNA",
#       grepl("snRNA", cat)~ "snRNA",
#       grepl("tRNA", cat)~ "tRNA",
#       grepl("HMSbeagle", cat)~ "HMSbeagle",
#       TRUE ~ "Other"
#     ))

#write.csv(total_interaction_df, "total_interaction_df.csv")

total_interaction_df =read.csv("total_interaction_df.csv")

#######################
#Average interaction df
##########################

total_ints1 = cdsBlood@InputFiles$all$all$s1
total_ints2 = cdsBlood@InputFiles$all$all$s2
total_ints3 = cdsBlood@InputFiles$all$all$s3

#Average inter-RNA interactions
inter_ints1_indices = which(total_ints1$V4 != total_ints1$V10)
inter_ints2_indices = which(total_ints2$V4 != total_ints2$V10)
inter_ints3_indices = which(total_ints3$V4 != total_ints3$V10)

inter1 = total_ints1[inter_ints1_indices,c(4,10)]
inter1_table = as.data.frame(table(c(inter1$V4, inter1$V10)))
colnames(inter1_table) = c("RowNames", "s1_inter")
inter2 = total_ints2[inter_ints2_indices,c(4,10)]
inter2_table = as.data.frame(table(c(inter2$V4, inter2$V10)))
colnames(inter2_table) = c("RowNames", "s2_inter")

inter3 = total_ints3[inter_ints3_indices,c(4,10)]
inter3_table = as.data.frame(table(c(inter3$V4, inter3$V10)))
colnames(inter3_table) = c("RowNames", "s3_inter")


inter_per_sample_df = full_join(inter1_table, inter2_table, by="RowNames")

inter_per_sample_df = full_join(inter_per_sample_df, inter3_table, by="RowNames")
inter_per_sample_df[is.na(inter_per_sample_df)]=0

#Average intra-RNA interactions
intra_ints1_indices = which(total_ints1$V4 == total_ints1$V10)
intra_ints2_indices = which(total_ints2$V4 == total_ints2$V10)
intra_ints3_indices = which(total_ints3$V4 == total_ints3$V10)

intra1 = total_ints1[intra_ints1_indices,c(4,10)]
intra1_table = as.data.frame(table(c(intra1$V4, intra1$V10)))
colnames(intra1_table) = c("RowNames", "s1_intra")
intra2 = total_ints2[intra_ints2_indices,c(4,10)]
intra2_table = as.data.frame(table(c(intra2$V4, intra2$V10)))
colnames(intra2_table) = c("RowNames", "s2_intra")

intra3 = total_ints3[intra_ints3_indices,c(4,10)]
intra3_table = as.data.frame(table(c(intra3$V4, intra3$V10)))
colnames(intra3_table) = c("RowNames", "s3_intra")


intra_per_sample_df = full_join(intra1_table, intra2_table, by="RowNames")

intra_per_sample_df = full_join(intra_per_sample_df, intra3_table, by="RowNames")
intra_per_sample_df[is.na(intra_per_sample_df)]=0





##########################
#### Connection Score ####
##########################
##### Connection score filtering 
#according to Cai e.a. 2020 a connection score is computed by dividing the total number of reads connecting 
#two different RNA fragments to the coverage of chimeric reads at these two RNA fragments (coverage A_B/√(coverage Axcoverage B))

connectionScore <- function(chim_AB, chimA, chimB) {
  log_chim_AB <- log(chim_AB)
  log_chimA <- log(chimA)
  log_chimB <- log(chimB)
  
  log_score <- log_chim_AB - 0.5 * (log_chimA + log_chimB)
  return(exp(log_score))
}







#Blood
blood_tot = total_interaction_df$total_inter_int[1]
blood_int_totalint = left_join(blood_int_df, total_interaction_df, by="RowNames")
head(blood_cs_df)
blood_cs_df = blood_int_totalint %>%
  mutate(blood_cs = connectionScore(blood_int, rep(blood_tot, nrow(.)), total_inter_int))

nrow(tS18_cs_df[tS18_cs_df$tS18_cs >0.002,])
+17 + 2


ggplot(blood_cs_df, aes(x=blood_cs)) + 
  geom_histogram(binwidth = 0.001)

#beagle
sum(beagle_int_df$beagle_int[-2])
beagle_tot = total_interaction_df$total_inter_int[2]
beagle_int_totalint = left_join(beagle_int_df, total_interaction_df, by="RowNames")

beagle_cs_df = beagle_int_totalint %>%
  mutate(beagle_cs = connectionScore(beagle_int, rep(beagle_tot, nrow(.)), total_inter_int))

#tS18
tS18_tot = total_interaction_df$total_inter_int[3]
tS18_int_totalint = left_join(tS18_int_df, total_interaction_df, by="RowNames")

tS18_cs_df = tS18_int_totalint %>%
  mutate(tS18_cs = connectionScore(tS18_int, rep(tS18_tot, nrow(.)), total_inter_int))

#Getting the total number of inter-RNA interactions for each transcript

#Two types of connection scores - using only inter-RNA or also intra+inter total chimera


######################
#### MC algorithm ####
######################
# 
# ### MC algorithm requires 3 components: 1. Observed interactions 2. Transcript abdundance 3. Number to be drawing
# 
# # Important: what do we pick as transcript abundance? Options: inter, intra, inter/intra, RNAseq not appropriate, singlets
# 
# #Plotting the three optionshttp://127.0.0.1:37295/graphics/c43ae70b-f7d2-46ce-a8b5-58f5fed213dc.png
# blood_obs_tot_df = left_join(blood_int_df, total_interaction_df, by="RowNames")
# ggplot(blood_obs_tot_df ,aes(x=log2(total_inter_int), y=log2(blood_int), colour=category)) +
#   geom_point()+
#   theme(legend.position = "none") +
#   ggplot(blood_obs_tot_df, aes(x=log2(total_intra_int), y=log2(blood_int), colour=category))+
#   geom_point() + 
#   theme(legend.position = "none") +
#   ggplot(blood_obs_tot_df, aes(x=log2(total_interintra_int), y=log2(blood_int), colour=category))+
#   geom_point()
# 
# 
# blood_obs_tot_d
# 
# ggplot(blood_obs_tot_df ,aes(x=category, y=total_inter_int)) +
#   geom_col()+
#   theme(legend.position = "none") +
#   ggplot(blood_obs_tot_df ,aes(x=category, y=total_intra_int)) +
#   geom_col()+
#   theme(legend.position = "none") +
#   ggplot(blood_obs_tot_df ,aes(x=category, y=total_interintra_int)) +
#   geom_col()
#   # ggplot(blood_obs_tot_df, aes(x=log2(total_intra_int), y=log2(blood_int), colour=category))+
#   # geom_point() + 
#   # theme(legend.position = "none") +
#   # ggplot(blood_obs_tot_df, aes(x=log2(total_interintra_int), y=log2(blood_int), colour=category))+
#   # geom_point()
# 
# 
# 
# #1. Observed interactions
# #2. 
# 
# 
# 
# 
# #requires 3 components
# 
# #1: the number of inter-RNA interactions we have found in the dataset

head(total_interaction_df)
sum(beagle_int_df$transcripts$[-2])
blood_real_inter_mc = blood_int_df[2:nrow(blood_int_df),]
beagle_real_inter_mc = beagle_int_df[c(1,3:nrow(beagle_int_df)),]
tS18_real_inter_mc = tS18_int_df[c(1,2,4:nrow(tS18_int_df)),]

sum
# 
#write.csv(blood_real_inter_mc,"real_inter_blood.csv")
 write.csv(beagle_real_inter_mc,"real_inter_beagle.csv")
# write.csv(tS18_real_inter_mc,"real_inter_tS18.csv")
# 
# colnames(blood_real_inter_mc)
# head(blood_real_inter_mc)
# 
# #2: The total distribution of available for interactions
# total_interaction_df
# 
# blood_total_distr = total_interaction_df[2:nrow(total_interaction_df),]
# blood_total_distr = blood_total_distr[blood_total_distr$RowNames %in% blood_real_inter_mc$RowNames,]
# #blood_total_distr = rep(blood_total_distr$total_inter_intras_lst, blood_total_distr$Freq)
# 
# beagle_total_distr = total_interaction_df[c(1,3:nrow(total_interaction_df)),]
# beagle_total_distr = beagle_total_distr[beagle_total_distr$RowNames %in% beagle_real_inter_mc$RowNames,]
# 
# #beagle_total_distr = rep(beagle_total_distr$total_inter_intras_lst, beagle_total_distr$Freq)
# 
# nrow(total_interaction_df)
# tS18_total_distr = total_interaction_df[c(1,2,4:nrow(total_interaction_df)),]
# tS18_total_distr = tS18_total_distr[tS18_total_distr$RowNames %in% tS18_real_inter_mc$RowNames,]
# 
# #tS18_total_distr = rep(tS18_total_distr$total_inter_intras_lst, tS18_total_distr$Freq)
# 
# 
# head(tS18_total_distr)
# 
# write.csv(blood_total_distr,"blood_total_distr.csv")
# write.csv(beagle_total_distr,"beagle_total_distr.csv")
# write.csv(tS18_total_distr,"tS18_total_distr.csv")
# 
# ?left_join
# 
# #3: The number of times to randomly pull from the total distribution
# head(tS18_real_inter_mc)
sum(blood_real_inter_mc$blood_int)
sum(beagle_real_inter_mc$beagle_int)
# #beagle_real_inter_mc[grep(".2",beagle_real_inter_mc$beagle_sum),]
sum(tS18_real_inter_mc$tS18_int)
# 
# 

########

# 
# 
# #N= 100000
# 
# total_intras_df_blood = total_intras_table[total_intras_table$total_intras_lst %in% rownames(sim_higher_than_observed),]
# total_intras_df_beagle = total_intras_table[total_intras_table$total_intras_lst %in% rownames(sim_higher_than_observed_beagle),]
# total_intras_df_ts18 = total_intras_table[total_intras_table$total_intras_lst %in% rownames(sim_higher_than_observed_ts18),]
# 
# 
# write.csv(total_intras_df_blood,"total_distr_blood.csv")
# write.csv(total_intras_df_beagle,"total_distr_beagle.csv")
# write.csv(total_intras_df_ts18,"total_distr_ts18.csv")
# 
# 
# 
# 
# connection_score_beagle$RowNames = rownames(connection_score_beagle)
# connection_score_beagle = connection_score_beagle[order(connection_score_beagle$RowNames),]
# connection_score_blood$RowNames = rownames(connection_score_blood)
# connection_score_blood = connection_score_blood[order(connection_score_blood$RowNames),]
# connection_score_tS18$RowNames = rownames(connection_score_tS18)
# connection_score_tS18 = connection_score_tS18[order(connection_score_tS18$RowNames),]
# 
# nrow()
# 
# 
# 
# # write.csv(connection_score_blood,"Connection_score_blood.csv")
# # write.csv(connection_score_beagle,"Connection_score_beagle.csv")
# # write.csv(connection_score_tS18,"Connection_score_tS18.csv")


###### Importing the simulation files######

TE_type = "beagle"
TE_type="blood"
import_Sim_fun = function(TE_type, cs_df, type="normal"){
  if(type=="normal"){
    sim0 = read.csv(paste("002COMRADES_project/output/HPC/MonteCarloSim/",TE_type,"0output_MC_simulation.csv",sep=""), row.names = 1)
    sim1 = read.csv(paste("002COMRADES_project/output/HPC/MonteCarloSim/",TE_type,"1output_MC_simulation.csv",sep=""), row.names = 1)
    sim2 = read.csv(paste("002COMRADES_project/output/HPC/MonteCarloSim/",TE_type,"2output_MC_simulation.csv",sep=""), row.names = 1)
    sim3 = read.csv(paste("002COMRADES_project/output/HPC/MonteCarloSim/",TE_type,"3output_MC_simulation.csv",sep=""), row.names = 1)
    sim4 = read.csv(paste("002COMRADES_project/output/HPC/MonteCarloSim/",TE_type,"4output_MC_simulation.csv",sep=""), row.names = 1)
    sim5 = read.csv(paste("002COMRADES_project/output/HPC/MonteCarloSim/",TE_type,"5output_MC_simulation.csv",sep=""), row.names = 1)
    sim6 = read.csv(paste("002COMRADES_project/output/HPC/MonteCarloSim/",TE_type,"6output_MC_simulation.csv",sep=""), row.names = 1)
    sim7 = read.csv(paste("002COMRADES_project/output/HPC/MonteCarloSim/",TE_type,"7output_MC_simulation.csv",sep=""), row.names = 1)
    sim8 = read.csv(paste("002COMRADES_project/output/HPC/MonteCarloSim/",TE_type,"8output_MC_simulation.csv",sep=""), row.names = 1)
    sim9 = read.csv(paste("002COMRADES_project/output/HPC/MonteCarloSim/",TE_type,"9output_MC_simulation.csv",sep=""), row.names = 1)
    sim_occ = sim0$ocurrences/10000 + sim1$ocurrences/10000 + sim2$ocurrences/10000 + sim3$ocurrences/10000 + sim4$ocurrences/10000 + sim5$ocurrences/10000 + sim6$ocurrences/10000+ sim7$ocurrences/10000 + sim8$ocurrences/10000 + sim9$ocurrences/10000
    sim_occ = sim_occ/10
    sim_p = sim0$higher_pval + sim1$higher_pval + sim2$higher_pval + sim3$higher_pval + sim4$higher_pval + sim5$higher_pval + sim6$higher_pval+ sim7$higher_pval + sim8$higher_pval + sim9$higher_pval
    sim_p = sim_p/100000
    sim =data.frame(sim_pval = sim_p, sim_occ = sim_occ)
    #sim[sim$sim ==0,] =0.00001
    head(sim)
    sim$RowNames = sim0$RowNames
    #cs_sim = full_join(real_inter, sim_p, by="RowNames")
    #ownames(cs_sim) = cs_sim$RowNames
    cs_sim = left_join(cs_df,sim, by="RowNames")
    #cs_sim = cs_sim[,c(1,3,4,5)]
    return(cs_sim)
  }
  if(type=="random"){
    sim0 = read.csv(paste("002COMRADES_project/output/HPC/MonteCarloSim/",TE_type,"0output_MC_simulation_rd.csv",sep=""), row.names = 1)
    sim1 = read.csv(paste("002COMRADES_project/output/HPC/MonteCarloSim/",TE_type,"1output_MC_simulation_rd.csv",sep=""), row.names = 1)
    sim2 = read.csv(paste("002COMRADES_project/output/HPC/MonteCarloSim/",TE_type,"2output_MC_simulation_rd.csv",sep=""), row.names = 1)
    sim3 = read.csv(paste("002COMRADES_project/output/HPC/MonteCarloSim/",TE_type,"3output_MC_simulation_rd.csv",sep=""), row.names = 1)
    sim4 = read.csv(paste("002COMRADES_project/output/HPC/MonteCarloSim/",TE_type,"4output_MC_simulation_rd.csv",sep=""), row.names = 1)
    sim5 = read.csv(paste("002COMRADES_project/output/HPC/MonteCarloSim/",TE_type,"5output_MC_simulation_rd.csv",sep=""), row.names = 1)
    sim6 = read.csv(paste("002COMRADES_project/output/HPC/MonteCarloSim/",TE_type,"6output_MC_simulation_rd.csv",sep=""), row.names = 1)
    sim7 = read.csv(paste("002COMRADES_project/output/HPC/MonteCarloSim/",TE_type,"7output_MC_simulation_rd.csv",sep=""), row.names = 1)
    sim8 = read.csv(paste("002COMRADES_project/output/HPC/MonteCarloSim/",TE_type,"8output_MC_simulation_rd.csv",sep=""), row.names = 1)
    sim9 = read.csv(paste("002COMRADES_project/output/HPC/MonteCarloSim/",TE_type,"9output_MC_simulation_rd.csv",sep=""), row.names = 1)
    sim_occ = sim0$ocurrences/10000 + sim1$ocurrences/10000 + sim2$ocurrences/10000 + sim3$ocurrences/10000 + sim4$ocurrences/10000 + sim5$ocurrences/10000 + sim6$ocurrences/10000+ sim7$ocurrences/10000 + sim8$ocurrences/10000 + sim9$ocurrences/10000
    sim_occ = sim_occ/10
    sim_p = sim0$higher_pval + sim1$higher_pval + sim2$higher_pval + sim3$higher_pval + sim4$higher_pval + sim5$higher_pval + sim6$higher_pval+ sim7$higher_pval + sim8$higher_pval + sim9$higher_pval
    sim_p = sim_p/100000
    sim =data.frame(sim_pval_rd = sim_p, sim_occ_rd = sim_occ)
    sim$RowNames = sim9$RowNames
    sim = left_join(cs_df, sim, by="RowNames")
    return(sim)
  }
  # if(type=="batches"){
    
    batchfile = read.csv("002COMRADES_project/output/HPC/MonteCarloSim/batch_pvalues.csv")
    
    # sim0 = read.csv(paste("002COMRADES_project/output/HPC/MonteCarloSim/",TE_type,"0output_MC_simulation.csv",sep=""), row.names = 1)
    # sim1 = read.csv(paste("002COMRADES_project/output/HPC/MonteCarloSim/",TE_type,"1output_MC_simulation.csv",sep=""), row.names = 1)
    # sim2 = read.csv(paste("002COMRADES_project/output/HPC/MonteCarloSim/",TE_type,"2output_MC_simulation.csv",sep=""), row.names = 1)
    # sim3 = read.csv(paste("002COMRADES_project/output/HPC/MonteCarloSim/",TE_type,"3output_MC_simulation.csv",sep=""), row.names = 1)
    # sim4 = read.csv(paste("002COMRADES_project/output/HPC/MonteCarloSim/",TE_type,"4output_MC_simulation.csv",sep=""), row.names = 1)
    # sim5 = read.csv(paste("002COMRADES_project/output/HPC/MonteCarloSim/",TE_type,"5output_MC_simulation.csv",sep=""), row.names = 1)
    # sim6 = read.csv(paste("002COMRADES_project/output/HPC/MonteCarloSim/",TE_type,"6output_MC_simulation.csv",sep=""), row.names = 1)
    # sim7 = read.csv(paste("002COMRADES_project/output/HPC/MonteCarloSim/",TE_type,"7output_MC_simulation.csv",sep=""), row.names = 1)
    # sim8 = read.csv(paste("002COMRADES_project/output/HPC/MonteCarloSim/",TE_type,"8output_MC_simulation.csv",sep=""), row.names = 1)
    # sim9 = read.csv(paste("002COMRADES_project/output/HPC/MonteCarloSim/",TE_type,"9output_MC_simulation.csv",sep=""), row.names = 1)
    # sim_p0 = sim0$higher_pval/10000
    # sim_p1 = sim1$higher_pval/10000
    # sim_p2 = sim2$higher_pval/10000
    # sim_p3 = sim3$higher_pval/10000
    # sim_p4 = sim4$higher_pval/10000
    # sim_p5 = sim5$higher_pval/10000
    # sim_p6 = sim6$higher_pval/10000
    # sim_p7 = sim7$higher_pval/10000
    # sim_p8 =  sim8$higher_pval/10000
    # sim_p9= sim9$higher_pval/10000
    # 
    # sim_occ0 = sim0$ocurrences/10000
    # sim_occ1 = sim1$ocurrences/10000
    # sim_occ2 = sim2$ocurrences/10000
    # sim_occ3 = sim3$ocurrences/10000
    # sim_occ4 = sim4$ocurrences/10000
    # sim_occ5 = sim5$ocurrences/10000
    # sim_occ6 = sim6$ocurrences/10000
    # sim_occ7 = sim7$ocurrences/10000
    # sim_occ8 =  sim8$ocurrences/10000
    # sim_occ9= sim9$ocurrences/10000
    # sim = data.frame(sim_p0, sim_p1, sim_p2, sim_p3, sim_p4, sim_p5, sim_p6, sim_p7, sim_p8, sim_p9,
    #                  sim_occ0, sim_occ1, sim_occ2, sim_occ3, sim_occ4, sim_occ5, sim_occ6, sim_occ7, sim_occ8, sim_occ9)
    # sim$RowNames = sim0$RowNames
    return(sim)
  }
  
  

blood_cs_sim = import_Sim_fun("blood", blood_cs_df)
#blood_cs_sim$observed_rd = blood_int_rd$Freq

beagle_cs_sim = import_Sim_fun("beagle", beagle_cs_df)

sum(beagle_cs_sim$sim_occ[-2])
sum(blood_cs_sim$sim_occ[-1])
sum(tS18_cs_sim$sim_occ[-3])

tS18_cs_sim = import_Sim_fun("tS18",tS18_cs_df)


blood_int_rd = read.csv("002COMRADES_project/output/HPC/MonteCarloSim/bloodsingle_sim.csv", row.names = 1)
beagle_int_rd = read.csv("002COMRADES_project/output/HPC/MonteCarloSim/beaglesingle_sim.csv", row.names = 1)
#colnames(beagle_int_rd) = c("RowNames", "Beagle_rd_obs")
tS18_int_rd = read.csv("002COMRADES_project/output/HPC/MonteCarloSim/tS18single_sim.csv", row.names = 1)
#colnames(tS18_int_rd) = c("RowNames", "tS18_rd_obs")

beagle_cs_sim = left_join(beagle_cs_sim, beagle_int_rd, by="RowNames")
tS18_cs_sim = left_join(tS18_cs_sim, tS18_int_rd, by="RowNames")

blood_cs_sim_2 = import_Sim_fun("blood", blood_int_rd, "random")
blood_cs_sim_2$observed_rd = blood_int_rd$Freq

beagle_cs_sim_2 = import_Sim_fun("beagle", beagle_int_rd, "random")
beagle_cs_sim_2$observed_rd = beagle_int_rd$Freq

tS18_cs_sim_2 = import_Sim_fun("tS18",tS18_int_rd, "random")
tS18_cs_sim_2$observed_rd = tS18_int_rd$Freq

blood_cs_sim_3 = import_Sim_fun("blood", blood_cs_df, "batches")
beagle_cs_sim_3 = import_Sim_fun("beagle", beagle_cs_df, "batches")
tS18_cs_sim_3 = import_Sim_fun("tS18", tS18_cs_df, "batches")

blood_cs_sim_3= left_join(blood_cs_sim_3, blood_int_rd, by="RowNames")


##### Adjusted p values #####

#Using FDR benjamini hochberg

blood_cs_sim$q_val_blood = p.adjust(blood_cs_sim$sim_pval, method="BH")


beagle_cs_sim$q_val_beagle = p.adjust(beagle_cs_sim$sim_pval, method="BH")
tS18_cs_sim$q_val_tS18 = p.adjust(tS18_cs_sim$sim_pval, method="BH")


# 
# ggplot(blood_cs_sim, aes(x=-log2(p), y=log2(blood_cs), colour=category)) + 
#   geom_point() + 
#   geom_vline(xintercept = -log2(0.01)) + 
#   geom_hline(yintercept = log2(0.001))
# 




#Change in p value progression
batchfile = read.csv("002COMRADES_project/output/HPC/MonteCarloSim/blood_pbatch.csv")

blood_pbatches_vector = c(batchfile[,2],0)
blood_pbatches_vector = as.integer(blood_pbatches_vector)
n = length(blood_pbatches_vector)

blood_pbatches_mat <- matrix(blood_pbatches_vector, nrow = n/100, ncol = 100, byrow = FALSE)
blood_pbatches_mat = blood_pbatches_mat/1
blood_pbatches_df = as.data.frame(blood_pbatches_mat)

blood_cs_sim_3 = blood_cs_sim_3[1:14973,]
tail(blood_pval_adj)
adj_pval_batch_df = data.frame("batch1"= rep(0,nrow(blood_cs_sim_3)),
                               "batch2"= rep(0,nrow(blood_cs_sim_3)),
                               "batch3"= rep(0,nrow(blood_cs_sim_3)),
                               "batch4"= rep(0,nrow(blood_cs_sim_3)),
                               "batch5"= rep(0,nrow(blood_cs_sim_3)),
                               "batch6"= rep(0,nrow(blood_cs_sim_3)),
                               "batch7"= rep(0,nrow(blood_cs_sim_3)),
                               "batch8"= rep(0,nrow(blood_cs_sim_3)),
                               "batch9"= rep(0,nrow(blood_cs_sim_3)),
                               "batch10"= rep(0,nrow(blood_cs_sim_3)))
blood_cs_sim_3[is.na(blood_cs_sim_3)] = 0

for(j in 1:10){
  current_occbatch = blood_cs_sim_3[,j+10]
  current_pvalbatch = blood_cs_sim_3[,j]
  for(i in 1:nrow(blood_cs_sim_3)){
    current_mean= (current_occbatch[i] + blood_pval_adj$blood_int[i])/2
    current_pval = current_pvalbatch[i]
    #current_diff = [i] - blood_pval_adj$blood_int[i]
    #Choose closest mean
    closest_mean_rd = which.min(abs(current_mean - blood_pval_adj$mean_rd))
    
    #Take 500 window - 250 below and 250 above
    begin_window = closest_mean_rd-250
    end_window = closest_mean_rd+250
    if(begin_window < 0){
      end_window = end_window - begin_window
      begin_window = 0
    }
    if(end_window > nrow(blood_pval_adj)){
      begin_window = begin_window - (end_window - nrow(blood_pval_adj))
      end_window = nrow(blood_pval_adj)
    }
    local_distr = blood_pval_adj$sim_pval_rd[begin_window:end_window]
    #local_distr = blood_pval_adj$diff_rd[begin_window:end_window]
    
    local_distr_sd = sd(local_distr)
    local_distr_mean = mean(local_distr)
    
    adj_pval = pnorm(current_pval, mean = local_distr_mean, local_distr_sd, lower.tail = TRUE)
    #adj_pval = pnorm(current_diff, mean = local_distr_mean, local_distr_sd, lower.tail = TRUE)
    
    adj_pval_lst = c(adj_pval_lst, adj_pval)
    #local_sd_lst = c(local_sd_lst, local_distr_sd)
  }
  adj_pval_batch_df[,j] = adj_pval_lst
  adj_pval_lst = c()
}


pval_logfc = 
  adj_pval_batch_df2 = data.frame("category" = c(rep("batch1", 14973), 
                                                 rep("batch2", 14973),
                                                 rep("batch3", 14973),
                                                 rep("batch4", 14973),
                                                 rep("batch5", 14973),
                                                 rep("batch6", 14973),
                                                 rep("batch7", 14973),
                                                 rep("batch8", 14973),
                                                 rep("batch9", 14973),
                                                 rep("batch10", 14973)),
                                  "adj_pval" = c(adj_pval_batch_df$batch1,adj_pval_batch_df$batch2,adj_pval_batch_df$batch3,adj_pval_batch_df$batch4,adj_pval_batch_df$batch5,adj_pval_batch_df$batch6,adj_pval_batch_df$batch7,adj_pval_batch_df$batch8,adj_pval_batch_df$batch9,adj_pval_batch_df$batch10))

previous_pbatch = adj_pval_batch_df[,1]
previous_occbatch = blood_cs_sim_3[,11]
log_fc_plist = c()
log_fc_occlist = c()

for(batch in 2:10){
  current_pbatch = adj_pval_batch_df[,batch]
  combined_pbatch = (previous_pbatch+current_pbatch)/2
  log_fc_p = log2(combined_pbatch/current_pbatch)
  
  #observed over simulated data
  current_occbatch = blood_cs_sim_3[,batch+10]
  combined_occbatch = (previous_occbatch+current_occbatch)/2
  log_fc_occ = log2(combined_occbatch/current_occbatch)
  
  #print(log_fc)
  previous_pbatch = combined_pbatch
  previous_occbatch = combined_occbatch
  log_fc_plist = c(log_fc_plist, log_fc_p)
  log_fc_occlist = c(log_fc_occlist, log_fc_occ)
  
}


log_fc_df = data.frame(category = as.factor(rep(1:9,each= 14973)),
                       log2_fc = log_fc_plist) %>%
  filter(log2_fc >= quantile(log2_fc, 0.25) - 1.5 * IQR(log2_fc) & 
           log2_fc <= quantile(log2_fc, 0.75) + 1.5 * IQR(log2_fc))

min(log_fc_df$log2_fc)

ggplot(log_fc_df, aes(category,log2_fc))+
  geom_boxplot(outlier.shape = NA) +
  labs(title = "Blood p-values progression in Monte Carlo simulation",
       x="batch",
       y="log2 Fold Change")+
  theme_minimal()






###############################
##### Investigating properties of interactions#########
#######################################################
#Plotting the replicates to test reproducibility

#rep1 vs rep2


ggplot(inter_per_sample_df, aes(x=log2(s1_inter+1), y=log2(s2_inter+1)))+
  geom_point(colour = "#C191ED", alpha=0.3) + 
  geom_point(data = intra_per_sample_df, aes(x=log2(s1_intra+1), y=log2(s2_intra+1)), 
             colour = "#BEED91", alpha = 0.3 ) + 
  theme_minimal() + 
  xlab("log2(freq RNA interactions), Rep 1") + 
  ylab("log2(freq RNA interactions), Rep 2")
cor(inter_per_sample_df$s1_inter, inter_per_sample_df$s2_inter, method="pearson")
?cor
cor
ggsave("s1s2_reproducibility.tiff", height=4, width=4)

ggplot(inter_per_sample_df, aes(x=log2(s1_inter+1), y=log2(s3_inter+1)))+
  geom_point(colour = "#C191ED", alpha=0.3) + 
  geom_point(data = intra_per_sample_df, aes(x=log2(s1_intra+1), y=log2(s3_intra+1)), 
             colour = "#BEED91", alpha = 0.3 ) + 
  theme_minimal() + 
  xlab("log2(freq RNA interactions), Rep 1") + 
  ylab("log2(freq RNA interactions), Rep 3")
cor(inter_per_sample_df$s1_inter, inter_per_sample_df$s3_inter, method="pearson")
ggsave("s1s3_reproducibility.tiff", height=4, width=4)

ggplot(inter_per_sample_df, aes(x=log2(s2_inter+1), y=log2(s3_inter+1)))+
  geom_point(colour = "#C191ED", alpha=0.1) + 
  geom_point(data = intra_per_sample_df, aes(x=log2(s2_intra+1), y=log2(s3_intra+1)), 
             colour = "#BEED91", alpha = 0.3 ) + 
  theme_minimal() + 
  xlab("log2(freq RNA interactions), Rep 2") +
  ylab("log2(freq RNA interactions), Rep 3")
cor(inter_per_sample_df$s2_inter, inter_per_sample_df$s3_inter, method="pearson")
ggsave("s2s3_reproducibility.tiff", height=4, width=4)


#Colours per category
colours_lst = c("Blood" = "#b04d48", "HMSbeagle" = "#e0aa24", "3S18" = "#2a4864","protein_coding" = "lightgrey","rRNA"= "#8dbd05","ncRNA"= "#00a1ae","snRNA"= "#5e36cc", "snoRNA"="#fe318e", "miRNA"="#d9b6db", "pseudogene"="#fb5607", "tRNA" = "darkgreen")



#Plotting inter vs intra chimeric reads

total_interaction_df_2 = data.frame(type = c(rep("intra-RNA", nrow(total_interaction_df)), rep("inter-RNA", nrow(total_interaction_df))),
                                    chimeric_reads = c(total_interaction_df$total_intra_int, total_interaction_df$total_inter_int),
                                    category = rep(total_interaction_df$category, 2))

total_interaction_df_2$RowNames = total_interaction_df$RowNames
#Printing total inter and intra RNA interactions

sum(total_interaction_df_2$chimeric_reads[total_interaction_df_2$type =="inter-RNA"])/2
sum(total_interaction_df_2$chimeric_reads[total_interaction_df_2$type =="intra-RNA"])/2

ggplot(total_interaction_df_2, aes(x=type, y=chimeric_reads/2000000, fill=category)) +
  geom_bar(position="dodge", stat="identity") + 
  scale_fill_manual(values = colours_lst) + 
  xlab("Type of interaction") + 
  ylab("Number of chimeric reads * 1e+06") + 
  theme_minimal() + 
  theme(legend.direction = "horizontal", legend.position = "bottom")

ggsave("002COMRADES_project/output/interintra_barplot.tiff", height = 4.5, width = 5.3)

#Plotting RNAseq data from ... to show original distribution

aubago_tots_RNAS = as.data.frame(read.csv("002COMRADES_project/output/quant.sf", sep="\t"))
colnames(aubago_tots_RNAS) = c("RowNames", colnames(aubago_tots_RNAS)[-1])
combined_aubagotots_comradesint = left_join(total_interaction_df_2, aubago_tots_RNAS, by="RowNames")

nrow(combined_aubagotots_comradesint[combined_aubagotots_comradesint$type =="intra-RNA",])
nrow(combined_aubagotots_comradesint[combined_aubagotots_comradesint$type =="inter-RNA",])


combined_aubagotots_comradesint = rbind(combined_aubagotots_comradesint, combined_aubagotots_comradesint[combined_aubagotots_comradesint$type =="inter-RNA",])

nrow(combined_aubagotots_comradesint[combined_aubagotots_comradesint$type =="intra-RNA",])
#Creating seperate ones of chimeric reads containing the TEs and ones not containing the chimeric reads
nrow(combined_aubagotots_comradesint)
# Part 1
logfc1 <- log2((combined_aubagotots_comradesint$chimeric_reads[1:22101] + 0.00001)/
                 (combined_aubagotots_comradesint$TPM[1:22101] + 0.00001))

# Part 2
logfc2 <- log2((total_interaction_df$exclTE_inter_int + 0.00001)/
                 (combined_aubagotots_comradesint$TPM[22102:44202] + 0.00001))

# Part 3
logfc3 <- log2((total_interaction_df$TE_and_others_int + 0.00001)/
                 (combined_aubagotots_comradesint$TPM[44203:66303] + 0.00001))

# Combine and assign
combined_aubagotots_comradesint$logFC <- c(logfc1, logfc2, logfc3)

combined_aubagotots_comradesint$logFC_type = c(rep("intra-RNA",22101), rep("inter-RNA excl TEs", 22101), rep("inter-RNA incl TEs", 22101))
combined_aubagotots_comradesint = combined_aubagotots_comradesint[-c(22102, 22103, 22104),]

combined_aubagotots_comradesint = na.omit(combined_aubagotots_comradesint)
#combined_aubagotots_comradesint$logFC_type[1:3] = c("intra-RNA_TE")

mean(combined_logfc_df$mean_logFC)
combined_logfc_df =  combined_aubagotots_comradesint %>%
  group_by(logFC_type, category) %>%
  summarise(mean_logFC = mean(logFC), .groups = "drop")
print(combined_logfc_df, n=30)
ggplot(combined_logfc_df, aes(x=logFC_type, y=mean_logFC/2, fill=category)) +
  geom_col(position = "dodge") + 
  scale_fill_manual(values = colours_lst) + 
  xlab("Type of interaction") + 
  ylab("Log2FC of Chimeric over RNAseq reads") + 
  theme_minimal() + 
  theme(legend.direction = "horizontal", legend.position = "bottom")
unique(combined_logfc_df$category)

ggsave("002COMRADES_project/output/pulled_downTE_int.tiff", heigh=4.5, width = 5.3)

#separating out pulled down TEs

sepTEs_interaction_df = data.frame(type = c(rep("inter- and with other TE interactions", nrow(total_interaction_df)), rep("non-TE RNA interactions", nrow(total_interaction_df))),
                                   chimeric_reads = c(total_interaction_df$TEself_int[1:3],total_interaction_df$inclTE_inter_int[4:nrow(total_interaction_df)], total_interaction_df$exclTE_inter_int),
                                   category = c(rep(total_interaction_df$category, 2)))

length(c(total_interaction_df$TEself_int[1:3],total_interaction_df$inclTE_inter_int[4:nrow(total_interaction_df)], total_interaction_df$exclTE_inter_int))
ggplot(sepTEs_interaction_df, aes(x=type, y=chimeric_reads/1000000, fill=category)) +
  geom_bar(position="dodge", stat="identity") + 
  labs(title = "Inter-RNA interactions including and excluding pulled-down TEs") + 
  scale_fill_manual(values = colours_lst) + 
theme_minimal()




###############################
##### Top 1000 inter-RNAs #####
###############################


#Top 1000 of all inter-RNA transcripts

total_interaction_pvalues = right_join(total_interaction_df[,c(2,3,5)], blood_cs_sim[,c("RowNames","q_val_blood")], by="RowNames")
total_interaction_pvalues = full_join(total_interaction_pvalues, beagle_cs_sim[,c("RowNames","q_val_beagle")], by="RowNames")
total_interaction_pvalues = full_join(total_interaction_pvalues, tS18_cs_sim[,c("RowNames","q_val_tS18")], by="RowNames")
total_interaction_pvalues[is.na(total_interaction_pvalues)] = 1
total_interactions_top1000 = total_interaction_pvalues[order(-total_interaction_pvalues$total_interintra_int),][1:1000,]
#total_interactions_top1000[total_interactions_top1000$q_val_blood<0.01,]


#Random inter-RNA interactions across all interactions

total_distr = c(rep(total_interaction_df$RowNames, total_interaction_df$total_interintra_int))
total_interRNA_interactions = sum(total_interaction_df$total_inter_int)
total_interaction_rd_df = data.frame(table(as.character(sample(total_distr, total_interRNA_interactions))))
colnames(total_interaction_rd_df) = c("RowNames", "total_interRNA_rd_sample")
#Random interactions

total_interaction_rd_top1000 = left_join(total_interactions_top1000,total_interaction_rd_df, by="RowNames")

total_interaction_rd_top1000[total_interaction_rd_top1000$q_val_blood <0.05,]



total_interaction_rd_top1000 <- total_interaction_rd_top1000 %>%
  mutate(category = case_when(
    #rowSums(cbind(q_val_blood < 0.05, q_val_beagle < 0.05, q_val_tS18 < 0.05)) >=2 ~ "2 TEs: p<0.05",
    #rowSums(cbind(q_val_blood < 0.05, q_val_beagle < 0.05, q_val_tS18 < 0.05)) ==3 ~ "At least 3 TEs: p<0.05",
    #rowSums(cbind(q_val_blood < 0.05, q_val_beagle < 0.05, q_val_tS18 > 0.05)) >=2 ~ "Blood and HMSbeagle: p<0.05",
    #rowSums(cbind(q_val_blood < 0.05, q_val_tS18 < 0.05, q_val_beagle > 0.05)) >=2 ~ "Blood and 3S18: p<0.05",
    #rowSums(cbind(q_val_beagle < 0.05, q_val_tS18 < 0.05, q_val_blood < 0.05)) >=2 ~ "Blood and HMSbeagle: p<0.05",
    
    q_val_blood < 0.05 & q_val_beagle < 0.05 & q_val_tS18 >= 0.05~ "Blood and HMSbeagle: p<0.05",
    q_val_blood < 0.05 & q_val_tS18 < 0.05 & q_val_beagle >= 0.05~ "Blood and 3S18: p<0.05",
    q_val_beagle < 0.05 & q_val_tS18 < 0.05 & q_val_blood >= 0.05~ "3S18 and HMSbeagle: p<0.05",
    
    q_val_blood < 0.05 &q_val_beagle>=0.05 & q_val_tS18 >=0.05~ "Blood: p<0.05",
    q_val_beagle < 0.05 &q_val_blood>=0.05 & q_val_tS18 >=0.05 ~ "HMSbeagle: p<0.05",
    q_val_tS18 < 0.05 &q_val_blood >=0.05 &q_val_beagle>=0.05 ~ "3S18: p<0.05",
    q_val_blood >=0.05 &q_val_beagle>=0.05 & q_val_tS18 >=0.05 ~"insignificant",

  ))


table(total_interaction_rd_top1000$category)
ggplot(total_interaction_rd_top1000) + 
  geom_point(aes(x = log2(total_interintra_int+1
                          ), y=log2(total_interRNA_rd_sample+1), colour="random", alpha="random")) + 
  geom_point(aes(x=log2(total_interintra_int),
                 y=log2(total_inter_int),
                 colour = category,
                 alpha=category,
             shape=category))+ 
  geom_point(data = total_interaction_rd_top1000[total_interaction_rd_top1000$category!="insignificant",],aes(x=log2(total_interintra_int),
                 y=log2(total_inter_int),
                 colour = category,
                 alpha=category,
                 shape=category))+ 
  geom_point(data = total_interaction_rd_top1000[total_interaction_rd_top1000$category!="insignificant" & total_interaction_rd_top1000$category!="2 TEs: p<0.05" ,],aes(x=log2(total_interintra_int),
                                                                                                              y=log2(total_inter_int),
                                                                                                              colour = category,
                                                                                                              alpha=category,
                                                                                                              shape=category))+ 
  
  scale_alpha_manual(values = c("random" = 0.5,
                                "insignificant" = 0.1,
                                  "2 TEs: p<0.05" = 0.8,
                                "At least 3 TEs: p<0.05" = 1,
                                 "Blood: p<0.05" = 1,
                                "HMSbeagle: p<0.05" = 1,
                                "3S18: p<0.05" = 1)) +
  
  scale_colour_manual(values = c("random" = "lightgrey",
                                 "insignificant" = "#afc3db",
                                 "2 TEs: p<0.05" = "#606c38",
                                 "At least 3 TEs: p<0.05" = "#606c38",
                                 "Blood: p<0.05" = "#b04d48",
                                 "HMSbeagle: p<0.05" = "#e0aa24",
                                 "3S18: p<0.05" = "#2a4864")) +
  scale_shape_manual(values = c("2 TEs: p<0.05" = 18,
                                "random" = 1,
                                "insignificant" = 16,
                                "At least 3 TEs: p<0.05" = 16,
                                "Blood: p<0.05" = 16,
                                "HMSbeagle: p<0.05" = 16,
                                "3S18: p<0.05" = 16)) + 
  labs(title = "Top 1000 inter-RNA interactions",
      x = "log2(pseudo counts)",
      y = "log2(inter-RNA interactions)") + 
  
  theme_minimal()

ggsave("002COMRADES_project/output/top1000.tiff", height=5.5, width=6.8)
#Verification of inter-RNA interactions using snoRNA



#Interaction with snoRNAs - verificatino of inter-RNA interaction
total_interaction_df[is.na(total_interaction_df)] = 0
table(total_interaction_df$category[total_interaction_df$U3_inter_interactions >0])
ggplot(total_interaction_df, aes(x = log2(total_interintra_int + 1), y = log2(snoRNA_int + 1), 
                                 colour=ifelse(category=="snoRNA", "snoRNA", "other"))) +
  geom_point(alpha =0.5) +
  scale_colour_manual(values = c("snoRNA" = "#E36414", other="#afc3db")) + 
  xlab("log2(proxy reads)")+
  ylab("log2(snoRNA interactions)") +
  theme_minimal() 

ggsave("002COMRADES_project/output/snoRNA_inter_validation.tiff", height = 4, width=11)

# between 18S and 28S rRNAs, MALAT1 and U1 (ref. 20) and between spliceosomal RNAs, respectively (q < 0.002). 



#Man-Whittney test
total_interaction_df$snoRNA_binary = ifelse(total_interaction_df$category == "snoRNA", 1, 0)
wilcox.test(snoRNA_int  ~ snoRNA_binary, data = total_interaction_df)
#Wilcoxon rank sum test with continuity correction
# data:  snoRNA_int by snoRNA_binary
# W = 14735, p-value < 2.2e-16
# alternative hypothesis: true location shift is not equal to 0

?wilcox.test

#Comparing replicates in a pairwise manner







#circos graphs
total_int_cat_table$total = rowSums(total_int_cat_table)
circos.par("track.height" = 0.1)
rownames(total_int_cat_table) = c("tS18" ,rownames(total_int_cat_table)[2:11])
circos.initialize(sectors = rownames(total_int_cat_table),xlim=as.matrix(data.frame(rep(0,11), log2(total_int_cat_table$total))))
circos.track(1)
total_int_cat_table
?circos.initialize




#########################
#### intron/exon and genomic distance #####
#########################
#?????

blood_transcript_ids <- sub("_.*", "", blood_cs_sim$RowNames)


flybase_mart <- useMart("ensembl", dataset = "dmelanogaster_gene_ensembl")

transcript_ids <- c("FBtr0070148", "FBtr0070149")  # Replace with your transcript IDs

# Query for exons
exon_data <- getBM(
  attributes = c("ensembl_transcript_id","exon_chrom_start", "exon_chrom_end"),
  filters = "flybase_transcript_id",
  values = blood_transcript_ids,
  mart = flybase_mart
)

# Query for introns (intron data may be inferred)
intron_data <- getBM(
  attributes = c("ensembl_transcript_id","ensembl_exon_id", "transcript_exon_intron"),
  filters = "flybase_transcript_id",
  values = blood_transcript_ids,
  mart = flybase_mart
)


listAttributes(mart=flybase_mart)[grep("intron",listAttributes(mart=flybase_mart)[,1]),]


############################
#### MC and CS combined ####
############################

#Plotting simulated occurrences vs random interactions with connection score and p-values 


#Observed data
ggplot() +
  geom_point(data = blood_cs_sim, aes(x = log2(sim_occ + 1), 
                               y = log2(blood_int + 1), 
                               color = ifelse(q_val_blood <= 0.05 & blood_cs > 0.002,"Blood significant","insignificant"), 
                               alpha = ifelse(q_val_blood <= 0.05 & blood_cs > 0.002,"Blood significant","insignificant")),
             size = 1) + 
  geom_point(data = beagle_cs_sim, aes(x=log2(sim_occ + 1), 
                                y=log2(beagle_int + 1),
                                color = ifelse(q_val_beagle <= 0.05 & beagle_cs > 0.002,"Beagle significant","insignificant"), 
                                alpha = ifelse(q_val_beagle <= 0.05 & beagle_cs > 0.002,"Beagle significant","insignificant")),
             size = 1) + 
  geom_point(data = tS18_cs_sim, aes(x=log2(sim_occ + 1), 
                                       y=log2(tS18_int + 1),
                                       color = ifelse(q_val_tS18 <= 0.05 & tS18_cs > 0.002,"3S18 significant","insignificant"),
                                     alpha = ifelse(q_val_tS18 <= 0.05 & tS18_cs > 0.002,"3S18 significant","insignificant")),
             size = 1) + 
  scale_color_manual("Interaction significance",
                     values = c("Blood significant" = "#b04d48",
                                "Beagle significant" = "#e0aa24",
                                "3S18 significant" = "#2a4864",
                                "insignificant" = "lightgrey"),
                    ) +
    scale_alpha_manual("Interaction significance", 
                       values = c("Blood significant" = 1,
                                  "Beagle significant" = 1,
                                  "3S18 significant" = 1,
                                  "insignificant" = 0.1)) +
  labs(x = "log2(inter-RNA chimera), simulated",
       y = "log2(inter-RNA chimera), observed") + 
  theme_minimal()


ggsave("002COMRADES_project/output/TE_obs_vs_sim.tiff", heigh = 4, width = 5.6)
#Random data



#maybe do connection score
ggplot() +
  geom_point(data = blood_cs_sim_2, aes(x = log2(sim_occ_rd+1), 
                                      y = log2(observed_rd), 
                                      color = ifelse(sim_pval_rd <= 0.05,"Blood significant","insignificant"),
                                      alpha =  ifelse(sim_pval_rd <= 0.05,"Blood significant","insignificant")),
             size = 1) + 
  geom_point(data = beagle_cs_sim_2, aes(x=log2(sim_occ_rd+1), 
                                       y=log2(observed_rd),
                                       color = ifelse(sim_pval_rd   <= 0.05,"Beagle significant","insignificant"), 
                                       alpha = ifelse(sim_pval_rd   <= 0.05,"Beagle significant","insignificant")),
             size = 1) + 
  geom_point(data = tS18_cs_sim_2, aes(x=log2(sim_occ_rd+1), 
                                     y=log2(observed_rd),
                                     color = ifelse(sim_pval_rd   <= 0.05,"3S18 significant","insignificant"), 
                                     alpha = ifelse(sim_pval_rd   <= 0.05,"3S18 significant","insignificant")),
             size = 1) + 
  scale_color_manual("Interaction significance", 
                     values = c("Blood significant" = "#b04d48",
                                "Beagle significant" = "#e0aa24",
                                "3S18 significant" = "#2a4864",
                                "insignificant" = "lightgrey")) +
  scale_alpha_manual("Interaction significance", 
                     values = c("Blood significant" = 1,
                                                            "Beagle significant" = 1,
                                                            "3S18 significant" = 1,
                                                            "insignificant" = 0.1)) +
  labs( x = "log2(inter-RNA chimera), simulated",
       y = "log2(inter-RNA chimera), random ") + 
  theme_minimal()
ggsave("002COMRADES_project/output/TE_rd_vs_sim.tiff", heigh = 4, width = 5.5)
#Random data

#####p-value stabilisation####

#Boxplot showing how the p value distribution changed as the Monte
# Carlo simulation progressed using n = 735,247 pairwise RNA co-barcoding events. The 100,000 simulations were divided
# into 20 batches of 5,000 simulations. Batches were sequentially added to one another and the log2 fold change of the p
# value for each feature before and after the addition of a new batch of simulations was calculated. The graph shows that as
# the simulation progressed the p values stabilised. The borders, bar and whiskers of the box plot represent the first (Q1)
# and third (Q3) quartiles, the median and the most extreme data points within 1.5x the interquartile range from Q1, Q3,
# respectively

#
blood_pbatch = import_Sim_fun("blood", blood_cs_df, "batches")
beagle_pbatch = import_Sim_fun("blood", blood_cs_df, "batches")
tS18_pbatch = import_Sim_fun("blood", blood_cs_df, "batches")

# blood_pbatch = blood_pbatch %>%
#   filter(sim_p0 >0) %>%
#   filter(sim_p0 <1)
blood_pbatch = left_join(blood_pbatch, blood_cs_sim)
previous_pbatch = blood_pbatch[,1]
previous_occbatch = blood_pbatch[,11]
vkid
log_fc_plist = c()
log_fc_occlist = c()

nrow(blood_pbatch)
observed_data = blood_pbatch$blood_int
for(batch in 2:10){
  current_pbatch = blood_pbatch[,batch]
  combined_pbatch = (previous_pbatch+current_pbatch)/2
  log_fc_p = log2(combined_pbatch/current_pbatch)
  
  #observed over simulated data
  current_occbatch = blood_pbatch[,batch+10]
  combined_occbatch = (previous_occbatch+current_occbatch)/2
  log_fc_occ = log2(combined_occbatch/current_occbatch)
  
  #print(log_fc)
  previous_pbatch = combined_pbatch
  previous_occbatch = combined_occbatch
  log_fc_plist = c(log_fc_plist, log_fc_p)
  log_fc_occlist = c(log_fc_occlist, log_fc_occ)
  
}
nrow(blood_pbatch)
log_fc_df = data.frame(category = as.factor(rep(1:9,each= 1007)),
                       log2_fc = log_fc_list) %>%
  filter(log2_fc >= quantile(log2_fc, 0.25) - 1.5 * IQR(log2_fc) & 
           log2_fc <= quantile(log2_fc, 0.75) + 1.5 * IQR(log2_fc))

min(log_fc_df$log2_fc)

ggplot(log_fc_df, aes(category,log2_fc))+
  geom_boxplot(aes(colour="lightblue"), outlier.shape = NA) +
  labs(title = "Blood p-values progression in Monte Carlo simulation",
       x="batch",
       y="log2 Fold Change")+
  theme_minimal()


###Plotting CS vs p-value#####
table(blood_cs_sim$category)
blood_pval_adj = blood_cs_sim[blood_cs_sim$category != "protein_coding",]
blood_pval_adj$category
blood_cs_sim$q_val_blood[1] = NA
tS18_cs_sim[tS18_cs_sim$category =="snRNA" &tS18_cs_sim$q_val_tS18 <0.05,]


ggplot(blood_cs_sim, aes(x=-log2(q_val_blood+0.001), y= log2(blood_cs +0.001), colour= category,
                         alpha = ifelse(q_val_blood <0.05 & blood_cs > 0.002, "significant", "insignificant"))) +
  geom_point(size = 1) +
  scale_alpha_manual(values = c("insignificant" =0.3))+
  geom_point(data = blood_cs_sim[2:3,], aes(x=-log2(q_val_blood+0.001), y= log2(blood_cs+0.001)), size =2, alpha =1) + 
  scale_colour_manual(values = colours_lst) + 
  geom_text_repel(aes(label = c(NA, "HMSbeagle", "3S18", rep(NA,nrow(blood_cs_sim)-3))), alpha =1) +  
  #scale_color_manual(values = c("cs>0.001 & p_val<0.01" = "red", 'insignificant'="grey"))+
  geom_vline(xintercept = -log2(0.05+0.001), colour = "red", linetype = 2) + 
  geom_hline(yintercept = log2(0.002+0.001), colour = "red", linetype = 2) +
  #annotate("rect", fill = "darkgreen", alpha = 0.08, xmin =-log2(0.01), ymin = log2(0.001), xmax = Inf, ymax = Inf) + 
  labs(title = "Blood RNA interaction p-values and connection scores",
       x="-log2(p-value)",
       y="log2(connection score)",
       color ="") +
  theme_minimal() +
  theme(legend.position = "none")


ggsave("002COMRADES_project/output/blood_cs_vs_p.tiff", height = 4, width = 3.5)

#Dividing up the squares

#blood_cs_sim[is.na(blood_cs_sim)] =0
blood_q1 = blood_cs_sim[blood_cs_sim$q_val_blood >0.05 & blood_cs_sim$blood_cs > 0.002,]
blood_q2_sig = blood_cs_sim[blood_cs_sim$q_val_blood <0.05 & blood_cs_sim$blood_cs > 0.002,][-1,]
blood_q3 = blood_cs_sim[blood_cs_sim$q_val_blood >0.05 & blood_cs_sim$blood_cs < 0.002,]
blood_q4 = blood_cs_sim[blood_cs_sim$q_val_blood <0.05 & blood_cs_sim$blood_cs < 0.002,]

total_abundance = sum(blood_cs_sim$total_interintra_int)
blood_quantiles = data.frame(quantile = c(rep("Q1", nrow(blood_q1)),rep("Q2", nrow(blood_q2_sig)),rep("Q3", nrow(blood_q3)),rep("Q4", nrow(blood_q4))),
                             relative_abundance = c(blood_q1$total_interintra_int/total_abundance, blood_q2_sig$total_interintra_int/total_abundance, blood_q3$total_interintra_int/total_abundance, blood_q4$total_interintra_int/total_abundance),
                             relative_int = c(blood_q1$blood_int/blood_q1$total_inter_int,blood_q2_sig$blood_int/blood_q2_sig$total_inter_int,blood_q3$blood_int/blood_q3$total_inter_int,blood_q4$blood_int/blood_q4$total_inter_int))

ggplot(blood_quantiles, aes(x=log2(relative_int), y = log2(relative_abundance), colour = quantile)) + 
  geom_point()+
  theme_minimal() +
  labs(title = "Blood interactors per quantile")


write.csv(blood_q2_sig, "blood_significant_int.csv")
#HMSBeagle

beagle_cs_sim$sim_pval[2] = NA
ggplot(beagle_cs_sim, aes(x=-log2(q_val_beagle+0.001), y= log2(beagle_cs+0.001), colour = category,
                          alpha = ifelse(q_val_beagle <0.05 & beagle_cs > 0.002, "significant", "insignificant"))) +
  geom_point( size = 1) +
  scale_alpha_manual(values = c("insignificant" =0.3))+
  geom_point(data = beagle_cs_sim[1:3,], aes(x=-log2(q_val_beagle+0.001), y= log2(beagle_cs +0.001), colour = category), size =2) + 
  #geom_text(label = beagle_cs_sim$RowNames, size= 2.5)+
  #geom_point(aes(alpha = ifelse(q_val_beagle <0.01 & beagle_cs>0.001, "cs>0.001 & p_val<0.01", "insignificant"))) +
  scale_colour_manual(values = colours_lst) + 
  geom_text_repel(aes(label = c("Blood", "HMSbeagle", "3S18", rep(NA,nrow(beagle_cs_sim)-3))), 
  alpha = 1) +  
  #scale_color_manual(values = c("cs>0.001 & p_val<0.01" = "red", 'insignificant'="grey"))+
  geom_vline(xintercept = -log2(0.05+0.001), colour = "red", linetype = 2) + 
  geom_hline(yintercept = log2(0.002+0.001), colour = "red", linetype = 2) +
  #annotate("rect", fill = "darkgreen", alpha = 0.08, xmin =-log2(0.01), ymin = log2(0.001), xmax = Inf, ymax = Inf) + 
  labs(title = "HMSbeagle RNA interaction p-values and connection scores",
       x="-log2(p-value)",
       y="log2(connection score)",
       color ="") +
  theme_minimal()  +
  theme(legend.position = "none")  
  

ggsave("002COMRADES_project/output/beagle_cs_vs_p.tiff", height = 4, width = 3.5)

?ggplot2
#Dividing up the squares
#beagle_cs_sim[is.na(beagle_cs_sim)] =0

beagle_q1 = beagle_cbeagle_cs_sim_3beagle_q1 = beagle_cs_sim[beagle_cs_sim$q_val_beagle >0.05 & beagle_cs_sim$beagle_cs > 0.002,]
beagle_q2_sig = beagle_cs_sim[beagle_cs_sim$q_val_beagle <0.05 & beagle_cs_sim$beagle_cs > 0.002,][-2,]
beagle_q3 = beagle_cs_sim[beagle_cs_sim$q_val_beagle >0.05 & beagle_cs_sim$beagle_cs < 0.002,]
beagle_q4 = beagle_cs_sim[beagle_cs_sim$q_val_beagle <0.05 & beagle_cs_sim$beagle_cs < 0.002,]

total_abundance = sum(beagle_cs_sim$total_interintra_int)
beagle_quantiles = data.frame(quantile = c(rep("Q1", nrow(beagle_q1)),rep("Q2", nrow(beagle_q2_sig)),rep("Q3", nrow(beagle_q3)),rep("Q4", nrow(beagle_q4))),
                             relative_abundance = c(beagle_q1$total_interintra_int/total_abundance, beagle_q2_sig$total_interintra_int/total_abundance, beagle_q3$total_interintra_int/total_abundance, beagle_q4$total_interintra_int/total_abundance),
                             relative_int = c(beagle_q1$beagle_int/beagle_q1$total_inter_int,beagle_q2_sig$beagle_int/beagle_q2_sig$total_inter_int,beagle_q3$beagle_int/beagle_q3$total_inter_int,beagle_q4$beagle_int/beagle_q4$total_inter_int))

ggplot(beagle_quantiles, aes(x=log2(relative_int), y = log2(relative_abundance), colour = quantile)) + 
  geom_point()+
  theme_minimal() +
  labs(title = "HMSBeagle interactors per quantile")

write.csv(beagle_q2_sig, "beagle_significant_int.csv")



#3S18

ggplot(tS18_cs_sim, aes(x=-log2(q_val_tS18+0.001), y= log2(tS18_cs+0.001), colour = category,
                        alpha = ifelse(q_val_tS18 <0.05 & tS18_cs > 0.002, "significant", "insignificant"))) +
  geom_point(size = 1) +
  scale_alpha_manual(values = c("insignificant" = 0.3)) +
  
  geom_point(data = tS18_cs_sim[1:3,], aes(x=-log2(q_val_tS18+0.001), y= log2(tS18_cs+0.001), colour = category), size =1) + 
  
  #geom_text(label = tS18_cs_sim$RowNames, size= 2.5)+
  #geom_point(aes(alpha = ifelse(q_val_tS18 <0.01 & tS18_cs>0.001, "cs>0.001 & p_val<0.01", "insignificant"))) +
  scale_colour_manual(values = colours_lst) + 
  geom_text_repel(aes(label = c("Blood", "HMSbeagle", "3S18", rep(NA,nrow(tS18_cs_sim)-3))), alpha=1) +  
  #scale_color_manual(values = c("cs>0.001 & p_val<0.01" = "red", 'insignificant'="grey"))+
  geom_vline(xintercept = -log2(0.05+0.001), colour = "red", linetype = 2) + 
  geom_hline(yintercept = log2(0.002+0.001), colour = "red", linetype = 2) +
  #annotate("rect", fill = "darkgreen", alpha = 0.08, xmin =-log2(0.01), ymin = log2(0.001), xmax = Inf, ymax = Inf) + 
  labs(title = "3S18 RNA interaction p-values and connection scores",
       x="-log2(p-value)",
       y="log2(connection score)",
       color ="") +
  
  theme_minimal() +
  theme(legend.position = "none")
  

ggsave("002COMRADES_project/output/tS18_cs_vs_p.tiff", height = 4, width = 3.5)


#Dividing up the squares
#tS18_cs_sim[is.na(tS18_cs_sim)] =0

tS18_q1 = tS18_cs_sim[tS18_cs_sim$q_val_tS18 >0.05 & tS18_cs_sim$tS18_cs > 0.002,]
tS18_q2_sig = tS18_cs_sim[tS18_cs_sim$q_val_tS18 <0.05 & tS18_cs_sim$tS18_cs > 0.002,][-1,]
tS18_q3 = tS18_cs_sim[tS18_cs_sim$q_val_tS18 >0.05 & tS18_cs_sim$tS18_cs < 0.002,]
tS18_q4 = tS18_cs_sim[tS18_cs_sim$q_val_tS18 <0.05 & tS18_cs_sim$tS18_cs < 0.002,]

total_abundance = sum(tS18_cs_sim$total_interintra_int)
tS18_quantiles = data.frame(quantile = c(rep("Q1", nrow(tS18_q1)),rep("Q2", nrow(tS18_q2_sig)),rep("Q3", nrow(tS18_q3)),rep("Q4", nrow(tS18_q4))),
                              relative_abundance = c(tS18_q1$total_interintra_int/total_abundance, tS18_q2_sig$total_interintra_int/total_abundance, tS18_q3$total_interintra_int/total_abundance, tS18_q4$total_interintra_int/total_abundance),
                              relative_int = c(tS18_q1$tS18_int/tS18_q1$total_inter_int,tS18_q2_sig$tS18_int/tS18_q2_sig$total_inter_int,tS18_q3$tS18_int/tS18_q3$total_inter_int,tS18_q4$tS18_int/tS18_q4$total_inter_int))

ggplot(tS18_quantiles, aes(x=log2(relative_int), y = log2(relative_abundance), colour = quantile)) + 
  geom_point()+
  theme_minimal() +
  labs(title = "3S18 interactors per quantile")

write.csv(tS18_q2_sig, "tS18_significant_int.csv")




