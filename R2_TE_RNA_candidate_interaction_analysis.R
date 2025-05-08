#Downstream analysis of significant interactions

library(devtools)
library(RColorBrewer)
library(seqinr)
#BiocManager::install("GenomicRanges")
library(GenomicRanges)
library(dplyr)
setwd("C:/Users/miabe/OneDrive/Desktop/Part III Systems/Project/")

devtools::load_all("C:/Users/miabe/OneDrive/Desktop/Part III Systems/Project/000Scripts/rnaCrosslinkOO")



#loading required files 

blood_q2_sig = read.csv("blood_significant_int.csv")
beagle_q2_sig = read.csv("beagle_significant_int.csv")
tS18_q2_sig = read.csv("tS18_significant_int.csv")

cdsBlood
cdsHMSbeagle
cdstS18



##########################################################################################################
####                                             FUNCTIONS                                            ####
##########################################################################################################



ceiling(cds3S18@rnaSize/20)

consistentInts2 = function(cds, sig_ids){
  
  #forward interaction df
  forward_int_trace_df = getInteractions(cds = cds, interactors = sig_ids)
  forward_int_trace_df = forward_int_trace_df[forward_int_trace_df$sample != "c3",]

    
  #Forward interactions divided in bins to check consistency among replicates
  prev_window = c(0,0)
  binned_forward_trace_df = data.frame()
  
  for(bin in 1:ceiling(cds@rnaSize/20)){
    current_window_TE = c(prev_window[2],prev_window[2]+20)
    if(current_window_TE[2]>cds@rnaSize){
      current_window_TE[2] = cds@rnaSize
    }
    #print(current_window_TE)

    current_trace_df = forward_int_trace_df %>%
      filter(Position >= current_window_TE[1] & Position <current_window_TE[2]) %>%
      group_by(rna, sample) %>%
      mutate(sum_bin = sum(depth)) %>%
      mutate(over_threshold = sum_bin >20) %>%
      ungroup()
    if(nrow(current_trace_df)==0){
      #print("skip")
      prev_window = current_window_TE
      next}
    current_trace_df = current_trace_df%>%
      group_by(rna, Position) %>%
      mutate(number_cons_replicates = sum(over_threshold)) %>%
      mutate(depth = ifelse(number_cons_replicates<3, 0, depth)) %>%
      ungroup()
    binned_forward_trace_df = rbind(binned_forward_trace_df, current_trace_df)
    prev_window = current_window_TE
  }
  
  #Reverse interactions creating dataframe and dividing in bins
  binned_rev_trace_df = data.frame()
  
  reverse_int_trace_df = data.frame()
  #print(consistent_peaks)
  for(transcript in unique(binned_forward_trace_df$rna)){
    #Reverse interaction creating the dataframe
    if(transcript=="FBtr0452214_FBgn0263847_LU_snRNA"){next}
    reverse_int_trace_current = getReverseInteractions(cds,transcript)
    #print(reverse_int_trace_current)
    if(nrow(reverse_int_trace_current)==0){
      next}
    reverse_int_trace_current = reverse_int_trace_current[reverse_int_trace_current$sample !="c3",]
    if(transcript =="FBtr0074208_FBgn0003920_14B_snRNA"){
    }
    reverse_int_trace_df = rbind(reverse_int_trace_df, reverse_int_trace_current)

    #reverse interaction dividing into bins per RNA
    prev_window = c(0,0)
    
    rna_size = max(reverse_int_trace_current$Position)
    for(bin in 1:ceiling(rna_size/20)){
      current_window_rna = c(prev_window[2],prev_window[2]+20)
      #print(current_window_rna)
      
      if(current_window_rna[2]>rna_size){
        current_window_rna[1] = current_window_rna[1] - (current_window_rna[2] - rna_size)
        current_window_rna[2] = rna_size
      }

      filt_reverse_int = reverse_int_trace_current %>%
        filter(Position >= current_window_rna[1] & Position <current_window_rna[2]) %>%
        group_by(rna, sample) %>%
        mutate(sum_bin = sum(depth)) %>%
        mutate(over_threshold = sum_bin >20) %>%
        ungroup() 
      if(nrow(filt_reverse_int)==0){
        #print("skip")
        prev_window = current_window_rna
        next}
      filt_reverse_int = filt_reverse_int%>%
        group_by(rna, Position) %>%
        mutate(number_cons_replicates = sum(over_threshold)) %>%
        mutate(depth = ifelse(number_cons_replicates<3, 0, depth)) %>%
        ungroup()
      binned_rev_trace_df = rbind(binned_rev_trace_df, filt_reverse_int)
      prev_window = current_window_rna
    }
  }

  #print("part2")
  #print(unique(binned_forward_trace_df$Position))
  #print(binned_forward_trace_df)
  
  
  #Adding mean depth to forward data frame
  if(nrow(binned_forward_trace_df >0)){

    consistent_peaks =  as.data.frame(binned_forward_trace_df) %>%
      group_by(rna, Position) %>%
      summarize(mean_depth = mean(depth, na.rm=TRUE), .groups ="drop") %>%
      ungroup() %>%
      group_by(rna) %>%
      filter(sum(mean_depth) >0)}
  else{
    consistent_peaks = binned_forward_trace_df
  }

  #adding mean depth to reverse data frame
  if(nrow(binned_rev_trace_df)>0){
    binned_rev_trace_df2 =  binned_rev_trace_df %>%
      group_by(rna, Position) %>%
      summarize(mean_depth = mean(depth, na.rm=TRUE), .groups ="drop") %>%
      ungroup()}
  
  #print("part3")
  
  #print(unique(consistent_peaks$mean_depth))
  
 return(list(forward_int_trace_df,consistent_peaks, reverse_int_trace_df,binned_rev_trace_df2))

}





RNA_types = c("rRNA", "miRNA","protein_coding","ncRNA","pseudogene","snoRNA","snRNA","tRNA")



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





##########################################################################################################
####                        PLOTTING THE INTERACTIONS FORWARD AND REVERSE                             ####
##########################################################################################################
consistentInts = function(cds, sig_ids){
  # forward interaction df
  forward_int_trace_df = getInteractions(cds = cds, interactors = sig_ids)
  forward_int_trace_df = forward_int_trace_df[forward_int_trace_df$sample %in% c("s1", "s2", "s3"),]
  
  # Forward interactions divided in bins to check consistency among replicates
  prev_window = c(0,0)
  binned_forward_trace_df = data.frame()
  
  for(bin in 1:ceiling(cds@rnaSize/20)){
    current_window_TE = c(prev_window[2], prev_window[2] + 20)
    if(current_window_TE[2] > cds@rnaSize){
      current_window_TE[2] = cds@rnaSize
    }
    
    current_trace_df = forward_int_trace_df %>%
      filter(Position >= current_window_TE[1] & Position < current_window_TE[2]) %>%
      group_by(rna, sample) %>%
      mutate(sum_bin = sum(depth)) %>%
      mutate(over_threshold = sum_bin > 20) %>%
      ungroup()
    
    if(nrow(current_trace_df) == 0){
      prev_window = current_window_TE
      next
    }
    
    current_trace_df = current_trace_df %>%
      group_by(rna, Position) %>%
      mutate(number_cons_replicates = sum(over_threshold)) %>%
      mutate(depth = ifelse(number_cons_replicates < 2, 0, depth)) %>%
      ungroup()
    
    binned_forward_trace_df = rbind(binned_forward_trace_df, current_trace_df)
    prev_window = current_window_TE
  }
  
  # Adding mean depth to forward data frame
  if(nrow(binned_forward_trace_df > 0)){
    consistent_peaks = as.data.frame(binned_forward_trace_df) %>%
      group_by(rna, Position) %>%
      summarize(mean_depth = mean(depth, na.rm = TRUE), .groups = "drop") %>%
      ungroup() %>%
      group_by(rna) %>%
      filter(sum(mean_depth) > 0)
  } else {
    consistent_peaks = binned_forward_trace_df
  }
  
  # Reverse interactions creating dataframe and dividing in bins
  reverse_int_trace_df = data.frame()
  rnas = unique(consistent_peaks$rna)
  for(transcript in rnas){
    # if(transcript == "FBtr0452214_FBgn0263847_LU_snRNA"){ next }
    # 
    reverse_int_trace_current = getReverseInteractions(cds, transcript)
    # if(nrow(reverse_int_trace_current) == 0){ next }
    # 
    # reverse_int_trace_current = reverse_int_trace_current[reverse_int_trace_current$sample %in% c("s1", "s2", "s3"),]
    # if(transcript == "FBtr0074208_FBgn0003920_14B_snRNA"){}
    # 
    reverse_int_trace_df = rbind(reverse_int_trace_df, reverse_int_trace_current)
    
    # reverse interaction dividing into bins per RNA
    prev_window = c(0, 0)
    rna_size = max(reverse_int_trace_current$Position)
  }
    
  #   for(bin in 1:ceiling(rna_size/20)){
  #     current_window_rna = c(prev_window[2], prev_window[2] + 20)
  #     if(current_window_rna[2] > rna_size){
  #       current_window_rna[1] = current_window_rna[1] - (current_window_rna[2] - rna_size)
  #       current_window_rna[2] = rna_size
  #     }
  #     
  #     filt_reverse_int = reverse_int_trace_current %>%
  #       filter(Position >= current_window_rna[1] & Position < current_window_rna[2]) %>%
  #       group_by(rna, sample) %>%
  #       mutate(sum_bin = sum(depth)) %>%
  #       mutate(over_threshold = sum_bin > 20) %>%
  #       ungroup()
  #     
  #     if(nrow(filt_reverse_int) == 0){
  #       prev_window = current_window_rna
  #       next
  #     }
  #     
  #     filt_reverse_int = filt_reverse_int %>%
  #       group_by(rna, Position) %>%
  #       mutate(number_cons_replicates = sum(over_threshold)) %>%
  #       mutate(depth = ifelse(number_cons_replicates < 3, 0, depth)) %>%
  #       ungroup()
  #     
  #     binned_rev_trace_df = rbind(binned_rev_trace_df, filt_reverse_int)
  #     prev_window = current_window_rna
  #   }
  # }
  
  # adding mean depth to reverse data frame
  if(nrow(reverse_int_trace_df) > 0){
    reverse_int_trace_df2 = reverse_int_trace_df %>%
      group_by(rna, Position) %>%
      summarize(mean_depth = mean(depth, na.rm = TRUE), .groups = "drop") %>%
      ungroup()
  }
  else{
    reverse_int_trace_df2 = reverse_int_trace_df
  }
  
  
  return(list(forward_int_trace_df, consistent_peaks, reverse_int_trace_df, reverse_int_trace_df2))
}




devtools::load_all("C:/Users/miabe/OneDrive/Desktop/Part III Systems/Project/000Scripts/rnaCrosslinkOO")


normalize_interaction_matrix <- function(M) {
  row_sums <- rowSums(M)
  col_sums <- colSums(M)
  
  # Avoid division by zero
  row_sums[row_sums == 0] <- 1
  col_sums[col_sums == 0] <- 1
  
  M_norm <- M / outer(row_sums, col_sums, "*")
  return(M_norm)
}
krNorm = function(TE, RNA_ids_lst){
  

isSymmetric.matrix(mat)
mat = plotInteractionsAverage(cdsBlood, rna = cdsBlood@rnas,RNA_ids[1], returnData= TRUE)
mat = as.matrix(mat)
View(mat)

aspectHeatmap(normalize_interaction_matrix(log2(mat+1)))
#############
####Blood####
#############

blood_sig_ids =c(blood_q2_sig$RowNames)


blood_domains <- data.frame(
  x_start = c(0, 966, 1863, 3749, 7010),
  x_end = c(400, 1271, 3166, 6733, 7410),
  domain = c("LTR1", "sORF", "ORF1", "Pol", "LTR2"),
  colours = c("#9B94B3","#ABC2C3","#ABC2C3","#ABC2C3","#9B94B3")
)


blood_q2_sig$category


#blood_sig_ids_snoRNA = blood_q2_sig$RowNames[blood_q2_sig$category=="snoRNA"]
blood_sig_ids_ncRNA = blood_q2_sig$RowNames[blood_q2_sig$category=="ncRNA"]
blood_sig_ids_tRNA = blood_q2_sig$RowNames[blood_q2_sig$category=="tRNA"]
blood_sig_ids_rRNA = blood_q2_sig$RowNames[blood_q2_sig$category=="rRNA"]

blood_sig_ids_miRNA = blood_q2_sig$RowNames[blood_q2_sig$category=="miRNA"]
blood_sig_ids_pseudogene = blood_q2_sig$RowNames[blood_q2_sig$category=="pseudogene"]
#blood_sig_ids_protein

#consistent_peaks = Blood_snRNAints[2]
#Blood_snoRNAints = consistentInts(cdsBlood, blood_sig_ids_snoRNA)
#blood_sig_ids_ncRNA = blood_q2_sig$RowNames[blood_q2_sig$category=="ncRNA"]

Blood_ncRNAints = consistentInts(cdsBlood, blood_sig_ids_ncRNA)
#Blood_tRNAints = consistentInts(cdsBlood, blood_sig_ids_tRNA)
Blood_rRNAints = consistentInts(cdsBlood, blood_sig_ids_rRNA)
#Blood_miRNAints = consistentInts(cdsBlood, blood_sig_ids_miRNA)
Blood_pseudogeneints = consistentInts(cdsBlood, blood_sig_ids_pseudogene)




###################
#### HMSbeagle ####
###################
beagle_domains <- data.frame(
  x_start = c(0, 1257, 3197, 6797),
  x_end = c(266, 2930, 6049, 7062),
  domain = c("LTR1", "GAG", "RT-pol", "LTR2"),
  colours = c("#9B94B3","#ABC2C3","#ABC2C3", "#9B94B3")
)


aubago


#beagle_sig_ids_snoRNA = beagle_q2_sig$RowNames[beagle_q2_sig$category=="snoRNA"]
beagle_sig_ids_ncRNA= beagle_q2_sig$RowNames[beagle_q2_sig$category=="ncRNA"][-1]
beagle_ncRNAints = consistentInts(cdsHMSbeagle, beagle_sig_ids_ncRNA)

beagle_sig_ids_rRNA = beagle_q2_sig$RowNames[beagle_q2_sig$category=="rRNA"][-1]

beagle_sig_ids_miRNA = beagle_q2_sig$RowNames[beagle_q2_sig$category=="miRNA"][-1]
beagle_sig_ids_pseudogene = beagle_q2_sig$RowNames[beagle_q2_sig$category=="pseudogene"][-1]


#beagle_snoRNAints = consistentInts(cdsHMSbeagle, beagle_sig_ids_snoRNA)
beagle_rRNAints = consistentInts(cdsHMSbeagle, beagle_sig_ids_rRNA)

beagle_miRNAints = consistentInts(cdsHMSbeagle, beagle_sig_ids_miRNA)
beagle_pseudogeneints = consistentInts(cdsHMSbeagle, beagle_sig_ids_pseudogene)

grep("FBgn0000003",tS18_cs_sim$RowNames)

tS18_cs_sim[4109,]




###################
#### 3S18 ####
###################
tS18_domains <- data.frame(
  x_start = c(0, 919, 5767),
  x_end = c(361, 5766, 6126),
  domain = c("LTR1","RT-pol + Integrase", "LTR2"),
  colours = c("#9B94B3","#ABC2C3", "#9B94B3")
)



tS18_sig_ids_snRNA = tS18_q2_sig$RowNames[tS18_q2_sig$category=="snRNA"]
#tS18_sig_ids_snoRNA = tS18_q2_sig$RowNames[tS18_q2_sig$category=="snoRNA"]
tS18_sig_ids_ncRNA = tS18_q2_sig$RowNames[tS18_q2_sig$category=="ncRNA"]
tS18_sig_ids_rRNA = tS18_q2_sig$RowNames[tS18_q2_sig$category=="rRNA"]

#tS18_sig_ids_miRNA = tS18_q2_sig$RowNames[tS18_q2_sig$category=="miRNA"]
tS18_sig_ids_pseudogene = tS18_q2_sig$RowNames[tS18_q2_sig$category=="pseudogene"]


tS18_snRNAints = consistentInts(cds3S18, tS18_sig_ids_snRNA)
#tS18_snoRNAints = consistentInts(cds3S18, tS18_sig_ids_snoRNA)

tS18_ncRNAints = consistentInts(cds3S18, tS18_sig_ids_ncRNA)

tS18_rRNAints = consistentInts(cds3S18, tS18_sig_ids_rRNA)
#tS18_miRNAints = consistentInts(cds3S18, tS18_sig_ids_miRNA)
tS18_pseudogeneints = consistentInts(cds3S18, tS18_sig_ids_pseudogene)




#############
####rRNA####
#############

#rRBA colours = purples
greens = brewer.pal(n=9, name="Greens")
purples = brewer.pal(n=9, name="Purples")
oranges = brewer.pal(n=9, name="Oranges")
blues = brewer.pal(n=9, name="Blues")

rRNA_colours = c(greens[5],greens[8],"turquoise",oranges[3],oranges[5], oranges[7],oranges[9], blues[4], blues[7], purples[5], purples[7], purples[9])
### Blood

#forward interactions
Blood_rRNA_consist_ints = Blood_rRNAints[[2]]

RNA_ids <- unique(Blood_rRNA_consist_ints$rna)
#RNA_ids = RNA_ids[-c(6,8,10,11,12,13)]
Blood_rRNA_consist_ints = Blood_rRNA_consist_ints[Blood_rRNA_consist_ints$rna %in% RNA_ids,]
blood_rRNA_original_ints = Blood_rRNAints[[1]]

transcript_types = c("primary_rRNA","18S_1/18S_2","2S_1","28S","primary_rRNA_2", "18S_1/18S_2","5.8S_1/5.8S_2/5.8S_3","2S_2","18S_3","5.8S_1/5.8S_2/5.8S_3","2S_3","5.8S_1/5.8S_2/5.8S_3","2S_4","primary_rRNA_3","5.8S_4")


#transcript_types = c("mit_LSU", "mt_SSU", "primary_rRNA","18S","2S","28S","primary_rRNA2", "18S2", "5.8S","2S2", "primary3", "18S3", "5.8S2", "2S3", "28S2", "5.8S3", "2S4", "primary4", "5.8S4")
#transcript_colours = c(setNames(oranges, RNA_ids))
transcript_names = setNames(transcript_types, RNA_ids)
Blood_rRNA_consist_ints$rna_name = sapply(Blood_rRNA_consist_ints$rna, function(x) transcript_names[x])


ggplot() +
  #geom_rect(data=blood_domains, aes(xmin = x_start, xmax = x_end, ymin = -Inf, ymax = Inf, fill = domain), alpha = 0.3) +
  #scale_fill_manual(values = structure(blood_domains$colours, names = blood_domains$domain)) +
  #geom_vline(xintercept = c(blood_domains$x_start, blood_domains$x_end), colour = "grey")+
  geom_line(data = Blood_rRNA_consist_ints[Blood_rRNA_consist_ints$rna_name =="2S_1",], aes(x = Position, y = mean_depth, color = rna_name),alpha = 0.7) +
  #facet_wrap(~rna_name) + 
  scale_colour_manual(values = rRNA_colours) + 
  labs(
       x = "Position on Blood",
       y = "Mean Read Depth") + 
  theme_minimal() + 
  theme(legend.position = "none")
ggsave(paste("002COMRADES_project/output/significant_trans_RNAint/blood_rRNA_filt.tiff", sep=""), height =4, width =4)

#reverse interactions
Blood_rRNA_revints = Blood_rRNAints[[4]]
#Blood_rRNA_revints = Blood_rRNA_revints[Blood_rRNA_revints$rna %in% RNA_ids,]
Blood_rRNA_revints$rna_name = sapply(Blood_rRNA_revints$rna, function(x) transcript_names[x])
rnatypes = setNames(c("primary_rRNA","18S","2S","28S","primary_rRNA", "18S","5.8S","2S","18S","5.8S","2S","5.8S","2S","primary_rRNA","5.8S")

                    ,transcript_types)
Blood_rRNA_revints$rnatype = sapply(Blood_rRNA_revints$rna_name, function(x) rnatypes[x])

ggplot() +
  #geom_rect(data=domains_df, aes(xmin = x_start, xmax = x_end, ymin = -Inf, ymax = Inf, fill = domain), alpha = 0.3) +
  #scale_fill_manual(values = domain_colours) +
  geom_line(data = Blood_rRNA_revints, aes(x = Position, y = mean_depth, color = rna_name),alpha = 0.7) +
  facet_wrap(~rnatype, scale="free") + 
  scale_colour_manual(values = rRNA_colours) + 
  labs(title = "Blood rRNA interactions",
       x = "Position on rRNA",
       y = "Mean Read Depth") + 
  theme_minimal()
ggsave(paste("002COMRADES_project/output/significant_trans_RNAint/blood_rRNA_filt_rev.tiff", sep=""), height =4, width =7.3)


#twoway plot - need to fix rnacrosslink package!!
for(rna in RNA_ids){
  print(rna)
  plotInteractionsAverage(cdsBlood, cdsBlood@rnas, rna, b="max",d="max", "C:/Users/miabe/OneDrive/Desktop/Part III Systems/Project/002COMRADES_project/output/significant_trans_RNAint/two_way_int_heatmaps/rRNA")
}




### HMSbeagle
beagle_rRNA_consist_ints = beagle_rRNAints[[2]]
ggplot(beagle_rRNAints[[1]], aes(x = Position, y=depth, colour = sample)) + 
  geom_line() + 
  facet_wrap(~rna, scale="free")
RNA_ids <- unique(beagle_rRNA_consist_ints$rna)

#Removing multimappers
#RNA_ids = RNA_ids[-c(6,8,10,11,12,13)]
beagle_rRNA_consist_ints = beagle_rRNA_consist_ints[beagle_rRNA_consist_ints$rna %in% RNA_ids,]
#transcript_types = c("primary_rRNA","18S","2S","28S","primary_rRNA_2", "5.8S", "18S_2","primary_rRNA_3","5.8S_2")

transcript_names = setNames(transcript_types, RNA_ids)

beagle_rRNA_consist_ints$rna_name = sapply(beagle_rRNA_consist_ints$rna, function(x) transcript_names[x])


ggplot() +
  #geom_vline(xintercept = c(beagle_domains$x_start, beagle_domains$x_end), colour = "grey")+
  
  geom_line(data = beagle_rRNA_consist_ints, aes(x = Position, y = mean_depth, color = rna_name),alpha = 0.7) +
  #facet_wrap(~rna_name, scale="free") + 
  scale_colour_manual(values = rRNA_colours) + 
  labs(
       x = "Position on HMSbeagle",
       y = "Mean Read Depth") + 
  theme_minimal() + 
  theme(legend.position = "none")
ggsave(paste("002COMRADES_project/output/significant_trans_RNAint/beagle_rRNA_filt.tiff", sep=""), height =4, width =4)




#reverse interactions
beagle_rRNA_revints = beagle_rRNAints[[4]]
beagle_rRNA_revints = beagle_rRNA_revints[beagle_rRNA_revints$rna %in% RNA_ids,]
beagle_rRNA_revints$rna_name = sapply(beagle_rRNA_revints$rna, function(x) transcript_names[x])
#rnatypes = setNames(c("primary_rRNA","18S","2S","28S","primary_rRNA", "5.8S", "18S","primary_rRNA","5.8S")
#,transcript_types)
beagle_rRNA_revints$rnatype = sapply(beagle_rRNA_revints$rna_name, function(x) rnatypes[x])

ggplot() +
  #geom_rect(data=domains_df, aes(xmin = x_start, xmax = x_end, ymin = -Inf, ymax = Inf, fill = domain), alpha = 0.3) +
  #scale_fill_manual(values = domain_colours) +
  geom_line(data = beagle_rRNA_revints, aes(x = Position, y = mean_depth, color = rna_name),alpha = 0.7) +
  facet_wrap(~rnatype, scale="free") + 
  scale_colour_manual(values = rRNA_colours) + 
  labs(title = "HMSbeagle rRNA interactions",
       x = "Position on rRNA",
       y = "Mean Read Depth") + 
  theme_minimal()
ggsave(paste("002COMRADES_project/output/significant_trans_RNAint/beagle_rRNA_filt_rev.tiff", sep=""), height =4, width =6.5)

#twoway plot - need to fix rnacrosslink package!!
for(rna in RNA_ids){
  print(rna)
  plotInteractionsAverage(cdsHMSbeagle, cdsHMSbeagle@rnas, rna, b="max",d="max", "C:/Users/miabe/OneDrive/Desktop/Part III Systems/Project/002COMRADES_project/output/significant_trans_RNAint/two_way_int_heatmaps/rRNA")
}




#3S18
tS18_rRNA_consist_ints = tS18_rRNAints[[2]]

RNA_ids <- unique(tS18_rRNA_consist_ints$rna)
transcript_names
#RNA_ids = RNA_ids[-c(6,8,10,11,12,13)]

#transcript_types = c("primary_rRNA","18S","2S","28S","primary_rRNA_2", "5.8S", "18S_2","primary_rRNA_3","5.8S_2")

transcript_names = setNames(transcript_types, RNA_ids)
tS18_rRNA_consist_ints$rna_name = sapply(tS18_rRNA_consist_ints$rna, function(x) transcript_names[x])
tS18_rRNA_consist_ints = tS18_rRNA_consist_ints[tS18_rRNA_consist_ints$rna %in% RNA_ids,]

ggplot() +
 # geom_vline(xintercept = c(tS18_domains$x_start, tS18_domains$x_end), colour = "grey")+
  
 geom_line(data = tS18_rRNA_consist_ints, aes(x = Position, y = mean_depth, color = rna_name),alpha = 0.7) +
  #facet_wrap(~rna) + 
  scale_colour_manual(values = rRNA_colours) + 
  labs(
       x = "Position on 3S18",
       y = "Mean Read Depth") + 
  theme_minimal() + 
  theme(legend.position = "none")


ggsave(paste("002COMRADES_project/output/significant_trans_RNAint/tS18_rRNA_filt.tiff", sep=""), height =4, width =4)
substr(blood_seq, 400,442)

#reverse interactions
tS18_rRNA_revints = tS18_rRNAints[[4]]
tS18_rRNA_revints = tS18_rRNA_revints[tS18_rRNA_revints$rna %in% RNA_ids,]
tS18_rRNA_revints$rna_name = sapply(tS18_rRNA_revints$rna, function(x) transcript_names[x])
tS18_rRNA_revints$rnatype = sapply(tS18_rRNA_revints$rna_name, function(x) rnatypes[x])
ggplot() +
  #geom_rect(data=domains_df, aes(xmin = x_start, xmax = x_end, ymin = -Inf, ymax = Inf, fill = domain), alpha = 0.3) +
  #scale_fill_manual(values = domain_colours) +
  geom_line(data = tS18_rRNA_revints, aes(x = Position, y = mean_depth, color = rna_name),alpha = 0.7) +
  facet_wrap(~rnatype, scale="free") + 
  scale_colour_manual(values = rRNA_colours) + 
  labs(title = "3S18 rRNA interactions",
       x = "Position on 3S18",
       y = "Mean Read Depth") + 
  theme_minimal()
ggsave(paste("002COMRADES_project/output/significant_trans_RNAint/tS18_rRNA_filt_rev.tiff", sep=""), height =4, width =6.5)




for(rna in RNA_ids){
  print(rna)
  plotInteractionsAverage(cds3S18, cds3S18@rnas, rna, b="max",d="max", "C:/Users/miabe/OneDrive/Desktop/Part III Systems/Project/002COMRADES_project/output/significant_trans_RNAint/two_way_int_heatmaps/rRNA")
}



##################
####pseudogene####
##################

pseudogene_colours = c(greens[4], greens[8],"#C0F6F7","#81EDEF", "#24DFE3","#15A0A2", "#0C5C5D",blues[8])


### Blood

#forward interactions
Blood_pseudogene_consist_ints = Blood_pseudogeneints[[2]]

RNA_ids <- unique(Blood_pseudogene_consist_ints$rna)
Blood_pseudogene_original_ints = Blood_pseudogeneints[[1]]

ggplot(Blood_pseudogene_original_ints, aes(x = Position, y=depth, colour = sample)) + 
  geom_line() + 
  facet_wrap(~rna, scale="free")


transcript_types = c("28S_pseudo", "28S_pseudo2", "5.8S_pseudo","28S_pseudo3","18S_pseudo","28S_pseudo4","18S_pseudo2", "28S_pseudo5")

transcript_names = setNames(transcript_types, RNA_ids)
Blood_pseudogene_consist_ints$rna_name = sapply(Blood_pseudogene_consist_ints$rna, function(x) transcript_names[x])

ggplot() +
  #geom_vline(xintercept = c(blood_domains$x_start, blood_domains$x_end), colour = "grey")+
  geom_line(data = Blood_pseudogene_consist_ints, aes(x = Position, y = mean_depth, color = rna_name),alpha = 0.7) +
  #facet_wrap(~rna) + 
  scale_colour_manual(values = pseudogene_colours) + 
  labs(
       x = "Position on Blood",
       y = "Mean Read Depth") + 
  theme_minimal() + 
  theme(legend.position = "none")
ggsave(paste("002COMRADES_project/output/significant_trans_RNAint/blood_pseudogene_filt.tiff", sep=""), height =4, width =4)

#reverse interactions
Blood_pseudogene_revints = Blood_pseudogeneints[[4]]
Blood_pseudogene_revints$rna_name = sapply(Blood_pseudogene_revints$rna, function(x) transcript_names[x])
rnatypes = setNames(c("28S_pseudo", "28S_pseudo", "5.8S_pseudo","28S_pseudo","18S_pseudo","28S_pseudo","18S_pseudo", "28S_pseudo"),
                    
                    transcript_types)
Blood_pseudogene_revints$rnatype = sapply(Blood_pseudogene_revints$rna_name, function(x) rnatypes[x])

ggplot() +
  #geom_rect(data=domains_df, aes(xmin = x_start, xmax = x_end, ymin = -Inf, ymax = Inf, fill = domain), alpha = 0.3) +
  #scale_fill_manual(values = domain_colours) +
  geom_line(data = Blood_pseudogene_revints, aes(x = Position, y = mean_depth, color = rna_name),alpha = 0.7) +
  facet_wrap(~rnatype, scale="free") + 
  scale_colour_manual(values = pseudogene_colours) + 
  labs(title = "Blood pseudogene interactions",
       x = "Position on Blood",
       y = "Mean Read Depth") + 
  theme_minimal()
ggsave(paste("002COMRADES_project/output/significant_trans_RNAint/blood_pseudogene_filt_rev.tiff", sep=""), height =3, width =6.2)


#twoway plot - need to fix rnacrosslink package!!
for(rna in RNA_ids){
  plotInteractionsAverage(cdsBlood, "AY180916_AY180916_AY180916_AY180916", "FBtr0074208_FBgn0003920_14B_pseudogene", b="max",d="max", "002COMRADES_project/output/significant_trans_RNAint/")
}



### HMSbeagle
pseudogene_colours = c(greens[4], greens[8],"#C0F6F7","#81EDEF", "#24DFE3","#15A0A2", "#0C5C5D",blues[5],blues[8])


beagle_pseudogene_consist_ints = beagle_pseudogeneints[[2]]
RNA_ids <- unique(beagle_pseudogene_consist_ints$rna)
transcript_names


transcript_types = c("28S_pseudo", "28S_pseudo2", "5.8S_pseudo","28S_pseudo3","18S_pseudo","28S_pseudo4","5.8S_pseudo2","18S_pseudo2", "28S_pseudo5")
transcript_names = setNames(transcript_types, RNA_ids)
beagle_pseudogene_consist_ints$rna_name = sapply(beagle_pseudogene_consist_ints$rna, function(x) transcript_names[x])


ggplot() +
  #geom_vline(xintercept = c(beagle_domains$x_start, beagle_domains$x_end), colour = "grey")+
  geom_line(data = beagle_pseudogene_consist_ints, aes(x = Position, y = mean_depth, color = rna_name),alpha = 0.7) +
  #facet_wrap(~rna) + 
  scale_colour_manual(values = pseudogene_colours) + 
  labs(
       x = "Position on HMSbeagle",
       y = "Mean Read Depth") + 
  theme_minimal() + 
  theme(legend.position = "none")
ggsave(paste("002COMRADES_project/output/significant_trans_RNAint/beagle_pseudogene_filt.tiff", sep=""), height =4, width =4)




#reverse interactions
beagle_pseudogene_revints = beagle_pseudogeneints[[4]]
beagle_pseudogene_revints$rna_name = sapply(beagle_pseudogene_revints$rna, function(x) transcript_names[x])
rnatypes = setNames(c("28S_pseudo", "28S_pseudo", "5.8S_pseudo","28S_pseudo","18S_pseudo","28S_pseudo","5.8S_pseudo","18S_pseudo", "28S_pseudo"),
                    
                    transcript_types)
beagle_pseudogene_revints$rnatype = sapply(beagle_pseudogene_revints$rna_name, function(x) rnatypes[x])

ggplot() +
  #geom_rect(data=domains_df, aes(xmin = x_start, xmax = x_end, ymin = -Inf, ymax = Inf, fill = domain), alpha = 0.3) +
  #scale_fill_manual(values = domain_colours) +
  
  geom_line(data = beagle_pseudogene_revints, aes(x = Position, y = mean_depth, color = rna_name),alpha = 0.7) +
  facet_wrap(~rnatype, scale="free") + 
  scale_colour_manual(values = pseudogene_colours) + 
  labs(title = "HMSbeagle pseudogene interactions",
       x = "Position on HMSbeagle",
       y = "Mean Read Depth") + 
  theme_minimal() + 
  theme(legend.position = "none")
ggsave(paste("002COMRADES_project/output/significant_trans_RNAint/beagle_pseudogene_filt_rev.tiff", sep=""), height =3, width =9)


#3S18
tS18_pseudogene_consist_ints = tS18_pseudogeneints[[2]]

RNA_ids <- unique(tS18_pseudogene_consist_ints$rna)

transcript_names = setNames(transcript_types, RNA_ids)
tS18_pseudogene_consist_ints$rna_name = sapply(tS18_pseudogene_consist_ints$rna, function(x) transcript_names[x])


ggplot() +
  #geom_vline(xintercept = c(tS18_domains$x_start, tS18_domains$x_end), colour = "grey")+
  geom_line(data = tS18_pseudogene_consist_ints, aes(x = Position, y = mean_depth, color = rna_name),alpha = 0.7) +
  #facet_wrap(~rna) + 
  scale_colour_manual(values = pseudogene_colours) + 
  labs(
       x = "Position on 3S18",
       y = "Mean Read Depth") + 
  theme_minimal() + 
  theme(legend.position = "none")


ggsave(paste("002COMRADES_project/output/significant_trans_RNAint/tS18_pseudogene_filt.tiff", sep=""), height =4, width =4)


#reverse interactions
tS18_pseudogene_revints = tS18_pseudogeneints[[4]]
tS18_pseudogene_revints$rna_name = sapply(tS18_pseudogene_revints$rna, function(x) transcript_names[x])
tS18_pseudogene_revints$rnatype = sapply(tS18_pseudogene_revints$rna_name, function(x) rnatypes[x])
ggplot() +
  #geom_rect(data=domains_df, aes(xmin = x_start, xmax = x_end, ymin = -Inf, ymax = Inf, fill = domain), alpha = 0.3) +
  #scale_fill_manual(values = domain_colours) +
  geom_line(data = tS18_pseudogene_revints, aes(x = Position, y = mean_depth, color = rna_name),alpha = 0.7) +
  facet_wrap(~rnatype, scale="free") + 
  scale_colour_manual(values = pseudogene_colours) + 
  labs(title = "3S18 pseudogene interactions",
       x = "Position on 3S18",
       y = "Mean Read Depth") + 
  theme_minimal()
ggsave(paste("002COMRADES_project/output/significant_trans_RNAint/tS18_pseudogene_filt_rev.tiff", sep=""), height =3, width =6.2)













#############
####tRNA####
#############
colours_lst = c("Blood" = "#b04d48", "HMSbeagle" = "#e0aa24", "3S18" = "#2a4864","protein_coding" = "lightgrey","rRNA"= "#8dbd05","ncRNA"= "#00a1ae","snRNA"= "#5e36cc", "snoRNA"="#fe318e", "miRNA"="#d9b6db", "pseudogene"="#fb5607", "tRNA" = "#9D0645")

colors = brewer.pal(8, "Dark2")
#tRNA_colours = c("#4F0322","#9D0642", "#EC0964","#F9629E")
### Blood

Blood_tRNAints = consistentInts(cdsBlood, blood_sig_ids_tRNA)


#forward interactions
Blood_tRNA_consist_ints = Blood_tRNAints[[2]]
RNA_ids = unique(Blood_tRNA_consist_ints$rna)

sum(Blood_tRNAints[[1]]$depth)/cdsBlood@rnaSize
sum(Blood_tRNAints[[3]]$depth)/73

nrow(blood_q2_sig[blood_q2_sig$category =="tRNA",])
# Corresponding transcript types
#transcript_types <- c("Arg_TCG_2-1","Arg_TCG_3-1","Arg_TCG_3-2","Arg_TCG_3-3","Arg_TCG_3-4","Arg_TCG_2-2","Arg_TCG_4-1","Arg_TCG_2-3","Arg_TCG_2-4")
#length(transcript_types)
transcript_types <- c("Arg_TCG_2-1/-2/-3", "Arg_TCG_3-1/-2/-3/-4","Arg_TCG_3-1/-2/-3/-4","Arg_TCG_3-1/-2/-3/-4","Arg_TCG_3-1/-2/-3/-4","Arg_TCG_2-1/-2/-3","Arg_TCG_4-1","Arg_TCG_2-1/-2/-3","Arg_TCG_2-4")


transcript_names = setNames(transcript_types, RNA_ids)
Blood_tRNA_consist_ints$rna_name = sapply(Blood_tRNA_consist_ints$rna, function(x) transcript_names[x])


ggplot() +
  #geom_rect(data=blood_domains, aes(xmin = x_start, xmax = x_end, ymin = -Inf, ymax = Inf, fill = domain), alpha = 0.3) +
  #scale_fill_manual(values = structure(blood_domains$colours, names = blood_domains$domain)) +
  geom_rect(data=blood_domains, aes(xmin = 350, xmax = 550, ymin = -0.3, ymax = -0.1), fill="lightgrey",alpha =0.1) +
  
  geom_rect(data=blood_domains, aes(xmin = 400, xmax = 418, ymin = -0.3, ymax = -0.1), fill="#9D0642",alpha =1) +
  

  geom_line(data = subset(Blood_tRNA_consist_ints, Position <550 &Position >350), aes(x = Position, y = mean_depth, color = rna_name),alpha = 0.7) +
  
  #facet_wrap(~rna) + 
  scale_colour_manual(values = colors) + 
  labs(
       x = "Position on Blood",
       y = "Mean Read Depth") + 
  theme_minimal() + theme(legend.position = "bottom")
ggsave(paste("002COMRADES_project/output/significant_trans_RNAint/blood_tRNA_filt.tiff", sep=""), height =4, width =5.3)


#reverse interactions
Blood_tRNA_revints = Blood_tRNAints[[4]]
Blood_tRNA_revints$rna_name = sapply(Blood_tRNA_revints$rna, function(x) transcript_names[x])
#rnatypes = c(setNames(c(rep("Arg_TCG", 4)), transcript_types))
#Blood_tRNA_revints$rnatype = sapply(Blood_tRNA_revints$rna_name, function(x) rnatypes[x])



ggplot() +
  #geom_rect(data=domains_df, aes(xmin = x_start, xmax = x_end, ymin = -Inf, ymax = Inf, fill = domain), alpha = 0.3) +
  #scale_fill_manual(values = domain_colours) +
  geom_rect(data=blood_domains, aes(xmin = 0, xmax = 63, ymin = -0.5, ymax = -0.1), fill="lightgrey",alpha =0.1) +
  
  geom_rect(data=blood_domains, aes(xmin = 58, xmax = 73, ymin = -0.5, ymax = -0.1), fill="#9D0642",alpha =1) +
  
  geom_line(data = Blood_tRNA_revints, aes(x = Position, y = mean_depth, color = rna_name),alpha = 0.7) +
  #facet_wrap(~rnatype, scale="free") + 
  scale_colour_manual(values = colors) + 
  #facet_wrap(~rna) +
  labs(
       x = "Position on tRNA",
       y = "Mean Read Depth") + 
  theme_minimal() + theme(legend.position = "none")
ggsave(paste("002COMRADES_project/output/significant_trans_RNAint/blood_tRNA_filt_rev.tiff", sep=""), height =4, width =5.3)
max(Blood_tRNA_revints$Position)



plotInteractionsAverage(cdsBlood, cdsBlood@rnas, "FBtr0073857_FBgn0011952_Arg-TCG-2-1_tRNA",a =350, b=550,d="max", "C:/Users/miabe/OneDrive/Desktop/Part III Systems/Project/002COMRADES_project/output/significant_trans_RNAint/two_way_int_heatmaps/tRNA")


#RNA folding

#Blood sequence
#Blood sequence associated with subset 
TE_fasta = read.fasta(file = "002COMRADES_project/input/TE_sequences.fa")
blood_seq = getSequence.list(TE_fasta, as.string = TRUE)[[3]][[1]]
substr(blood_seq, 375,575)
#tRNA motif: Arg-TCG	cgaccgtgaca - not found the whole thing, but cgaccgtga
grep("cgaccgtga",blood_tRNA_int_seq)

#tRNA sequence

annotation_fasta = read.fasta(file = "002COMRADES_project/input/Dmel_hyb_full_TE_cdna_ncrna.fasta")

tRNA_seq = getSequence.list(annotation_fasta["FBtr0073857_FBgn0011952_Arg-TCG-2-1_tRNA"], as.string = TRUE)[[1]][[1]]
substr(tRNA_seq, 35,73)
paste(substr(blood_seq, 350,550),"&",substr(tRNA_seq, 1,73),"cca", sep="" )





###
### HMSbeagle
###


beagle_sig_ids_tRNA = beagle_q2_sig$RowNames[beagle_q2_sig$category=="tRNA"][-1]
beagle_tRNAints = consistentInts(cdsHMSbeagle, beagle_sig_ids_tRNA)


beagle_tRNA_consist_ints = beagle_tRNAints[[2]]

beagle_tRNA_original_ints = beagle_tRNAints[[1]]
RNA_ids = unique(beagle_tRNA_original_ints$rna)
ggplot(beagle_tRNA_original_ints, aes(x = Position, y=depth, colour = sample)) + 
  geom_line() +
  facet_wrap(~rna)

transcript_types <- c("Lys_TTT2-1/2/3/4","Lys_TTT2-1/2/3/4","Lys_TTT2-1/2/3/4","Lys_TTT2-1/2/3/4","Lys_TTT2-5")
beagle_tRNA_consist_ints =beagle_tRNA_consist_ints[beagle_tRNA_consist_ints$rna %in% RNA_ids,]

transcript_names = setNames(transcript_types, RNA_ids)
beagle_tRNA_consist_ints$rna_name = sapply(beagle_tRNA_consist_ints$rna, function(x) transcript_names[x])

#Motif in HMSbeagle for PBS find: cgcccgaacagggac
beagle_seq = getSequence.list(TE_fasta, as.string = TRUE)[[1]][[1]]
#tRNA motif: Arg-TCG	cgaccgtgaca - not found the whole thing, but cgaccgtga
grep("cgcccaacgtggggc",beagle_seq)

beagle_tRNA_added = beagle_tRNA_consist_ints[1,]
beagle_tRNA_added$Position = 272
beagle_tRNA_added$mean_depth = 0

beagle_tRNA_consist_ints = rbind(beagle_tRNA_consist_ints, beagle_tRNA_added)
ggplot() +
  geom_rect(aes(xmin = 250, xmax = 450, ymin = -0.18, ymax = -0.1), fill="lightgrey",alpha =0.3) +
  
  geom_rect(aes(xmin = 272, xmax = 290, ymin = -0.175, ymax = -0.1), fill="#9D0642",alpha =1) +
  geom_line(data = subset(beagle_tRNA_consist_ints, Position <450 &Position >250), aes(x = Position, y = mean_depth, color = rna_name),alpha = 0.7) +
  #facet_wrap(~rna) + 
  scale_colour_manual(values = colors[5:6]) + 
  labs(x = "Position on HMSbeagle",
       y = "Mean Read Depth") + 
  theme_minimal() + theme(legend.position = "bottom")
ggsave(paste("002COMRADES_project/output/significant_trans_RNAint/beagle_tRNA_filt.tiff", sep=""), height =4, width =5.3)
#reverse interactions
beagle_tRNA_consist_revints = beagle_tRNAints[[4]]
beagle_tRNA_consist_revints = beagle_tRNA_consist_revints[beagle_tRNA_consist_revints$rna %in% RNA_ids,]
beagle_tRNA_consist_revints = beagle_tRNA_consist_revints[beagle_tRNA_consist_revints$rna %in% RNA_ids,]
beagle_tRNA_consist_revints$rna_name = sapply(beagle_tRNA_consist_revints$rna, function(x) transcript_names[x])


#Find gtccctgttcgggcg
tRNA_lysttt = getSequence.list(annotation_fasta[2/], as.string = TRUE)[[1]][[1]]
grep("gtccctgttcgggcg", substr(tRNA_lysttt, 59,73))

ggplot() +
  geom_rect(aes(xmin = 0, xmax = 59, ymin = -0.18, ymax = -0.1), fill="lightgrey",alpha =0.3) +
  geom_rect(aes(xmin = 58, xmax = 73, ymin = -0.18, ymax = -0.1), fill="#9D0642",alpha =1) +
  geom_line(data = beagle_tRNA_consist_revints, aes(x = Position, y = mean_depth, color = rna_name),alpha = 0.7) +
  #facet_wrap(~rnatype, scale="free") + 
  scale_colour_manual(values = colors[5:6]) + 
  #facet_wrap(~rna) +
  labs(
       x = "Position on tRNA",
       y = "Mean Read Depth") + 
  theme_minimal() + theme(legend.position = "none")
ggsave(paste("002COMRADES_project/output/significant_trans_RNAint/beagle_tRNA_filt_rev.tiff", sep=""), height =4, width =5.3)


#twoway plot - need to fix rnacrosslink package!!
for(rna in RNA_ids){
  print(rna)
  plotInteractionsAverage(cdsBlood, cdsBlood@rnas, rna, b="max",d="max", "C:/Users/miabe/OneDrive/Desktop/Part III Systems/Project/002COMRADES_project/output/significant_trans_RNAint/two_way_int_heatmaps/tRNA")
}

#Find gtccctgttcgggcg
paste(substr(beagle_seq, 250,425),"&", tRNA_lysttt, sep="")



#3S18

tS18_sig_ids_tRNA = tS18_q2_sig$RowNames[tS18_q2_sig$category=="tRNA"]

tS18_tRNAints = consistentInts(cds3S18, tS18_sig_ids_tRNA)


tS18_tRNA_consist_ints = tS18_tRNAints[[2]]


tS18_tRNA_original_ints = tS18_tRNAints[[1]]

ggplot(tS18_tRNA_original_ints, aes(x = Position, y=depth, colour = sample)) + 
  geom_line() +
  facet_wrap(~rna)



#finding sequence ctctgctattg

tS18_seq = getSequence.list(TE_fasta, as.string = TRUE)[[2]][[1]]


grep("tccttcgagccggat", substr(tS18_seq, 372,386))


RNA_ids <- unique(tS18_tRNA_consist_ints$rna)
transcript_types <- c("Tyr_GTA_1-2/3/3/5/3/8/3-1","Tyr_GTA_1-2/3/3/5/3/8/3-1","Tyr_GTA_1-6","Tyr_GTA_1-2/3/3/5/3/8/3-1","Tyr_GTA_1-2/3/3/5/3/8/3-1","Tyr_GTA_1-2/3/3/5/3/8/3-1","Tyr_GTA_1-2/3/3/5/3/8/3-1","Tyr_GTA_1-2/3/3/5/3/8/3-1","Tyr_GTA_1-2/3/3/5/3/8/3-1")
transcript_names = setNames(transcript_types, RNA_ids)
tS18_tRNA_consist_ints$rna_name = sapply(tS18_tRNA_consist_ints$rna, function(x) transcript_names[x])

ggplot() +
  geom_rect(aes(xmin = 250, xmax = 530, ymin = -0.225, ymax = -0.1), fill="lightgrey",alpha =0.3) +
  geom_rect(aes(xmin = 372, xmax = 386, ymin = -0.225, ymax = -0.1), fill="#9D0642",alpha =1) +
  geom_line(data = subset(tS18_tRNA_consist_ints,Position >250 &Position<530), aes(x = Position, y = mean_depth, color = rna_name),alpha = 1) +
  #facet_wrap(~rna) + 
  scale_colour_manual(values = c(colors[7], "darkblue")) + 
  labs(x = "Position on 3S18",
       y = "Mean Read Depth") + 
  theme_minimal() + theme(legend.position = "bottom")


ggsave(paste("002COMRADES_project/output/significant_trans_RNAint/tS18_tRNA_filt.tiff", sep=""), height =4, width =5.3)



#reverse interactions
tS18_tRNA_revints = tS18_tRNAints[[4]]
tS18_tRNA_revints = tS18_tRNA_revints[tS18_tRNA_revints$rna %in% RNA_ids,]
tS18_tRNA_revints$rna_name = sapply(tS18_tRNA_revints$rna, function(x) transcript_names[x])

#Find atccggctcgaagga
tRNA_tyrGTA = getSequence.list(annotation_fasta[RNA_ids[1]], as.string = TRUE)[[1]][[1]]
grep("atccggctcgaagga", substr(tRNA_tyrGTA, 59,73))



ggplot() +
  geom_rect(aes(xmin = 0, xmax = 59, ymin = -0.4, ymax = -0.1), fill="lightgrey",alpha =0.3) +
  geom_rect(aes(xmin = 58, xmax = 73, ymin = -0.4, ymax = -0.1), fill="#9D0642",alpha =1) +
  geom_line(data = tS18_tRNA_revints, aes(x = Position, y = mean_depth, color = rna_name),alpha = 0.7) +
  #facet_wrap(~rna, scale="free") + 
  scale_colour_manual(values =  c(colors[7], "darkblue")) + 
  labs(
       x = "Position on tRNA",
       y = "Mean Read Depth") + 
  theme_minimal() + theme(legend.position = "none")
ggsave(paste("002COMRADES_project/output/significant_trans_RNAint/tS18_tRNA_filt_rev.tiff", sep=""), height =4, width =5.3)




#twoway plot - need to fix rnacrosslink package!!
for(rna in RNA_ids){
  print(rna)
  plotInteractionsAverage(cds3S18, cds3S18@rnas, rna, b="max",d="max", "C:/Users/miabe/OneDrive/Desktop/Part III Systems/Project/002COMRADES_project/output/significant_trans_RNAint/two_way_int_heatmaps/tRNA")
}



#For RNA cofold
paste(substr(tS18_seq, 280,560), "&", tRNA_tyrGTA)





#############
####snRNA####
#############

#There are some reads which are definitely multimapped:
#Combined: U6_96Aa 96Ab and 96Ac

grep("atgatga", substr(blood_seq, 5724,5730))

### Blood
blood_sig_ids_snRNA = blood_q2_sig$RowNames[blood_q2_sig$category=="snRNA"]
Blood_snRNAints = consistentInts(cdsBlood, blood_sig_ids_snRNA)
#splice site prediction

SUB


#forward interactions
Blood_snRNA_consist_ints = Blood_snRNAints[[2]]
#Blood_snRNA_ncRNA_consist_ints = rbind(Blood_snRNA_consist_ints, Blood_ncRNA_consist_ints)
#Blood_snRNA_original_ints= Blood_snRNAints[[1]]
ggplot(data=Blood_snRNA_original_ints, aes(x=Position, y=depth, colour = sample)) +
  geom_line() + facet_wrap(~rna)

RNA_ids <- unique(Blood_snRNA_consist_ints$rna)


# Corresponding transcript types
#RNA_ids = RNA_ids[c(1:8,11,15)]
transcript_types <- c("U2_14B/34ABa/38ABa snRNA", "U1_21D/95Ca/95Cb snRNA", "U2_34ABb/34ABc snRNA","U2_14B/34ABa/38ABa snRNA","U2_38ABb snRNA","U2_14B/34ABa/38ABa snRNA","U1_95Cc snRNA","U1_21D/95Ca/95Cb snRNA","U1_21D/95Ca/95Cb snRNA","U6_96Aa/96Ab/96Ac snRNA","U6_96Aa/96Ab/96Ac snRNA","U6_96Aa/96Ab/96Ac snRNA","U2_34ABb/34ABc snRNA","U7 snRNA")
length(transcript_types)
greens = brewer.pal(n=9, name="Greens")
purples = brewer.pal(n=9, name="Purples")
oranges = brewer.pal(n=9, name="Oranges")
blues = brewer.pal(n=9, name="Blues")


# 
# colors <- c(
#   "U5_63Bc" = blues[9],
#   "U2_14B" = greens[4],
#   "U5_14B" = blues[3],
#   "U5_23D" = blues[4],
#   "U1_21D" = purples[6], 
#   "U1_82Eb" = purples[7],# U1 group - purple
#   "U2_34ABb" = greens[6],
#   "U5_34A" = blues[5], # U2 group - green
#   "U2_34ABa" = greens[5],
#   "U5_35D"= blues[6],
#   "U2_38ABb" = greens[8],
#   "U5_38ABb" = blues[7],
#   "U2_38ABa" = greens[7],
#   "U5_38ABa" = blues[8],
#   #"U1_95Cb" = purples[3],
#   "U1_96Cc" = purples[9],
#   "U1_95Cc" = purples[8], 
#   "U1_95Ca" = purples[5],  
#   "U6_96A" = oranges[8],  
#   #"U6_96Ab" = oranges[2], 
#   #"U6_96Ac" = oranges[3],  
#   "U2_34ABc" = greens[6],  
#   "U7" = "#998833",
#   "Like-U" = "magenta"
# )

colors = c(blues[4], blues[6], greens[3], greens[5], greens[8], purples[8], oranges[4], "brown")

#transcript_colours = c(setNames(colors, RNA_ids))
transcript_names = setNames(transcript_types, RNA_ids)
Blood_snRNA_consist_ints$rna_name = sapply(Blood_snRNA_consist_ints$rna, function(x) transcript_names[x])
unique(Blood_snRNA_ncRNA_consist_ints$rna_name)
ggplot() +
  #geom_rect(data=Blood_donor_sites, aes(xmin = Start, xmax = End, ymin = -0.5, ymax =0), fill =  "turquoise") + 
  #geom_rect(data=Blood_acc_sites, aes(xmin = Start, xmax = End, ymin = -0.5, ymax = 0), fill = "orange") + 
  #geom_rect(data=blood_domains, aes(xmin = x_start, xmax = x_end, ymin = -Inf, ymax = Inf, fill = domain), alpha = 0.2) +
  #geom_vline(xintercept = c(249,907,1034,1043, 1854,2330, 2358, 2769, 4012, 5823, 6119, 7260), colour = "turquoise", linetype = 3) + 
  #geom_vline(xintercept = 3282, colour = "orange", linetype = 3) + 
  geom_line(data = Blood_snRNA_consist_ints, aes(x = Position, y = mean_depth, color = rna_name),alpha = 0.7) +
  #facet_wrap(~rna) + 
  scale_colour_manual(values = colors) + 
  labs(
       x = "Position on Blood",
       y = "Mean Read Depth") + 
  theme_minimal()
ggsave(paste("002COMRADES_project/output/significant_trans_RNAint/blood_snRNA_filt.pdf", sep=""), height =4, width =6.5)

#reverse interactions
Blood_snRNA_revints = Blood_snRNAints[[4]]
#Blood_snRNA_ncRNA_revints = rbind(Blood_snRNA_revints, Blood_ncRNA_revints)
Blood_snRNA_revints$rna_name = sapply(Blood_snRNA_revints$rna, function(x) transcript_names[x])
rnatypes = c(setNames(c("U2 snRNA", "U1 snRNA","U2 snRNA","U2 snRNA","U2 snRNA","U2 snRNA", "U1 snRNA","U1 snRNA", "U1 snRNA",rep("U6 snRNA", 3), "U2 snRNA","U7 snRNA"), transcript_types))
Blood_snRNA_revints$rnatype = sapply(Blood_snRNA_revints$rna_name, function(x) rnatypes[x])

ggplot() +
  #geom_rect(data=domains_df, aes(xmin = x_start, xmax = x_end, ymin = -Inf, ymax = Inf, fill = domain), alpha = 0.3) +
  #scale_fill_manual(values = domain_colours) +
  geom_line(data = Blood_snRNA_revints, aes(x = Position, y = mean_depth, color = rna_name),alpha = 0.7) +
  facet_wrap(~rnatype, scale="free") + 
  scale_colour_manual(values = colors) + 
  labs(title = "Blood snRNA interactions",
       x = "Position on snRNA",
       y = "Mean Read Depth") + 
  theme_minimal()
ggsave(paste("002COMRADES_project/output/significant_trans_RNAint/blood_snRNA_filt_rev.tiff", sep=""), height =4, width =5.5)
unique(Blood_snRNA_revints$rna)

Blood_snRNA_revints_sep = Blood_snRNAints[[4]]
ggplot(Blood_snRNA_revints_sep, aes(x=Position, y=depth, colour=sample)) + 
  geom_line() + 
  facet_wrap(~rna)

#twoway plot - need to fix rnacrosslink package!!
for(rna in RNA_ids){
  plotInteractionsAverage(cdsBlood, "AY180916_AY180916_AY180916_AY180916", rna, b="max",d="max", "C:/Users/miabe/OneDrive/Desktop/Part III Systems/Project/002COMRADES_project/output/significant_trans_RNAint/two_way_int_heatmaps")
}





### HMSbeagle
beagle_sig_ids_snRNA = beagle_q2_sig$RowNames[beagle_q2_sig$category=="snRNA"][-1]
beagle_snRNAints = consistentInts(cdsHMSbeagle, beagle_sig_ids_snRNA)


beagle_snRNA_consist_ints = beagle_snRNAints[[2]]

beagle_original_ints = beagle_snRNAints[[1]]


RNA_ids <- unique(beagle_snRNA_consist_ints$rna)
reds = brewer.pal(8, "Reds")

transcript_types <- c("U5_63Bc","U2_14B","U5_14B/38ABb","U5_23D", "U2_34ABb","U5_34A","U2_34ABa","U5_35D", "U5_14B/38ABb" ,"U2_38ABa","U5_38ABa","U6_96Aa/96Ab/96Ac","U6_96Aa/96Ab/96Ac","U6_96Aa/96Ab/96Ac","U2_34ABc", "U7")

transcript_names = setNames(transcript_types, RNA_ids)
beagle_snRNA_consist_ints$rna_name = sapply(beagle_snRNA_consist_ints$rna, function(x) transcript_names[x])

colors = c(greens[4], greens[5],greens[6],greens[7], greens[8],reds[2],reds[3],reds[5],reds[6],reds[7],reds[8],purples[8], oranges[4])


beagle_donor_sites = read.delim("002COMRADES_project/input/beagle_donor_splice_sites.txt", sep = "")
beagle_acc_sites = read.delim("002COMRADES_project/input/beagle_acc_splice_sites.txt", sep="")



beagle_seq
ggplot() +
  geom_vline(aes(xintercept = c(958,2538,3905,411),colour = "Donor_sites"),, linetype = 3) + 
  geom_vline( aes(xintercept = c(514,2333,3372,3991),colour = "Acceptor_sites"), linetype = 3) + 
  scale_colour_manual(values = c(Donor_sites = "turquoise", Acceptor_sites = "orange")) + 
  geom_line(data = beagle_snRNA_consist_ints, aes(x = Position, y = mean_depth, color = rna_name),alpha = 0.7) +
  #facet_wrap(~rna) + 
  scale_colour_manual(values = colors) + 
  labs(
       x = "Position on HMSbeagle",
       y = "Mean Read Depth") + 
  theme_minimal()
ggsave(paste("002COMRADES_project/output/significant_trans_RNAint/beagle_snRNA_filt.tiff", sep=""), height =4, width =6.5)




#reverse interactions
beagle_snRNA_revints = beagle_snRNAints[[4]]
beagle_snRNA_revints = beagle_snRNA_revints[beagle_snRNA_revints$rna %in% RNA_ids,]
beagle_snRNA_revints$rna_name = sapply(beagle_snRNA_revints$rna, function(x) transcript_names[x])
transcript_types <- c("U5_63Bc","U2_14B","U5_14B/38ABb","U5_23D", "U2_34ABb","U5_34A","U2_34ABa","U5_35D", "U5_14B/38ABb" ,"U2_38ABa","U5_38ABa","U6_96Aa/96Ab/96Ac","U6_96Aa/96Ab/96Ac","U6_96Aa/96Ab/96Ac","U2_34ABc", "U7")

rnatypes = c(setNames(c("U5","U2","U5","U5","U2", "U5","U2", "U5","U5","U2","U5","U6","U6","U6", "U2", "U7"), transcript_types))
beagle_snRNA_revints$rnatype = sapply(beagle_snRNA_revints$rna_name, function(x) rnatypes[x])

ggplot() +
  #geom_rect(data=domains_df, aes(xmin = x_start, xmax = x_end, ymin = -Inf, ymax = Inf, fill = domain), alpha = 0.3) +
  #scale_fill_manual(values = domain_colours) +
  geom_line(data = beagle_snRNA_revints, aes(x = Position, y = mean_depth, color = rna_name),alpha = 0.7) +
  facet_wrap(~rnatype, scale="free") + 
  scale_colour_manual(values = colors) + 
  labs(title = "HMSbeagle snRNA interactions",
       x = "Position on snRNA",
       y = "Mean Read Depth") + 
  theme(legend.position = "none") + 
  theme_minimal() 
ggsave(paste("002COMRADES_project/output/significant_trans_RNAint/beagle_snRNA_filt_rev.tiff", sep=""), height =4, width =5.5)

RNA_ids <- unique(beagle_snRNA_consist_ints$rna)


for(rna in RNA_ids){
  plotInteractionsAverage(cdsHMSbeagle, cdsHMSbeagle@rnas, rna, b="max",d="max", "C:/Users/miabe/OneDrive/Desktop/Part III Systems/Project/002COMRADES_project/output/significant_trans_RNAint/two_way_int_heatmaps/snRNA")
}
 #grep("FBtr0077658_FBgn0003934_23D_snRNA",cdsHMSbeagle@InputFiles$FBgn0001207_FBgn0001207_TEHMSbeagle_HMSbeagle$host$s1$V10)


#First one HMSbeagle 
print(beagle_snRNA_consist_ints[beagle_snRNA_consist_ints$Position >5500 &
                                  beagle_snRNA_consist_ints$mean_depth>0 &
                                  beagle_snRNA_consist_ints$Position <6000,],n=1000)




#U6 in LTRs 20-87
beagle_seq
RNA_ids
U6_seq = getSequence.list(annotation_fasta["FBtr0084650_FBgn0004188_96Aa_snRNA"], as.string = TRUE)[[1]][[1]]
beagle_U6_20_87 = paste(substr(beagle_seq, 0,100),"&", U6_seq, sep="")
substr(beagle_seq, 0,40)
substr(U6_seq, 30, 70)
#6811 - 6879
beagle_U6_6811_6879 = paste(substr(beagle_seq, 6790,6900),"&", U6_seq, sep="")

#Repetitive region binding - U563BC 668-799, U514B 621 799, U523D 668 799,  U534A 668-799, U5 23D 668-799, U5 35D 688 799
#U234ABb  621 959,U2-34ABa 621 959, U2 38ABa 668 799,  U2 34Abc 621 959

#U563BC 668-799, U2 38ABa 668 799

U563BC = getSequence.list(annotation_fasta["FBtr0073017_FBgn0003938_63BC_snRNA"], as.string = TRUE)[[1]][[1]]
beagle_U563bc_668_799 = paste(substr(beagle_seq, 620,860),"&", U563BC, sep="")

U238ABa = getSequence.list(annotation_fasta["FBtr0081313_FBgn0003922_38ABa_snRNA"], as.string = TRUE)[[1]][[1]]
beagle_U238ABa_668_799 = paste(substr(beagle_seq, 620,860),"&", U238ABa, sep="")

grep("gagcaacacgga", substr(beagle_seq, 664,728))
#U6 96Aa 1240 1378

beagle_U6_1240_1378 = paste(substr(beagle_seq, 1200,1410),"&", U6_seq, sep="")


#U563BC 4180-4219, U2 38ABa 4140 4218
beagle_U563bc_4180_4219 = paste(substr(beagle_seq, 4130,4270),"&", U563BC, sep="")

beagle_U238ABa_4140_4218 = paste(substr(beagle_seq, 4130,4270),"&", U238ABa, sep="")

#U2 38Aba 4520 4539

beagle_U238ABa_4520_4539 = paste(substr(beagle_seq, 4500,4560),"&", U238ABa, sep="")

#U6 5060 5091

beagle_U6_5060_5090 = paste(substr(beagle_seq, 5000,5130),"&", U6_seq, sep="")
substr(beagle_seq, 5060,5130)


#U6 5410 5419
beagle_U6_5410_5419 = paste(substr(beagle_seq, 5380,5450),"&", U6_seq, sep="")


#U7 5844 5913,
U7_seq = getSequence.list(annotation_fasta["FBtr0100848_FBgn0053504_U7_snRNA"], as.string = TRUE)[[1]][[1]]
beagle_U7_5844_5913 = paste(substr(beagle_seq, 5820,5930),"&", U7_seq, sep="")

substr(beagle_seq, 5820,5860)


beagle_U7_668_799= paste(substr(beagle_seq, 5820,5930),"&", U7_seq, sep="")

#U2 5780 5799

beagle_U2_5780_5799 = paste(substr(beagle_seq, 5750,5820),"&", U238ABa, sep="")



write.fasta(sequences = list(beagle_U238ABa_4140_4218, beagle_U238ABa_4520_4539, beagle_U238ABa_668_799, beagle_U2_5780_5799, beagle_U6_1240_1378, beagle_U6_20_87, beagle_U6_5060_5090, beagle_U6_5410_5419, beagle_U6_6811_6879, beagle_U563bc_4180_4219, beagle_U563bc_668_799, beagle_U7_5844_5913),
                             names = c("beagle_U238ABa_4140_4218", "beagle_U238ABa_4520_4539", "beagle_U238ABa_668_799", "beagle_U2_5780_5799", "beagle_U6_1240_1378", "beagle_U6_20_87", "beagle_U6_5060_5090", "beagle_U6_5410_5419", "beagle_U6_6811_6879", "beagle_U563bc_4180_4219", "beagle_U563bc_668_799", "beagle_U7_5844_5913"),
                             file.out = "beagle_snRNA_int_forfold.fasta")














#3S18
tS18_snRNA_consist_ints = tS18_snRNAints[[2]]
tS18_ncRNA_consist_ints = tS18_ncRNAints[[2]]
zero_row = tS18_ncRNA_consist_ints[tS18_ncRNA_consist_ints$Position >4000,][1,]
zero_row$Position = zero_row$Position -1
zero_row$mean_depth = 0
tS18_snRNA_ncRNA_consist_ints = rbind(tS18_snRNA_consist_ints,tS18_ncRNA_consist_ints, zero_row)

RNA_ids <- unique(tS18_snRNA_ncRNA_consist_ints$rna)

transcript_types <- c("U2_38ABb snRNA","U6_96Aa/96Ab/96Ac snRNA","U6_96Aa/96Ab/96Ac snRNA","U6_96Aa/96Ab/96Ac snRNA","Like-U snRNA", "unknown long ncRNA")
transcript_names = setNames(transcript_types, RNA_ids)
tS18_snRNA_ncRNA_consist_ints$rna_name = sapply(tS18_snRNA_ncRNA_consist_ints$rna, function(x) transcript_names[x])

colors = c("#F5C000", greens[5], purples[5], "brown")


tS18_donor_sites = read.delim("002COMRADES_project/input/tS18_donor_splice_sites.txt", sep = "")
tS18_acc_sites = read.delim("002COMRADES_project/input/tS18_acc_sites.txt", sep="")

tS18_domains

ggplot() +
  geom_vline(xintercept = c(2299,3258,4120,4251,4887,4911,5762), colour = "turquoise", linetype = 3) + 
  geom_vline(xintercept = c(72,73,2348, 3903, 4516, 5837,5838), colour = "orange", linetype = 3) + 
    geom_line(data = tS18_snRNA_ncRNA_consist_ints, aes(x = Position, y = mean_depth, color = rna_name),alpha = 0.7) +
  #facet_wrap(~rna) + 
  scale_colour_manual(values = colors) + 
  labs(title = "3S18 snRNA interactions",
       x = "Position on 3S18",
       y = "Mean Read Depth") + 
  theme_minimal()


ggsave(paste("002COMRADES_project/output/significant_trans_RNAint/tS18_snRNA_filt.tiff", sep=""), height =4, width =6.5)


#reverse interactions
tS18_snRNA_revints = tS18_snRNAints[[4]]
tS18_ncRNA_revints = tS18_ncRNAints[[4]]
tS18_snRNA_ncRNA_revints = rbind(tS18_snRNA_revints,tS18_ncRNA_revints)

tS18_snRNA_ncRNA_revints$rna_name = sapply(tS18_snRNA_ncRNA_revints$rna, function(x) transcript_names[x])
rnatypes = c(setNames(c("U2 snRNA","U6 snRNA","U6 snRNA", "U6 snRNA","likeU snRNA", "unknown lncRNA"), transcript_types))
tS18_snRNA_ncRNA_revints$rnatype = sapply(tS18_snRNA_ncRNA_revints$rna_name, function(x) rnatypes[x])
ggplot() +
  #geom_rect(data=domains_df, aes(xmin = x_start, xmax = x_end, ymin = -Inf, ymax = Inf, fill = domain), alpha = 0.3) +
  #scale_fill_manual(values = domain_colours) +
  geom_line(data = tS18_snRNA_ncRNA_revints, aes(x = Position, y = mean_depth, color = rna_name),alpha = 0.7) +
  facet_wrap(~rnatype, scale="free") + 
  scale_colour_manual(values = colors) + 
  labs(title = "3S18 snRNA interactions",
       x = "Position on 3S18",
       y = "Mean Read Depth") + 
  theme_minimal()
ggsave(paste("002COMRADES_project/output/significant_trans_RNAint/tS18_snRNA_filt_rev.tiff", sep=""), height =4, width =5.5)


for(rna in RNA_ids){
  print(rna)
  plotInteractionsAverage(cds3S18, cds3S18@rnas, rna, b="max",d="max", "C:/Users/miabe/OneDrive/Desktop/Part III Systems/Project/002COMRADES_project/output/significant_trans_RNAint/two_way_int_heatmaps/snRNA")
}

RNA_ids
plotInteractionsAverage(cds3S18, cds3S18@rnas, "FBtr0346558_FBgn0267303_CR45739_ncRNA", b="max",d="max", "C:/Users/miabe/OneDrive/Desktop/Part III Systems/Project/002COMRADES_project/output/significant_trans_RNAint/two_way_int_heatmaps/snRNA")
plotInteractionsAverage(cds3S18, cds3S18@rnas, "FBtr0081293_FBgn0003923_38ABb_snRNA" , b="max",d="max", "C:/Users/miabe/OneDrive/Desktop/Part III Systems/Project/002COMRADES_project/output/significant_trans_RNAint/two_way_int_heatmaps/snRNA")


#############
#Folding RNAs
print(tS18_snRNA_ncRNA_consist_ints[tS18_snRNA_ncRNA_consist_ints$Position >5760 &
                                      tS18_snRNA_ncRNA_consist_ints$mean_depth>0 &
                                      tS18_snRNA_ncRNA_consist_ints$Position <6000,],n=300)



RNA_ids
#LTR + unknown lncRNA
#First one 3S18 0-39, mostly last bit of the ncRNA
tS18_seq = getSequence.list(TE_fasta, as.string = TRUE)[[2]][[1]]
TE_fasta = read.fasta(file = "002COMRADES_project/input/TE_sequences.fa")

lncRNA = getSequence.list(annotation_fasta["FBtr0346558_FBgn0267303_CR45739_ncRNA"], as.string = TRUE)[[1]][[1]]

tS18_lncRNA_0_39 = paste(substr(tS18_seq, 0,49),"&", substr(lncRNA, 343,370), sep="")
gtgacatttacagg
#Second LTR

tS18_lncRNA_5766_5813 = paste(substr(tS18_seq, 5767,5867),"&", substr(lncRNA, 290,400), sep="")

#Second LTR + U6
U6_snRNA = getSequence.list(annotation_fasta["FBtr0084650_FBgn0004188_96Aa_snRNA"], as.string = TRUE)[[1]][[1]]
tS18_U6_5780_5793 = paste(substr(tS18_seq, 5756,5813),"&", U6_snRNA, sep="")


#U6 snRNA - 2560, 2639

tS18_U6_2560_2639 = paste(substr(tS18_seq, 2550,2629),"&", U6_snRNA, sep="")
substr(U6_snRNA, 40,107)

#U6 snRNA 3500 3559

tS18_U6_3500_3559 = paste(substr(tS18_seq, 3490,3569),"&", U6_snRNA, sep="")


#LikeU 1560-1619

tS18_LU_1560_1619 = paste(substr(tS18_seq, 1550,1629),"&", likeU_snRNA, sep="")


#LikeU 2700 2736
likeU_snRNA = getSequence.list(annotation_fasta["FBtr0452214_FBgn0263847_LU_snRNA"], as.string = TRUE)[[1]][[1]]

tS18_LU_2700_2736 = paste(substr(tS18_seq, 2690,2746),"&", likeU_snRNA, sep="")

#likeU 5020-5139

likeU_snRNA = getSequence.list(annotation_fasta["FBtr0452214_FBgn0263847_LU_snRNA"], as.string = TRUE)[[1]][[1]]

tS18_LU_5020_5139 = paste(substr(tS18_seq, 5010,5149),"&", likeU_snRNA, sep="")

grep("attgccc", substr(tS18_seq, 5010,5070))

#U2 snRNA 720 799
U2_snRNA = getSequence.list(annotation_fasta["FBtr0081293_FBgn0003923_38ABb_snRNA"], as.string = TRUE)[[1]][[1]]
tS18_U2_720_799 = paste(substr(tS18_seq, 710,809),"&", substr(U2_snRNA,0,110), sep="")

#U2 snRNA 1740 1759
tS18_U2_1740_1759 = paste(substr(tS18_seq, 1730,1769),"&", U2_snRNA, sep="")



#U2 snRNA 4580 4738
transcript_names

U2_snRNA = getSequence.list(annotation_fasta["FBtr0081293_FBgn0003923_38ABb_snRNA"], as.string = TRUE)[[1]][[1]]
tS18_U2_4580_4738 = paste(substr(tS18_seq, 4570,4748),"&", U2_snRNA, sep="")

substr(U2_snRNA, 110, 135)

write.fasta(sequences = list(tS18_lncRNA_0_39, tS18_lncRNA_5766_5813, tS18_U6_2560_2639, tS18_U6_3500_3559, tS18_U6_5780_5793, tS18_LU_1560_1619, tS18_LU_2700_2736, tS18_LU_5020_5139, tS18_U2_720_799, tS18_U2_1740_1759, tS18_U2_4580_4738),
            names = c("tS18_lncRNA_0_39", "tS18_lncRNA_5766_5813", "tS18_U6_2560_2639", "tS18_U6_3500_3559", "tS18_U6_5780_5793", "tS18_LU_1560_1619", "tS18_LU_2700_2736", "tS18_LU_5020_5139", "tS18_U2_720_799", "tS18_U2_1740_1759", "tS18_U2_4580_4738"),
            file.out = "tS18_snRNA_int_forfold.fasta")


ggplot() +
  geom_vline(xintercept = c(2299,3258,4120,4251,4887,4911,5762), colour = "turquoise", linetype = 3) + 
  geom_vline(xintercept = c(72,73,2348, 3903, 4516, 5837,5838), colour = "orange", linetype = 3) + 
  


####snoRNA####


### Blood


### Blood

#forward interactions
Blood_snoRNA_consist_ints = Blood_snoRNAints[[2]]

ggplot() +
  geom_line(data=Blood_snoRNAints[[1]], aes(x=Position, y=depth))

#No consistent ones



### HMSbeagle
beagle_snoRNA_consist_ints = beagle_snoRNAints[[2]]

ggplot() +
  geom_line(data=beagle_snoRNAints[[1]], aes(x=Position, y=depth))

#No consistent ones


#3S18
tS18_snoRNA_consist_ints = tS18_snoRNAints[[2]]

ggplot() +
  geom_line(data=beagle_snoRNAints[[1]], aes(x=Position, y=depth))

#no significant ones




#############
####ncRNA####
#############

### Blood

#forward interactions
Blood_ncRNA_consist_ints = Blood_ncRNAints[[2]]
Blood_ncRNA_original_ints = Blood_ncRNAints[[1]]
ggplot(Blood_ncRNA_consist_ints, aes(x=Position, y=mean_depth, colour=rna)) + 
  geom_line() + facet_wrap(~rna)
RNA_ids <- unique(Blood_ncRNA_consist_ints$rna)


# Corresponding transcript types
Blood_ncRNA_consist_ints = Blood_ncRNA_consist_ints[Blood_ncRNA_consist_ints$rna %in% RNA_ids,]
transcript_types <- c("Uhg4_ncRNA")

reds = brewer.pal(n=8, name= "Reds")
colors <- c(
  #"7SL_signalling_particle_1" = reds[3],
  "7SL_signalling_particle" = reds[5])

transcript_colours = c(setNames(colors, RNA_ids))
transcript_names = setNames(transcript_types, RNA_ids)
Blood_ncRNA_consist_ints$rna_name = sapply(Blood_ncRNA_consist_ints$rna, function(x) transcript_names[x])

ggplot() +
  #geom_rect(data=blood_domains, aes(xmin = x_start, xmax = x_end, ymin = -Inf, ymax = Inf, fill = domain), alpha = 0.3) +
  #scale_fill_manual(values = structure(blood_domains$colours, names = blood_domains$domain)) +
  geom_line(data = subset(Blood_ncRNA_consist_ints, Position <600 &Position >100), aes(x = Position, y = mean_depth, color = rna_name),alpha = 0.7) +
  #facet_wrap(~rna) + 
  scale_colour_manual(values = colors) + 
  labs(title = "Blood ncRNA interactions",
       x = "Position on Blood",
       y = "Mean Read Depth") + 
  theme_minimal()
ggsave(paste("002COMRADES_project/output/significant_trans_RNAint/blood_ncRNA_filt.tiff", sep=""), height =4, width =5.5)

#reverse interactions
Blood_ncRNA_revints = Blood_ncRNAints[[4]]
Blood_ncRNA_revints = Blood_ncRNA_revints[Blood_ncRNA_revints$rna %in% RNA_ids,]
Blood_ncRNA_revints$rna_name = sapply(Blood_ncRNA_revints$rna, function(x) transcript_names[x])
rnatypes = c(setNames(c("7SL_signalling_particle","7SL_signalling_particle"), transcript_types))
Blood_ncRNA_revints$rnatype = sapply(Blood_ncRNA_revints$rna_name, function(x) rnatypes[x])

ggplot() +
  #geom_rect(data=domains_df, aes(xmin = x_start, xmax = x_end, ymin = -Inf, ymax = Inf, fill = domain), alpha = 0.3) +
  #scale_fill_manual(values = domain_colours) +
  geom_line(data = Blood_ncRNA_revints, aes(x = Position, y = mean_depth, color = rna_name),alpha = 0.7) +
  facet_wrap(~rnatype, scale="free") + 
  scale_colour_manual(values = colors) + 
  labs(title = "Blood ncRNA interactions",
       x = "Position on Blood",
       y = "Mean Read Depth") + 
  theme_minimal()
ggsave(paste("002COMRADES_project/output/significant_trans_RNAint/blood_ncRNA_filt_rev.tiff", sep=""), height =3, width =4.5)


#twoway plot - need to fix rnacrosslink package!!
for(rna in RNA_ids){
  plotInteractionsAverage(cdsBlood, "AY180916_AY180916_AY180916_AY180916", RNA_ids, b="max",d="max", "002COMRADES_project/output/significant_trans_RNAint/two_way_int_heatmaps/ncRNA/")
}


### HMSbeagle
beagle_ncRNA_consist_ints = beagle_ncRNAints[[2]]
beagle_ncRNA_original_ints = beagle_ncRNAints[[1]]

ggplot(beagle_ncRNA_original_ints, aes(x = Position, y=depth, colour = sample)) + 
  geom_line()
RNA_ids <- unique(beagle_ncRNA_consist_ints$rna)

ggplot() +
  geom_rect(data=beagle_domains, aes(xmin = x_start, xmax = x_end, ymin = -Inf, ymax = Inf, fill = domain), alpha = 0.3) +
  scale_fill_manual(values = structure(blood_domains$colours, names = blood_domains$domain)) +
  geom_line(data = beagle_ncRNA_consist_ints, aes(x = Position, y = mean_depth, color = rna),alpha = 0.7) +
  #facet_wrap(~rna) + 
  scale_colour_manual(values = colors) + 
  labs(title = "Blood ncRNA interactions",
       x = "Position on Blood",
       y = "Mean Read Depth") + 
  theme_minimal()
ggsave(paste("002COMRADES_project/output/significant_trans_RNAint/blood_ncRNA_filt.tiff", sep=""), height =4, width =5.5)

 
#3S18
tS18_ncRNA_consist_ints = tS18_ncRNAints[[2]]
tS18_ncRNA_original_ints = tS18_ncRNAints[[1]]

ggplot(tS18_ncRNA_original_ints, aes(x = Position, y=depth, colour = sample)) + 
  geom_line()



#LTR1 
ggplot() +
  #geom_rect(data=tS18_domains, aes(xmin = x_start, xmax = x_end, ymin = -Inf, ymax = Inf, fill = domain), alpha = 0.3) +
  #scale_fill_manual(values = structure(blood_domains$colours, names = blood_domains$domain)) +
  geom_line(data = subset(tS18_ncRNA_consist_ints, Position >4500), aes(x = Position, y = mean_depth, color = rna),alpha = 0.7) +
  #facet_wrap(~rna) + 
  #scale_colour_manual(values = colors) + 
  labs(title = "3S18 ncRNA interactions",
       x = "Position on 3S18",
       y = "Mean Read Depth") + 
  theme_minimal()

ggsave(paste("002COMRADES_project/output/significant_trans_RNAint/blood_ncRNA_filt.tiff", sep=""), height =4, width =5.5)

#LTR2



#############
####miRNA####
#############

### Blood

#forward interactions
Blood_miRNA_consist_ints = Blood_miRNAints[[2]]

ggplot() +
  geom_line(data=Blood_miRNAints[[1]], aes(x=Position, y=depth, color=sample))

#No consistent ones



### HMSbeagle
beagle_miRNA_consist_ints = beagle_miRNAints[[2]]

ggplot() +
  geom_line(data=beagle_snoRNAints[[1]], aes(x=Position, y=depth))

#No consistent ones

### 3S18
#no significant ones




##########################################################################################################
####                   EXTRACTING THE SEQUENCES OF INTERACTION LOCATION FOR FOLDING                   ####
##########################################################################################################
getIntSeqs = function(TE, consist_ints, consistent_rints){
  # Remove 0 positions to get ranges
  consist_ints = consist_ints[order(consist_ints$Position),]
  consist_ints = consist_ints[order(consist_ints$rna),]
  int_peaks = consist_ints[consist_ints$mean_depth != 0,]
  
  # Forward sequence intervals
  int_peaks_filt_forward <- int_peaks %>%
    arrange(rna, Position) %>%
    group_by(rna) %>%
    mutate(diff = Position - lag(Position, default = first(Position)),
           group_id = cumsum(diff > 2)) %>%
    group_by(rna, group_id) %>%
    summarise(start = min(Position),
              end = max(Position),
              total_depth = sum(mean_depth),
              rna_name = rna_name[1],
              .groups = "drop") %>%
    mutate(length = end - start) %>%
    mutate(start2 = start - 10,
           end2 = end + 10) %>%
    ungroup()
  
  # Reverse sequence intervals
  consistent_rints = consistent_rints[consistent_rints$mean_depth >0,]
  int_peaks_filt_rev <- consistent_rints %>%
    arrange(rna, Position) %>%
    group_by(rna) %>%
    mutate(diff = Position - lag(Position, default = first(Position)),
           group_id = cumsum(diff > 2)) %>%
    group_by(rna, group_id) %>%
    summarise(start_rev = min(Position),
              end_rev = max(Position),
              total_depth_rev = sum(mean_depth),
              rna_name_rev = rna_name[1],
              .groups = "drop") %>%
    mutate(length_rev = end_rev - start_rev) %>%
    mutate(start2_rev = ifelse(start_rev > 10, start_rev - 10, 1),
           end2_rev = end_rev + 10) %>%
    ungroup()
  
  int_peaks_forward_rev = left_join(int_peaks_filt_forward, int_peaks_filt_rev, by = "rna")
  
  # Get forward sequences (TE)
  if(TE == "Blood"){
    fasta_number = 3
  }
  if(TE == "HMSbeagle"){
    fasta_number = 1
  }
  if(TE == "tS18"){
    fasta_number = 2
  }
  
  seq = getSequence.list(TE_fasta, as.string = TRUE)[[fasta_number]][[1]]
  
  # Get reverse sequences from fasta file
  transcripts_peaks = int_peaks_forward_rev$rna
  sequence_rna = c()
  sequence_TE = c()
  
  for(i in 1:nrow(int_peaks_forward_rev)){
    # Forward sequence
    sequence_TE = c(sequence_TE, substr(seq, int_peaks_forward_rev$start[i], int_peaks_forward_rev$end[i]))
    
    # Reverse sequence
    current_transcript = transcripts_peaks[i]
    current_seq = getSequence(annotation_fasta[[current_transcript]], as.string = TRUE)[[1]]
    print(current_seq)
    current_seq_rna = substr(current_seq, int_peaks_forward_rev$start_rev[i], int_peaks_forward_rev$end_rev[i])
    #print(current_seq_rna)
    sequence_rna = c(sequence_rna, current_seq_rna)
  }
  
  int_peaks_forward_rev_seq = int_peaks_forward_rev[, c(1, 3, 4, 5, 6, 7, 11, 12, 13)]
  int_peaks_forward_rev_seq$TE_seq = sequence_TE
  int_peaks_forward_rev_seq$rna_seq = sequence_rna
  
  # Creating RNA sequences for folding
  combined_RNA_seqs = c()
  for(i in 1:nrow(int_peaks_forward_rev_seq)){
    current_combined_RNA_seqs = paste(int_peaks_forward_rev_seq$TE_seq[i],
                                      #strrep("a", 50),
                                      "&",
                                      int_peaks_forward_rev_seq$rna_seq[i],
                                      sep = "")
    combined_RNA_seqs = c(combined_RNA_seqs, current_combined_RNA_seqs)
  }
  
  int_peaks_forward_rev_seq$combined_rna_seq = combined_RNA_seqs
  int_peaks_forward_rev_seq$combined_rna_name = paste(int_peaks_forward_rev_seq$rna, int_peaks_forward_rev_seq$start, int_peaks_forward_rev_seq$end, sep="")
  
  return(int_peaks_forward_rev_seq)
}

#Reading in the FASTA file
TE_fasta = read.fasta(file = "002COMRADES_project/input/TE_sequences.fa")
annotation_fasta = read.fasta(file = "002COMRADES_project/input/Dmel_hyb_full_TE_cdna_ncrna.fasta")


#Creating a GRanges class of the TE

TE_sizes = c(cdsHMSbeagle@rnaSize, cds3S18@rnaSize, cdsBlood@rnaSize)



# snRNA

View(blood_snRNA_int_seqs)
blood_snRNA_int_seqs = getIntSeqs("Blood", Blood_snRNA_consist_ints, Blood_snRNA_revints)
print(blood_snRNA_int_seqs,n = 22)

ggplot(subset(Blood_snRNA_consist_ints, Position<1300 &Position >800) , aes(x=Position, y=mean_depth, colour=rna_name)) + geom_line() + geom_vline(xintercept = c(907,1034,1043) )
#extracting peaks: 
#Blood 910-1169

substr(blood_seq, 1034,1043)

#Blood-U1 95Cc
substr(blood_seq,910, 1169)
RNA_ids
U1_95Cc_seq = getSequence(annotation_fasta[["FBtr0084487_FBgn0004187_95Cc_snRNA"]], as.string = TRUE)[[1]]
U1_95CC_Blood_1000 = paste(substr(blood_seq,910, 1169), "&",substr(U1_95Cc_seq, 1,100), sep="")

#Blood-U1 21D/95Ca/cb
U1_21_seq = getSequence(annotation_fasta[["FBtr0078028_FBgn0003916_21D_snRNA"]], as.string = TRUE)[[1]]

U1_21D_Blood_1000 = paste(substr(blood_seq,910, 1169), "&",U1_21_seq, sep="")
substr(U1_21_seq,160,170)

#Blood-U7
U7_seq = getSequence(annotation_fasta[["FBtr0100848_FBgn0053504_U7_snRNA"]], as.string = TRUE)[[1]]
U7_Blood_1000 = paste(substr(blood_seq,910, 1169), "&",substr(U7_seq,10,50), sep="")

#Blood 3300 region, acceptor

ggplot(subset(Blood_snRNA_consist_ints, Position<3800 &Position >3100) , aes(x=Position, y=mean_depth, colour=rna_name)) + geom_line() + geom_vline(xintercept = c(3282) )
print(subset(Blood_snRNA_consist_ints, Position<6000 &Position >5000 &mean_depth >0), n=400)


#3250 - 3350

U1_95CC_Blood_3300 = paste(substr(blood_seq,3250, 3350), "&",U1_95Cc_seq, sep="")

U1_21D_Blood_3300 = paste(substr(blood_seq,3250, 3350), "&",U1_21_seq, sep="")


#Non overlapping highest region, 3470 -3770
print(subset(Blood_snRNA_consist_ints, Position<3800 &Position >3400 &mean_depth >0 &rna =="FBtr0074208_FBgn0003920_14B_snRNA"), n=400)
U2_14B_seq = getSequence(annotation_fasta[["FBtr0074208_FBgn0003920_14B_snRNA"]], as.string = TRUE)[[1]]
U2_34Abb_seq = getSequence(annotation_fasta[["FBtr0080443_FBgn0004192_34ABb_snRNA"]], as.string = TRUE)[[1]]
U2_38Abb_seq = getSequence(annotation_fasta[["FBtr0081293_FBgn0003923_38ABb_snRNA"]], as.string = TRUE)[[1]]

U2_14B_blod_3500 = paste(substr(blood_seq,3470, 3770), "&",substr(U2_14B_seq,1,120), sep="")
U2_34ABb_blood_3500 = paste(substr(blood_seq,3470, 3770),"&",substr(U2_34Abb_seq,1,100), sep="")
U2_38ABb_blood_3500 = paste(substr(blood_seq,3470, 3770), "&",substr(U2_38Abb_seq,1,130), sep="")

write.fasta(sequences = list(U1_95CC_Blood_1000, U1_21D_Blood_1000, U7_Blood_1000, U1_95CC_Blood_3300, U1_21D_Blood_3300, U2_14B_blod_3500,U2_34ABb_blood_3500, U2_38ABb_blood_3500),
            names = c("U1_95CC_Blood_1000", "U1_21D_Blood_1000", "U7_Blood_1000","U1_95CC_Blood_3300", "U1_21D_Blood_3300","U2_14B_blod_3500","U2_34ABb_blood_3500","U2_38ABb_blood_3500"),
            file.out = "Blood_snRNA_int_forfold.fasta")


grep("gccgatagc", substr(blood_seq,3650, 3670))
#3S18


ggplot(subset(tS18_snRNA_consist_ints, Position<3000 &Position >2300) , aes(x=Position, y=mean_depth, colour=rna_name)) + geom_line() 

#U2 5240 - 5335

U2_34ABb_blood_5240_5335 = paste(substr(blood_seq,5210, 5345),"&",substr(U2_34Abb_seq,1,100), sep="")

#U7 LTR 7300 - 7359 
print(subset(Blood_snRNA_consist_ints, Position>7000 &mean_depth>0), n=12000)

U7_Blood_7300_7359 = paste(substr(blood_seq,7290, 7369), "&",U7_seq, sep="")

#U6 6440-6710
print(subset(Blood_snRNA_consist_ints, Position<7000 &Position >6000 &mean_depth >0 &rna =="FBtr0084650_FBgn0004188_96Aa_snRNA"), n=400)
U6_seq = getSequence(annotation_fasta[["FBtr0084650_FBgn0004188_96Aa_snRNA"]], as.string = TRUE)[[1]]
U6_Blood_6440_6710 = paste(substr(blood_seq,6430, 6720), "&",U6_seq, sep="")


write.fasta(sequences = list(U2_34ABb_blood_5240_5335,U7_Blood_7300_7359, U6_Blood_6440_6710),
            names = c("U2_34ABb_blood_5240_5335","U7_Blood_7300_7359", "U6_Blood_6440_6710"),
            file.out = "blood_snRNA_int_forfold2.fasta")




#HMSbeagle


ggplot(subset(Blood_snRNA_consist_ints, Position<1300 &Position >800) , aes(x=Position, y=mean_depth, colour=rna_name)) + geom_line() + geom_vline(xintercept = c(907,1034,1043) )



























