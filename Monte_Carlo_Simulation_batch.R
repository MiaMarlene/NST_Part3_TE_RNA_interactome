

library(dplyr, lib.loc= "/mnt/beegfs6/home1/miska/mb2453/R/x86_64-pc-linux-gnu-library/4.2")



#Using command line variables
TE_type_sim <- as.character(commandArgs(trailingOnly = TRUE))
TE_type = substr(TE_type_sim, 1, nchar(TE_type_sim) -1)
print(TE_type_sim)
#Importing all the files
input = "/mnt/home3/miska/mb2453/002COMRADES_project/input/MonteCarloSim/"
output = "/mnt/home3/miska/mb2453/002COMRADES_project/output/MonteCarloSim/"

real_inter_df = read.csv(paste(input,"real_inter_", TE_type, ".csv", sep=""), row.names=1)
colnames(real_inter_df) = c("RowNames", "sum")

total_intras_df = read.csv(paste(input,TE_type,"_total_distr.csv", sep=""), row.names = 1)
colnames(total_intras_df) = c("RowNames", "Freq")
#head(total_intras_df)
total_distr = c(rep(total_intras_df$RowNames, total_intras_df$Freq))

sim_higher_than_observed = data.frame(higher_pval = c(rep(0,nrow(real_inter_df))), ocurrences = c(rep(0,nrow(real_inter_df))))
sim_higher_than_observed$RowNames = real_inter_df$RowNames


#batch_df = data.frame(batch1 = c(rep(0,nrow(real_inter_df))), batch2 = c(rep(0,nrow(real_inter_df))), batch3 = c(rep(0,nrow(real_inter_df))), batch4 = c(rep(0,nrow(real_inter_df))), batch5 = c(rep(0,nrow(real_inter_df))), batch6 = c(rep(0,nrow(real_inter_df))), batch7 = c(rep(0,nrow(real_inter_df))), batch8 = c(rep(0,nrow(real_inter_df))), batch9 = c(rep(0,nrow(real_inter_df))), batch10 = c(rep(0,nrow(real_inter_df))))

#batch_pval = c(rep(0, nrow(real_inter_df)))

#head(sim_higher_than_observed)
#head(total_intras_df)
#head(real_inter_df)

if(TE_type == "blood"){tot = 1138739}
if(TE_type == "beagle"){tot =1246076}
if(TE_type == "tS18"){tot =1075812}


#print(tot)
#Running the MC simulation
N=1000
M = 10

for(m in 1:M){
  TE_type_sim = paste(TE_type_sim, m, sep="")
  batch_pval = c(rep(0, nrow(real_inter_df)))
  sim_higher_than_observed = data.frame(higher_pval = c(rep(0,nrow(real_inter_df))), ocurrences = c(rep(0,nrow(real_inter_df))))
  sim_higher_than_observed$RowNames = real_inter_df$RowNames
  
  for(n in 1:N){
    sim_interactions = sample(total_distr, tot) #From the distribution of all samples present in this dataset, sample n times, where n is the total number of chimera with blood in the original dataset
    sim_table = as.data.frame(table(as.character(sim_interactions))) #Create a table of the simulated interactions
    colnames(sim_table) = c("RowNames", "Freq")
  #compare simulated interactions to bloodxtranscript interactions
   
    comparing_df = left_join(real_inter_df,sim_table, by="RowNames")
    comparing_df[is.na(comparing_df)] = 0
    comparing_df$higher = ifelse(comparing_df$Freq >= comparing_df$sum, 1, 0)
    sim_higher_than_observed$higher_pval = sim_higher_than_observed$higher_pval + comparing_df$higher
    sim_higher_than_observed$ocurrences = sim_higher_than_observed$ocurrences + comparing_df$Freq 
  #batch_pval = batch_pval + sim_higher_than_observed$higher_pval
#  if(n %%1000== 0){
 #   batch = as.integer(n/1000)
  #  batch_df[,batch] = batch_pval
   # batch_pval = c(rep(0, nrow(real_inter_df)))
    }

  write.csv(sim_higher_than_observed, paste(output, TE_type_sim, "batch_pvalues.csv", sep=""))
} 
#batch = as.integer(n/1000)
#batch_df[,batch] = batch_pval
#batch_pval = c(rep(0, nrow(real_inter_df)))

#write.csv(sim_higher_than_observed, paste(output, TE_type_sim,"output_MC_simulation.csv", sep=""))
#write.csv(batch_pval, paste(output, TE_type_sim, "batch_pvalues.csv", sep=""))

#write.csv(sim_table, paste(output, TE_type,"single_sim.csv",sep=""))

print("finished simulation")
