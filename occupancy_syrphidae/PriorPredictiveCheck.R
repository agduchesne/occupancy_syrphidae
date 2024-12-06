# Script to analyze and plot the Prior Predictive Check.
# Recommended to run in a separate session from the main script, as this script uses many of the same object names.
# Requires running the Prior Predictive Check model from the main script. See occupancy_syrphidae_main.R
# Please set working directory to occupancy_syrphidae folder


# Load Prior Predictive Check Model Summary and species list
load("FinalModel/splist.RData")
load("PriorPredCheck/PriorCheck.ModelSummary.RData")


# Libraries
library(dplyr)
library(ggplot2)
library(scales)
library(grid)
library(gridExtra)
library(rjags)
library(R2jags)
library(runjags)


#load("PriorPredCheck/PriorCheck.RData")
# Create Model Summary, if not already done
#PriorCheck.ModelSummary <- summary(PriorCheck$jags.out)
# Save the model summary on your computer
#save(PriorCheck.ModelSummary,file="PriorPredCheck/PriorCheck.ModelSummary.RData")  #SAVE


# Create column in ModelSummary that lists the species name associated with each parameter.
PriorCheck.ModelSummary <- mutate(data.frame(PriorCheck.ModelSummary), Species = NA, .before = "Lower95")
for(i in 1:nrow(PriorCheck.ModelSummary)){
  if(grepl("psi.sp|p.sp|psi.era", rownames(PriorCheck.ModelSummary[i,]))){ #If the parameter is any of the following parameters, find the associated species name: psi.sp[sp], p.sp[sp], psi.era[sp,era]
    if(grepl("psi.sp|p.sp", rownames(PriorCheck.ModelSummary[i,]))){ #if it is psi.sp or p.sp, the species number will be within square brackets
      PriorCheck.ModelSummary$Species[i] <- splist[as.numeric(regmatches(rownames(PriorCheck.ModelSummary[i,]), regexec("\\[(\\d+)\\]", rownames(PriorCheck.ModelSummary[i,])))[[1]][2])]} #regexec finds the position of the parts of the string that satisfy the condition (i.e., are within square brackets). regmatches then obtains the substrings at those positions. We only want the first match (there should only be one), so we use [[1]]. We want the number only, without the square brackets, so we use [2].
    if(grepl("psi.era", rownames(PriorCheck.ModelSummary[i,]))){ #if it is psi.era, the species number will be between a square bracket and a comma
      PriorCheck.ModelSummary$Species[i] <- splist[as.numeric(regmatches(rownames(PriorCheck.ModelSummary[i,]), regexec("\\[(\\d+)\\,", rownames(PriorCheck.ModelSummary[i,])))[[1]][2])]}}} #regexec finds the position of the parts of the string that satisfy the condition (i.e., are within square brackets). regmatches then obtains the substrings at those positions. We only want the first match (there should only be one), so we use [[1]]. We want the number only, without the square brackets, so we use [2].

PriorCheck.MCMC <- PriorCheck$jags.out$mcmc #from here on, use MCMC as shorthand for quick/easy access to the mcmc samples in the model output




################################
#### Estimates of Occupancy ####
################################

# Function to compute occupancy (as a mcmc.list object containing each sampled chain) for a given species and era
compute_occupancy <- function(sp,era){
  #Occupancy is baseline occupancy + effect of species + effect of era
  mu.psi <- matrix(as.matrix(PriorCheck.MCMC[,"mu.psi"]), nrow=PriorCheck$jags.out$sample, ncol=length(PriorCheck.MCMC)) #Obtain the necessary parameter distributions
  psi.sp <- matrix(as.matrix(PriorCheck.MCMC[,paste("psi.sp[",sp,"]",sep="")]), nrow=PriorCheck$jags.out$sample, ncol=length(PriorCheck.MCMC))
  psi.era <- matrix(as.matrix(PriorCheck.MCMC[,paste("psi.era[",sp,",",era,"]",sep="")]), nrow=PriorCheck$jags.out$sample, ncol=length(PriorCheck.MCMC))
  occupancy <- as.matrix(plogis(as.matrix(mu.psi+psi.sp+psi.era))) #Compute occupancy
  occupancy <- as.mcmc.list(lapply(seq_len(ncol(occupancy)), function(i) as.mcmc(occupancy[,i]))) #convert matrix back into mcmc object (each chain i is converted into an mcmc object by as.mcmc(), which are combined into an mcmc.list object)
  return(occupancy)}

# Function to compute across-species mean occupancy in each era
compute_cross.sp_occupancy <- function(era){
  mu.psi <- matrix(as.matrix(PriorCheck.MCMC[,"mu.psi"]), nrow=PriorCheck$jags.out$sample, ncol=length(PriorCheck.MCMC)) #Obtain mu.psi distribution
  if (era==1){ #If era is 1, across species occupancy is simply the expit of mu.psi. The effect of era 1 (the "default" era) is already contained in mu.psi.
    mu.psi.era <- 0
  }else{ #If era is not 1, must obtain the mu.psi.era parameter (which represents the mean effect of this era across species)
    mu.psi.era <- matrix(as.matrix(PriorCheck.MCMC[,paste("mu.psi.era[",era,"]",sep="")]), nrow=PriorCheck$jags.out$sample, ncol=length(PriorCheck.MCMC))}
  occupancy <- as.matrix(plogis(as.matrix(mu.psi+mu.psi.era))) #Compute occupancy
  occupancy <- as.mcmc.list(lapply(seq_len(ncol(occupancy)), function(i) as.mcmc(occupancy[,i]))) #convert matrix back into mcmc object
  return(occupancy)}

# Create a list object containing all occupancy posteriors
Occupancy_list <- list()
meanoccupancies <- list()
for (era in 1:PriorCheck$jags.data$nera){ #Compute cross-species occupancies for each era
  meanoccupancies[[paste(era)]] <- compute_cross.sp_occupancy(era)}
Occupancy_list[["cross.sp.mean"]] <- meanoccupancies #add the cross-species occupancy posteriors to the master list
rm(meanoccupancies)
for (sp in 1:PriorCheck$jags.data$nsp) { #Compute species-specific occupancies
  spoccupancies <- list()
  for (era in 1:PriorCheck$jags.data$nera) {
    spoccupancies[[paste(era)]] <- compute_occupancy(sp,era)}
  Occupancy_list[[splist[sp]]] <- spoccupancies #add this species' occupancy posteriors to the master list
  rm(spoccupancies)}

# Create a dataframe of summary statistics (mean, sd, 95% intervals) for all occupancy estimates.
get_occ_stats <- function(row){ #Function to get the summary stats for each era's occupancy for a given species
  dataframe <- matrix(nrow=PriorCheck$jags.data$nera,ncol=6) #create matrix to hold data.
  for (era in 1:PriorCheck$jags.data$nera){ #Get the species name, era name, lower95, mean, upper95, and SD.
    dataframe[era,] <- c(names(Occupancy_list[row]), names(Occupancy_list[[row]][era]), HPDinterval(as.mcmc(unlist(Occupancy_list[[row]][[era]])), prob=0.95)[[1]] , summary(Occupancy_list[[row]][[era]])$statistics[[1]], HPDinterval(as.mcmc(unlist(Occupancy_list[[row]][[era]])), prob=0.95)[[2]] , summary(Occupancy_list[[row]][[era]])$statistics[[2]])}
  return(dataframe)}
OccupancySummary <- as.data.frame(do.call(rbind,lapply(1:length(Occupancy_list), get_occ_stats))) #Create the dataframe
colnames(OccupancySummary) <- c("Species","Era","Lower95","Mean","Upper95","SD")

# Make sure columns are treated as numeric values
OccupancySummary$Lower95 <- as.numeric(OccupancySummary$Lower95)
OccupancySummary$Mean <- as.numeric(OccupancySummary$Mean)
OccupancySummary$Upper95 <- as.numeric(OccupancySummary$Upper95)
OccupancySummary$SD <- as.numeric(OccupancySummary$SD)

# Save Occupancy Summary
#save(OccupancySummary,file="PriorPredCheck/PriorCheck.OccupancySummary.RData")
# Load Occupancy Summary
load("PriorPredCheck/PriorCheck.OccupancySummary.RData")






#############################################
#### Occupancy Changes (delta occupancy) ####
#############################################

# Function to compute the raw occupancy change between two eras
compute_delta.occ <- function(row,era1,era2) { 
  Occ1 <- matrix(as.matrix(Occupancy_list[[row]][[era1]]), nrow=PriorCheck$jags.out$sample, ncol=length(PriorCheck.MCMC)) #Occupancy estimate of era1
  Occ2 <- matrix(as.matrix(Occupancy_list[[row]][[era2]]), nrow=PriorCheck$jags.out$sample, ncol=length(PriorCheck.MCMC)) #Occupancy estimate of era2
  delta.occ <- as.matrix(Occ1-Occ2) #subtract the distributions
  delta.occ <- as.mcmc.list(lapply(seq_len(ncol(delta.occ)), function(i) as.mcmc(delta.occ[,i]))) 
  return(delta.occ)}

# Create a list object containing all occupancy changes between eras (delta occupancy)
delta.occ_list <- list()
for (row in 1:length(Occupancy_list)) {
  occ.changes <- list()
  for (era in PriorCheck$jags.data$nera:2) {
    for (nextera in (era-1):1) {
      occ.changes[[paste(era,"-",nextera, sep="")]] <- compute_delta.occ(row,era,nextera)}} #Occupancy change is computed between the era in question and every era below it. Example with 4 eras: 4-3,4-2,4-1,3-2,3-1,2-1
  delta.occ_list[[names(Occupancy_list[row])]] <- occ.changes #add this species' contrasts to the master list
  rm(occ.changes,row,era,nextera)}

# Create a dataframe of summary statistics (mean, sd, 95% intervals) for all occupancy change estimates.
get_delta_occ_stats <- function(row){ #Function to get the summary stats for all occupancy changes for a given species
  dataframe <- matrix(nrow=PriorCheck$jags.data$nera*(PriorCheck$jags.data$nera-1)/2,ncol=6) #create matrix to hold data. The number of pairwise contrasts is always nera*(nera-1)/2
  for (contrast in 1:(PriorCheck$jags.data$nera*(PriorCheck$jags.data$nera-1)/2)){ #For each contrast, get the species name, contrast name, lower95, mean, upper95, and SD.
    dataframe[contrast,] <- c(names(Occupancy_list[row]), names(delta.occ_list[[row]][contrast]), HPDinterval(as.mcmc(unlist(delta.occ_list[[row]][[contrast]])), prob=0.95)[[1]] , summary(delta.occ_list[[row]][[contrast]])$statistics[[1]], HPDinterval(as.mcmc(unlist(delta.occ_list[[row]][[contrast]])), prob=0.95)[[2]] , summary(delta.occ_list[[row]][[contrast]])$statistics[[2]])}
  return(dataframe)}
delta.occSummary <- as.data.frame(do.call(rbind,lapply(1:length(Occupancy_list), get_delta_occ_stats))) #Create the dataframe
colnames(delta.occSummary) <- c("Species","EraContrast","Lower95","Mean","Upper95","SD")

# Make sure columns are treated as numeric values
delta.occSummary$Lower95 <- as.numeric(delta.occSummary$Lower95)
delta.occSummary$Mean <- as.numeric(delta.occSummary$Mean)
delta.occSummary$Upper95 <- as.numeric(delta.occSummary$Upper95)
delta.occSummary$SD <- as.numeric(delta.occSummary$SD)

# Save delta Occupancy Summary
#save(delta.occSummary,file="PriorPredCheck/PriorCheck.delta.occSummary.RData")
# Load delta Occupancy Summary
load("PriorPredCheck/PriorCheck.delta.occSummary.RData")






##################
### PLOT GRAPH ###
##################

load("PriorPredCheck/PriorCheck.OccupancySummary.RData")
load("PriorPredCheck/PriorCheck.delta.occSummary.RData")


Begin_year <- 1900
Final_year <- 2020
EraSize <- 30
EraLabels <- c(paste(Begin_year,"-",Begin_year+EraSize,sep=""))
for(i in seq(Begin_year+EraSize,Final_year-EraSize,EraSize)){
  EraLabels <- append(EraLabels, paste(i+1,"-",i+EraSize,sep=""))}
caterpillar_plot_delta.occ <- function(era1,era2){
  contrast.to.plot <- ifelse(era1 > era2, paste(era1,"-",era2,sep=""), paste(era2,"-",era1,sep="")) #Find the name of the contrast to plot. Make sure it is in the format [larger era]-[lesser era]
  delta.occ.to.plot <- filter(delta.occSummary, EraContrast==contrast.to.plot) %>% arrange(Mean) #Filter the delta.occupancy values to only the desired era contrast. Sort by ascending Mean
  delta.occ.to.plot$Species <- factor(delta.occ.to.plot$Species, levels=unique(delta.occ.to.plot$Species),ordered=T) #Do this to prevent ggplot from automatically sorting the species on the x axis alphabetically
  #Create caterpillar plot using ggplot2
  ggplot(delta.occ.to.plot[-which(delta.occ.to.plot$Species=="cross.sp.mean"),], aes(x = Species, y = Mean)) + #do not plot the cross.sp.mean as a vertical bar (take it out of the dataset)
    geom_errorbar(aes(ymin = Lower95, ymax = Upper95, color=ifelse(Lower95<0 & Upper95>0,"nochange",ifelse(Upper95<=0,"decline","increase"))), width = 0, position = position_dodge(width = 0), linewidth=2) +
    geom_point(position = position_dodge(width = 0), size=0.45, color="red") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    #geom_ribbon(aes(x= 1:PriorCheck$jags.data$nsp, ymin = delta.occ.to.plot[which(delta.occ.to.plot$Species=="cross.sp.mean"),]$Lower95, ymax = delta.occ.to.plot[which(delta.occ.to.plot$Species=="cross.sp.mean"),]$Upper95), fill = "red", alpha = 0.15) + #Shading for 95% CI of cross species mean
    #geom_hline(yintercept = delta.occ.to.plot[which(delta.occ.to.plot$Species=="cross.sp.mean"),]$Mean, linetype = "dashed", color = "red2", linewidth = 0.5) + #cross-species mean
    #scale_y_continuous(labels = percent_format()) +
    labs(x = paste("Species Identity (n=",PriorCheck$jags.data$nsp,")",sep=""), y = gsub("-","\u2013",paste("Change in Occupancy between ",EraLabels[ifelse(era1>era2,era2,era1)]," and ",EraLabels[ifelse(era1>era2,era1,era2)],sep=""))) +
    scale_color_manual(values=c("nochange" = "#B1B1B1", "decline" = "#AC4F55", "increase" = "#2F9BBF")) +
    guides(color=F) +
    theme_minimal() +
    theme(panel.grid.major.x=element_blank(), axis.text.x=element_blank(), axis.text.y = element_text(size = 22), plot.title = element_text(face = "bold"), axis.title.y = element_text(margin = margin(r = 10), size=14), axis.title.x = element_text(size=22))
}


p1 <- caterpillar_plot_delta.occ(3,4)
p2 <- caterpillar_plot_delta.occ(1,4)
grid.arrange(p1, p2, ncol=1)










