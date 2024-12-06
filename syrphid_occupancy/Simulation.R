# Script for running the model on simulated data created with known occupancy values, to verify its ability to accurately estimate occupancy and occupancy change.
# If you are only looking to graph existing simulation data, skip to "Compute Occupancy and Delta Occupancy" section around line 383
# Please set working directory to syrphid_occupancy folder


# Identify the simulation number (i.e., is this the first replicate simulation conducted? The second? The third?)
# This will be used to store the simulation results in a named folder
SimulationNumber <- 1
dir.create("Simulations")
dir.create(paste("Simulations/Sim",SimulationNumber,sep=""))



########################
#### Load Libraries ####
########################

library(ggplot2) #plotting tools
library(scales) #percent_format function for plotting
library(grid)
library(gridExtra)
library(tidyr) #drop_na function
library(abind) #abind function
library(dplyr) #Creating tables and other simple functions
options(dplyr.summarise.inform = FALSE) #Suppress the info messages when using dyplr::summarise()

# Libraries for Spatial Analysis:
library(sf)

# Libraries for JAGS model
library(parallel) #parallel.seeds function
library(rjags)
library(R2jags)
library(runjags)
library(coda) #For visualizing posteriors





#######################################
### Simulate Known Parameter Values ###
#######################################

sim.gridsize <- 30 #Grid size in simulation (the number of grid cells will be gridsize^2)
sim.nsp <- 50 #Number of species in simulation
sim.nsite <- sim.gridsize^2 #Number of total sites in simulation
# There are 4 eras and 6 intervals
# All parameter means will be 0 and all SD will be 0.1, with the exception of the 3 mu.psi.era parameters and 3 sigma.psi.era parameters, which will be more varied in order to create a wide range of different occupancy values for each species/era.


#### Simulated occupancy for each sp*era ####
# 1. Simulated Baseline. expit(0)=0.5, therefore the baseline occupancy probability is 50%
sim.mu.psi <- 0

# 2. Simulated psi.sp
sim.psi.sp <- numeric()
sim.sigma.psi.sp <- 0.1
for (sp in 1:sim.nsp){
  sim.psi.sp[sp] <- rnorm(1, 0, sim.sigma.psi.sp)}
#summary(sim.psi.sp) #Examine Simulated data

# 3. Simulated psi.era
sim.psi.era <- matrix(nrow=sim.nsp,ncol=4)
sim.mu.psi.era <- numeric()
sim.sigma.psi.era <- numeric()
for (era in 2:4){
  sim.mu.psi.era[era] <- runif(1, -0.5, 0.5) #Create variation in mean effect of era
  sim.sigma.psi.era[era] <- runif(1, 0.3, 0.8) #Create variation in species-specific effects of era. It doesn't seem realistic that all species will be similar in effect of era, so I would make this sd at least 0.3.
}
for (sp in 1:sim.nsp){
  sim.psi.era[sp,1] <- 0
  for (era in 2:4){
    sim.psi.era[sp,era] <- rnorm(1, sim.mu.psi.era[era], sim.sigma.psi.era[era])}}
#summary(sim.psi.era)


# Resulting simulated occupancy probability:
sim.occupancy <- matrix(nrow=sim.nsp,ncol=4)
for (sp in 1:sim.nsp){
  for (era in 1:4){
    sim.occupancy[sp,era] <- plogis(sim.mu.psi+sim.psi.sp[sp]+sim.psi.era[sp,era])
  }}
#summary(sim.occupancy)


#### Simulated detection for each sp*site*era ####
# 1. Simulated Baseline
sim.mu.p <- 0

# 2. Simulated p.sp
sim.p.sp <- numeric()
sim.sigma.p.sp <- 0.1
for (sp in 1:sim.nsp){
  sim.p.sp[sp] <- rnorm(1, 0, sim.sigma.p.sp)}
#summary(sim.p.sp)

# 3. Simulated p.site
sim.p.site <- matrix(nrow=sim.nsite, ncol=4)
sim.sigma.p.site <- 0.1
for (site in 1:sim.nsite){
  for (era in 1:4){
    sim.p.site[site,era] <- rnorm(1, 0, sim.sigma.p.site)}}
#summary(sim.p.site)

# 4. Simulated p.era
sim.p.era <- numeric()
sim.sigma.p.era <- 0.1
sim.p.era[1] <- 0
for (era in 2:4){
  sim.p.era[era] <- rnorm(1, 0, sim.sigma.p.era)}
#summary(sim.p.era)

# Resulting simulated detection probability:
sim.detection <- array(dim=c(sim.nsp,sim.nsite,4))
for (sp in 1:sim.nsp){
  for (site in 1:sim.nsite){
    for (era in 1:4){
      sim.detection[sp,site,era] <- plogis(sim.mu.p+sim.p.sp[sp]+sim.p.site[site,era]+sim.p.era[era])}}}
#View(sim.detection[,300,])


#### Simulated Visitation history for each site*era*interval ####
# Assume visits increase over time, reflecting the bias in my data towards recent time periods and improvements in methods to collect flower flies.
mu.v <- -0.5
v.era <- numeric()
v.era[1] <- 0
v.era[2] <- 0.2
v.era[3] <- 0.4
v.era[4] <- 0.6

#Save parameters
sim.params <- list(sim.mu.psi=sim.mu.psi,
                   sim.sigma.psi.sp=sim.sigma.psi.sp,
                   sim.psi.sp=sim.psi.sp,
                   sim.mu.psi.era=sim.mu.psi.era,
                   sim.sigma.psi.era=sim.sigma.psi.era,
                   sim.psi.era=sim.psi.era,
                   sim.occupancy=sim.occupancy,
                   sim.mu.p=sim.mu.p,
                   sim.sigma.p.sp=sim.sigma.p.sp,
                   sim.p.sp=sim.p.sp,
                   sim.sigma.p.site=sim.sigma.p.site,
                   sim.p.site=sim.p.site,
                   sim.p.era=sim.p.era,
                   sim.detection=sim.detection,
                   mu.v=mu.v,
                   v.era=v.era)
save(sim.params,file=paste("Simulations/Sim",SimulationNumber,"/sim.params.RData",sep=""))

#Load parameters
#load(paste("Simulations/Sim",SimulationNumber,"/sim.params.RData",sep=""))









##############################################
### Simulating Ranges and Observation Data ###
##############################################

# Create grid
sim.grid <- st_make_grid(st_bbox(c(xmin = 0, xmax = sim.gridsize, ymin = 0, ymax = sim.gridsize)), n = sim.gridsize)

# Convert the grid to a data frame
sim.grid_df <- st_sf(id = 1:length(sim.grid), geometry = sim.grid)

# Function to simulate a species' range
simulate.range <- function(seed){
  # Select 7 random cells
  #set.seed(seed) #If you want complete reproducibility
  random.cells <- sample_n(sim.grid_df, 7)
  
  # Find all gridcells within a convex hull of the random cells
  geometry.random.cells <- random.cells$geometry #Obtain the geometry of the cells where a species was found
  centroids.random.cells <- st_multipoint(st_coordinates(st_centroid(geometry.random.cells))) #Create multipoint object that contains all the centroid points of these cells
  sim.convex.hull <- st_convex_hull(centroids.random.cells) #Create a convex hull polygon around those points
  sim.range <- sim.grid_df[st_intersects(sim.grid, sim.convex.hull, sparse = FALSE),]  #Find all gridcells within the convex hull
  sim.logical.range <- ifelse(1:nrow(sim.grid_df) %in% sim.range$id, TRUE, FALSE)
  
  # In case you wish to plot the range:
  #plot(sim.grid); plot(sim.convex.hull,border=rgb(1,0,0, alpha=0.6),lwd=2,add=T); plot(geometry.random.cells, col=rgb(1,0,0, alpha=0.7), border=rgb(0,0,0, alpha=0.05), add=T); plot(sim.range$geometry,col=rgb(1,0,0, alpha=0.2),border=rgb(0,0,0, alpha=0.05), add=T)

  return(sim.logical.range)
}

# sim.sp.range is a matrix of [species x site] with logical data (using function above). Shows, for each species, which sites are within their range.
sim.sp.range <- do.call(rbind, lapply(1:sim.nsp,simulate.range))
colnames(sim.sp.range) <- paste("grid_",1:nrow(sim.grid_df),sep='')

# sim.vis.arr is an empty array of each "survey" or "visit" [site x era x interval]. Each interval in a given era is assumed to be a different "visit" to a site. In this analysis, we confine the model to only the relevant visits for each species (more on that below).
sim.vis.arr <-array(data=NA, dim=c(nrow(sim.grid_df), 4, 6), dimnames=list(paste("grid_",1:nrow(sim.grid_df),sep=''),as.vector(1:4),paste("interval_",1:6,sep='')))

# X is [species x site x era x interval]. The observation data. 0 for undetected, 1 for detected.
sim.X <- array(dim=c(sim.nsp,sim.nsite,4,6)) #Create observation array
sim.Z <- array(dim=c(sim.nsp,sim.nsite,4)) #Create Z array
sim.V <- array(dim=c(sim.nsite,4,6)) #Create Visit array
#Visit history V (whether the site was visited in this site*era*interval)
for (site in 1:sim.nsite){
  for (era in 1:4){
    for (interval in 1:6){
      
      v <- plogis(mu.v+v.era[era])
      
      sim.V[site,era,interval] <- rbinom(1,1,v)
    }}}      

for (sp in 1:sim.nsp){
  for (site in 1:sim.nsite){
    for (era in 1:4){
      
      psi <- plogis(sim.mu.psi+sim.psi.sp[sp]+sim.psi.era[sp,era])
      p <- plogis(sim.mu.p+sim.p.sp[sp]+sim.p.site[site,era]+sim.p.era[era])
      
      #Latent State Z (whether the site occupied)
      if(sim.sp.range[sp,site]==TRUE){ #Site may only be occupied if it is within species range
        Z <- rbinom(1,1,psi)
      }else{
        Z <- 0
      }
      sim.Z[sp,site,era] <- Z #Record Z State
      
      for (interval in 1:6){
        #Probability of a given sp*site*era*interval yielding a detection
        p.eff <- Z*sim.V[site,era,interval]*p #(if either Z or V is 0, detection is impossible)
        
        #Observation X
        sim.X[sp,site,era,interval] <- as.double(rbinom(1,1,p.eff))
}}}}

dimnames(sim.X) <- list(sp=1:sim.nsp, site=paste("grid_",1:nrow(sim.grid_df),sep=''), era=as.vector(1:4), visit=paste("interval_",1:6,sep=''))
dimnames(sim.Z) <- list(sp=1:sim.nsp, site=paste("grid_",1:nrow(sim.grid_df),sep=''), era=as.vector(1:4))


# Combine data objects into a single list object
sim.OccData <- list(sp.range=sim.sp.range,vis.arr=sim.vis.arr,X=sim.X, Z=sim.Z, nsp=sim.nsp, nsite=nrow(sim.grid_df), nera=4, nvisit=6)

# Keep only sites that yielded a detection of at least one species.
sim.site.keep <- which(apply(sim.OccData$X, 'site', sum)>0)
sim.OccData$X <- sim.OccData$X[,sim.site.keep,,,drop=FALSE] #Do not drop dimensions (ex. 4D object must stay 4D)
sim.OccData$sp.range <- sim.OccData$sp.range[,sim.site.keep,drop=FALSE]
sim.OccData$vis.arr <- sim.OccData$vis.arr[sim.site.keep,,,drop=FALSE]
sim.OccData$nsite <- length(sim.site.keep)

# Function to determine indexes of which visits are relevant for a given species
sim.get.indices <- function(sp.num) {
  vis.arr <- sim.OccData$vis.arr #Begin with a copy of vis.arr
  vis.arr[TRUE] <- 1 #Set all visits to 1
  nsp.detected <- apply(sim.OccData$X, 2:4, sum) #object with same dimensions as vis.arr, recording the number of species detected during each visit (visit = a certain site, era and interval... aka dimensions 2:4 of X)
  vis.arr[nsp.detected==0] <- 0 #Set a visit to 0 ("irrelevant") for this species if no species were detected during that visit/interval
  vis.arr[!sim.OccData$sp.range[sp.num,],,] <- 0 #Set a visit to 0 if the site is not in species range (i.e. sp.range==FALSE)
  
  tmp <- which(vis.arr==1, arr.ind=TRUE) #Find which cells in vis.arr are 1 (i.e., "relevant") for this species. arr.ind=TRUE means array indices are returned, because vis.arr is an array. Each row will direct you to a 1 within vis.arr[site,era,interval]. dim1 is site, dim2 is era and dim3 is interval.
  
  cbind(rep(sp.num,nrow(tmp)),tmp) #return the array indices found above, but add a column denoting the species number that these relevant visits are associated with.
}
sim.master.index <- do.call(rbind, lapply(1:sim.OccData$nsp, sim.get.indices)) #create matrix showing indexes of relevant visits for each species.
colnames(sim.master.index) <- c('sp','site','era','visit') #rename columns

# Create the dataset, containing observations for relevant visits only
sim.relevant.dat <- list(X=sim.OccData$X[sim.master.index],             #X becomes the observations (0 or 1) at relevant visits only, going through each species in series (one after the other). Sorted by species, then visit, then era, then site
                         Z=sim.OccData$Z,
                         sp.range.unaltered=sim.sp.range,
                         master.index=sim.master.index,                     #the master index, allowing the model to know the sp, site, era and interval of each observation in X
                         nsp=length(unique(sim.master.index[,'sp'])),       #number of species
                         nsite=length(unique(sim.master.index[,'site'])),   #number of relevant sites (across all spp.)
                         nera=length(unique(sim.master.index[,'era'])),     #number of eras
                         nvisit=length(unique(sim.master.index[,'visit'])), #number of intervals in an era (may be fewer than expected if some intervals are never relevant)
                         nind=nrow(sim.master.index))                       #number of relevant visits/observations total across all spp.

# Code to examine range sizes
#range.sizes <- as.numeric()
#for(sp in 1:50){
#  range.sizes[sp] <- length(which(sim.relevant.dat$sp.range.unaltered[sp,]==TRUE))}
# range.sizes


# Save dataset
save(sim.relevant.dat,file=paste("Simulations/Sim",SimulationNumber,"/sim.relevant.dat.RData",sep=""))
# Load dataset
#load(paste("Simulations/Sim",SimulationNumber,"/sim.relevant.dat.RData",sep=""))










###########################
#### Running the Model ####
###########################

# Specify the parameters to be monitored
get.params <- function() {
  c('mu.psi',         #Baseline/mean occupancy of syrphids (includes default effect of era 1)
    'mu.p',           #Mean detection probability of syrphids (includes default effect of era 1)
    'psi.sp',         #Effect of species on occupancy
    'sigma.psi.sp',   #Hyperprior sd for effects of species on occupancy. Represents variation amongst the effects of different species on occupancy.
    'p.sp',           #Effect of species on detection
    'sigma.p.sp',     #Hyperprior sd for effects of species on detection. Represents variation amongst the effects of different species on detection.
    'p.site',         #Effect of site*era on detection. Includes effect of sampling effort between sites and eras.
    'sigma.p.site',   #Hyperprior sd for effects of site*era on detection. Represents variation amongst the effects of different site*eras on detection.
    'psi.era',        #Species-specific effect of era on occupancy
    'mu.psi.era',     #Hyperprior mean for effects of era on occupancy. Represents the mean value amongst the different species' effects of era on occupancy.
    'sigma.psi.era',  #Hyperprior sd for effects of era on occupancy. Represents variation amongst the different species' effects of era on occupancy.
    'p.era') }        #Effect of era on detection.

# Function to run the model
sim.run.model <- function(sim.relevant.dat, n.iter, n.burnin, n.adapt, n.thin) { 

  model.txt <- sprintf('SyrphidModel.txt') #Bring in the BUGS-syntax model from the external source file

  # Create data object that will be passed to the model
  sim.jags.data <- list(X=sim.relevant.dat$X,                        #Observations during all relevant visits
                      era=sim.relevant.dat$master.index[,'era'],     #Eras associated with each observation
                      site=sim.relevant.dat$master.index[,'site'],   #Sites associated with each observation
                      sp=sim.relevant.dat$master.index[,'sp'],       #Species number associated with each observation
                      nsp=sim.relevant.dat$nsp,                      #Number of sp
                      nsite=sim.relevant.dat$nsite,                  #Number of sites
                      nera=sim.relevant.dat$nera,                    #Number of eras
                      nind=sim.relevant.dat$nind)                    #Number of total observations

  # Prepare initial values for MCMC
  Zst <- array(1,dim=c(sim.relevant.dat$nsp, sim.relevant.dat$nsite, sim.relevant.dat$nera)) #Creates array of [sp x site x era] filled with 1's
  make.inits <- function() {  #Function that helps initialize the random number generator of each MCMC chain when running parallel chains (on separate cores).
    RNG <- parallel.seeds("base::BaseRNG", 1)
    c(list(Z=Zst), RNG[[1]])}
  inits1 <- make.inits() #Initialize random number generator for each parallel chain you are prepared to run
  inits2 <- make.inits()
  inits3 <- make.inits()
  inits4 <- make.inits()
    
  # Then, run the model
  jags.out <- run.jags(model=model.txt,
                        monitor=get.params(),
                        data=sim.jags.data,
                        inits=list(inits1,inits2,inits3,inits4),
                        n.chains=4, #Four parallel chains being processed at a time
                        burnin=n.burnin,
                        sample=floor(n.iter/n.thin),
                        thin=n.thin,
                        adapt=n.adapt,
                        method='rjags')
  
  # Return the data object used and the model output
  list(sim.jags.data=sim.jags.data, jags.out=jags.out)
}

# Run the model! (using the run.model function defined above)
sim.model.out <- sim.run.model(sim.relevant.dat=sim.relevant.dat,
                       n.iter=30e3,   #Number of iterations of each MCMC chain
                       n.burnin=4e3, #Number of burn-in/warmup iterations
                       n.adapt=2e3,  #These iterations are used to tune the behaviour of the MCMC sampler, making it more efficient
                       n.thin=1e1)   #Thinning interval to increase independence of samples within chain. Every 10th sample is kept.

# Save Model outputs
save(sim.model.out,file=paste("Simulations/Sim",SimulationNumber,"/sim.model.out.RData",sep=""))     #SAVE
# Create Model Summary and save it
sim.ModelSummary <- summary(sim.model.out$jags.out) #Create model summary
save(sim.ModelSummary,file=paste("Simulations/Sim",SimulationNumber,"/sim.ModelSummary.RData",sep=""))  #SAVE











###############################################
#### Compute Occupancy and Delta Occupancy ####
###############################################

SimulationNumber <- 1 #Choose a simulation number to inspect/analyze
load(paste("Simulations/Sim",SimulationNumber,"/sim.model.out.RData",sep=""))
load(paste("Simulations/Sim",SimulationNumber,"/sim.ModelSummary.RData",sep=""))
load(paste("Simulations/Sim",SimulationNumber,"/sim.params.RData",sep=""))
load(paste("Simulations/Sim",SimulationNumber,"/sim.relevant.dat.RData",sep=""))


sim.MCMC <- sim.model.out$jags.out$mcmc #from here on, use MCMC as shorthand for quick/easy access to the mcmc samples in the model output



######## OCCUPANCY #########

# Function to compute estimated occupancy (as a mcmc.list object containing each sampled chain) for a given species and era
sim.compute_occupancy <- function(sp,era){
  #Occupancy is baseline occupancy + effect of species + effect of era
  mu.psi <- matrix(as.matrix(sim.MCMC[,"mu.psi"]), nrow=sim.model.out$jags.out$sample, ncol=length(sim.MCMC)) #Obtain the necessary parameter distributions
  psi.sp <- matrix(as.matrix(sim.MCMC[,paste("psi.sp[",sp,"]",sep="")]), nrow=sim.model.out$jags.out$sample, ncol=length(sim.MCMC))
  psi.era <- matrix(as.matrix(sim.MCMC[,paste("psi.era[",sp,",",era,"]",sep="")]), nrow=sim.model.out$jags.out$sample, ncol=length(sim.MCMC))
  occupancy <- as.matrix(plogis(as.matrix(mu.psi+psi.sp+psi.era))) #Compute occupancy
  occupancy <- as.mcmc.list(lapply(seq_len(ncol(occupancy)), function(i) as.mcmc(occupancy[,i]))) #convert matrix back into mcmc object (each chain i is converted into an mcmc object by as.mcmc(), which are combined into an mcmc.list object)
  return(occupancy)}

# Function to compute across-species mean occupancy in each era
sim.compute_cross.sp_occupancy <- function(era){
  mu.psi <- matrix(as.matrix(sim.MCMC[,"mu.psi"]), nrow=sim.model.out$jags.out$sample, ncol=length(sim.MCMC)) #Obtain mu.psi distribution
  if (era==1){ #If era is 1, across species occupancy is simply the expit of mu.psi. The effect of era 1 (the "default" era) is already contained in mu.psi.
    mu.psi.era <- 0
  }else{ #If era is not 1, must obtain the mu.psi.era parameter that represents the mean effect of this era across species
    mu.psi.era <- matrix(as.matrix(sim.MCMC[,paste("mu.psi.era[",era,"]",sep="")]), nrow=sim.model.out$jags.out$sample, ncol=length(sim.MCMC))
  }
  occupancy <- as.matrix(plogis(as.matrix(mu.psi+mu.psi.era))) #Compute occupancy
  occupancy <- as.mcmc.list(lapply(seq_len(ncol(occupancy)), function(i) as.mcmc(occupancy[,i]))) #convert matrix back into mcmc object
  return(occupancy)}

# Create a list object containing all occupancy posteriors (including cross-species occupancies)
sim.Occupancy_list <- list()
meanoccupancies <- list()
for (era in 1:4){ #Compute cross-species occupancies
  meanoccupancies[[paste(era)]] <- sim.compute_cross.sp_occupancy(era)}
sim.Occupancy_list[["cross.sp.mean"]] <- meanoccupancies
rm(meanoccupancies)
for (sp in 1:sim.relevant.dat$nsp) { #Compute species-specific occupancies
  spoccupancies <- list()
  for (era in 1:4) {
    spoccupancies[[paste(era)]] <- sim.compute_occupancy(sp,era)}
  sim.Occupancy_list[[paste(sp)]] <- spoccupancies #add this species' posteriors to the master list
  rm(spoccupancies)} 

# Create a dataframe of summary statistics (mean, sd, 95% intervals) for all occupancy estimates.
sim.get_occ_stats <- function(row){ #Function to get the summary stats for each era's occupancy for every species
  dataframe <- matrix(nrow=4,ncol=6) #create matrix to hold summaries
      for (era in 1:4){ #For each era, get the species name, era name, lower95, mean, upper95, and SD.
        dataframe[era,] <- c(names(sim.Occupancy_list[row]), names(sim.Occupancy_list[[row]][era]), HPDinterval(as.mcmc(unlist(sim.Occupancy_list[[row]][[era]])), prob=0.95)[[1]] , summary(sim.Occupancy_list[[row]][[era]])$statistics[[1]], HPDinterval(as.mcmc(unlist(sim.Occupancy_list[[row]][[era]])), prob=0.95)[[2]] , summary(sim.Occupancy_list[[row]][[era]])$statistics[[2]])}
  return(dataframe)}
sim.OccupancySummary <- as.data.frame(do.call(rbind,lapply(1:length(sim.Occupancy_list), sim.get_occ_stats))) #Create the dataframe
colnames(sim.OccupancySummary) <- c("Species","Era","Lower95","Mean","Upper95","SD")

# Add column showing true occupancy
#Calculate true cross-species occupancies
true.cross.sp <- numeric()
true.cross.sp[1] <- plogis(as.numeric(sim.params["sim.mu.psi"]))
for (era in 2:4){
  true.cross.sp[era] <- plogis(as.numeric(sim.params["sim.mu.psi"])+as.numeric(sim.params["sim.mu.psi.era"][[1]][era]))}
sim.occupancy.vector <- append(true.cross.sp, c(t(sim.params["sim.occupancy"][[1]])))
sim.OccupancySummary <- mutate(sim.OccupancySummary, TrueOccupancy=sim.occupancy.vector)

# Make sure columns are treated as numeric values
sim.OccupancySummary$Lower95 <- as.numeric(sim.OccupancySummary$Lower95)
sim.OccupancySummary$Mean <- as.numeric(sim.OccupancySummary$Mean)
sim.OccupancySummary$Upper95 <- as.numeric(sim.OccupancySummary$Upper95)
sim.OccupancySummary$SD <- as.numeric(sim.OccupancySummary$SD)
sim.OccupancySummary$TrueOccupancy <- as.numeric(sim.OccupancySummary$TrueOccupancy)



########## DELTA OCCUPANCY ############

# Function to compute vector of true delta occupancies between two eras
sim.compute.true.delta.occ <- function(era1,era2){
  if(era1 > era2){tmp <- era1; era1 <- era2; era2 <- tmp} #Make sure era 2 is the most recent of the two
  #Obtain cross-species delta occ
  sim.cross.sp.delta.occ <- (true.cross.sp[era2] - true.cross.sp[era1])/true.cross.sp[era1] #Percent change in occupancy
  #Obtain species-specific delta occ
  sim.occupancy <- sim.params["sim.occupancy"][[1]]
  sim.sp.delta.occ <- (sim.occupancy[,era2] - sim.occupancy[,era1])/sim.occupancy[,era1] #Percent change in occupancy
  sim.delta.occ <- append(sim.cross.sp.delta.occ, sim.sp.delta.occ)
  return(sim.delta.occ)
}

# Function to compute the estimated percent delta occupancy between two eras
sim.compute_delta.occ <- function(row,era1,era2) { 
  Occ1 <- matrix(as.matrix(sim.Occupancy_list[[row]][[era1]]), nrow=sim.model.out$jags.out$sample, ncol=length(sim.MCMC)) #Occupancy estimate of era1
  Occ2 <- matrix(as.matrix(sim.Occupancy_list[[row]][[era2]]), nrow=sim.model.out$jags.out$sample, ncol=length(sim.MCMC)) #Occupancy estimate of era2
  delta.occ <- as.matrix((Occ1-Occ2)/Occ2) #subtract the distributions and divide by the older era
  delta.occ <- as.mcmc.list(lapply(seq_len(ncol(delta.occ)), function(i) as.mcmc(delta.occ[,i]))) 
  return(delta.occ)}

# Create a list object containing all delta occupancy
sim.delta.occ_list <- list()
for (row in 1:length(sim.Occupancy_list)) {
  occ.changes <- list()
  for (era in 4:2) {
    for (nextera in (era-1):1) {
      occ.changes[[paste(era,"-",nextera, sep="")]] <- sim.compute_delta.occ(row,era,nextera)}} #Occupancy change is computed between the era in question and every era below it. Example with 4 eras: 4-3,4-2,4-1,3-2,3-1,2-1
  sim.delta.occ_list[[names(sim.Occupancy_list[row])]] <- occ.changes #add this species' contrasts to the master list
  rm(occ.changes,row,era,nextera)}

# Create a dataframe of summary statistics (mean, sd, 95% intervals) for all occupancy change estimates.
sim.get_delta_occ_stats <- function(row){ #Function to get the summary stats for all occupancy changes for a given species
  dataframe <- matrix(nrow=4*(4-1)/2,ncol=6) #create matrix to hold data. The number of pairwise contrasts is always nera*(nera-1)/2
  for (contrast in 1:(4*(4-1)/2)){ #Get the species name, contrast name, lower95, mean, upper95, and SD.
    dataframe[contrast,] <- c(names(sim.Occupancy_list[row]), names(sim.delta.occ_list[[row]][contrast]), HPDinterval(as.mcmc(unlist(sim.delta.occ_list[[row]][[contrast]])), prob=0.95)[[1]] , summary(sim.delta.occ_list[[row]][[contrast]])$statistics[[1]], HPDinterval(as.mcmc(unlist(sim.delta.occ_list[[row]][[contrast]])), prob=0.95)[[2]] , summary(sim.delta.occ_list[[row]][[contrast]])$statistics[[2]])}
  return(dataframe)}
sim.delta.occSummary <- as.data.frame(do.call(rbind,lapply(1:length(sim.Occupancy_list), sim.get_delta_occ_stats))) #Create the dataframe
colnames(sim.delta.occSummary) <- c("Species","EraContrast","Lower95","Mean","Upper95","SD")

# Make sure columns are treated as numeric values
sim.delta.occSummary$Lower95 <- as.numeric(sim.delta.occSummary$Lower95)
sim.delta.occSummary$Mean <- as.numeric(sim.delta.occSummary$Mean)
sim.delta.occSummary$Upper95 <- as.numeric(sim.delta.occSummary$Upper95)
sim.delta.occSummary$SD <- as.numeric(sim.delta.occSummary$SD)











#############################################################
#### Plotting Modelled Occupancy Against Known Occupancy ####
#############################################################

#Caterpillar plot comparing occupancy values
custom_labels <- c("cross.sp.mean-1" = "era1", "cross.sp.mean-2" = "era2", "cross.sp.mean-3" = "era3", "cross.sp.mean-4" = "era4")
caterpillar_plot_sim.occ <- function(){
  occ.to.plot <- sim.OccupancySummary %>% arrange(Mean) #Sort by ascending Mean
  occ.to.plot$Label <- paste(occ.to.plot$Species,"-",occ.to.plot$Era,sep="")
  occ.to.plot$Label <- factor(occ.to.plot$Label, levels=unique(occ.to.plot$Label),ordered=T) #Do this to prevent ggplot from automatically sorting the species on the x axis alphabetically
  #Create caterpillar plot using ggplot2
  ggplot(occ.to.plot, aes(x = Label, y = Mean)) +
    geom_errorbar(aes(ymin = Lower95, ymax = Upper95), width = 0, color="grey", position = position_dodge(width = 0), linewidth=2) +
    geom_vline(xintercept = which(occ.to.plot$Species=="cross.sp.mean"), linetype = "dashed", color = "grey") +
    geom_point(position = position_dodge(width = 0), size=1.5, color="red") +
    geom_point(data=occ.to.plot, aes(x=Label,y=TrueOccupancy), position = position_dodge(width = 0), size=1.5, color="blue") +
    #labs(title = paste("Estimated vs. True Occupancy",sep=""), subtitle= "(blue dot= true occupancy, red dot = mean estimated occupancy, bars = 95% credible interval)", x = paste("Species & era",sep=""), y = "Occupancy" ) +
    labs(x = "Species & era", y = "Occupancy") +
    scale_x_discrete(labels = function(x) ifelse(x %in% names(custom_labels), custom_labels[x], "")) +
    theme_minimal() +
    theme(panel.grid.major.x=element_blank(), plot.title = element_text(face = "bold", size=24), plot.subtitle = element_text(size=24), 
          axis.text.x = element_text(size=24), 
          axis.title = element_text(size=24), axis.text.y=element_text(size=24))
}
caterpillar_plot_sim.occ()

#What percent of true values lie within the estimated 95% CI of occupancy?
outside.95 <- 0
for (i in 1:nrow(sim.OccupancySummary)){
  if(sim.OccupancySummary$TrueOccupancy[i] < sim.OccupancySummary$Lower95[i] | sim.OccupancySummary$TrueOccupancy[i] > sim.OccupancySummary$Upper95[i]){
    outside.95 <- outside.95 + 1}}
(nrow(sim.OccupancySummary)-outside.95)/nrow(sim.OccupancySummary)

#Caterpillar plot comparing recent delta occupancy
caterpillar_plot_sim.delta.occ <- function(era1,era2){
  if(era1 > era2){tmp <- era1; era1 <- era2; era2 <- tmp}
  true.delta.occ <- as.vector(sim.compute.true.delta.occ(era1,era2))
  delta.occ.to.plot <- filter(sim.delta.occSummary,EraContrast==paste(era2,"-",era1,sep=""))
  delta.occ.to.plot <- mutate(delta.occ.to.plot, true.delta.occ=true.delta.occ) %>% arrange(Mean) #Add column with true delta occ values then sort by ascending mean
  delta.occ.to.plot$Species <- factor(delta.occ.to.plot$Species, levels=unique(delta.occ.to.plot$Species),ordered=T) #Do this to prevent ggplot from automatically sorting the species on the x axis
  assign(paste("sim.delta.occSummary.",era1,".",era2,sep=""), delta.occ.to.plot, envir=.GlobalEnv) #save this dataframe with the true delta occ values to the global environment
  #Create caterpillar plot using ggplot2
  ggplot(delta.occ.to.plot, aes(x = Species, y = Mean)) +
    geom_errorbar(aes(ymin = Lower95, ymax = Upper95), width = 0, color="grey", position = position_dodge(width = 0), linewidth=2) +
    geom_vline(xintercept = which(delta.occ.to.plot$Species=="cross.sp.mean"), linetype = "dashed", color = "grey") +
    geom_point(position = position_dodge(width = 0), size=1.5, color="red") +
    geom_point(data=delta.occ.to.plot, aes(x=Species,y=true.delta.occ), position = position_dodge(width = 0), size=1.5, color="blue") +
    scale_y_continuous(labels = percent_format(scale = 100)) +
    #ylim(-0.2,0.4) +
    #labs(title = paste("Estimated vs. True Delta Occupancy between eras ",era1," and ",era2,sep=""), subtitle= "(blue dot = true delta occupancy, red dot = mean estimated delta occupancy,\n  bars = 95% credible interval)", x = paste("Species",sep=""), y = "Delta Occupancy") + #Full title, in case plotting only one contrast
    labs(x = paste("Species Identity",sep=""), y = paste("Change in occupancy between 1961\u20131990 and 1991\u20132020 (%)")) + #Removes titles in order to plot multiple graphs in a grid with the function found below
    theme_minimal() +
    theme(panel.grid.major.x=element_blank(), plot.title = element_text(face = "bold"), axis.text.x= element_blank(), axis.title = element_text(size=12), axis.text.y = element_text(size=10) ) #If you want to bring back x axis labels, use axis.text.x = element_text(size= 3,angle = 0, vjust = 0.5, hjust=1)
}
caterpillar_plot_sim.delta.occ(3,4)

#What percent of true values lie within the estimated 95% CI of recent occ change?
outside.95 <- 0
for (i in 1:nrow(sim.delta.occSummary.3.4)){
  if(sim.delta.occSummary.3.4$true.delta.occ[i] < sim.delta.occSummary.3.4$Lower95[i] | sim.delta.occSummary.3.4$true.delta.occ[i] > sim.delta.occSummary.3.4$Upper95[i]){
    outside.95 <- outside.95 + 1}}
(nrow(sim.delta.occSummary.3.4)-outside.95)/nrow(sim.delta.occSummary.3.4)

#Caterpillar plot comparing historic delta occupancy
caterpillar_plot_sim.delta.occ <- function(era1,era2){
  if(era1 > era2){tmp <- era1; era1 <- era2; era2 <- tmp}
  true.delta.occ <- as.vector(sim.compute.true.delta.occ(era1,era2))
  delta.occ.to.plot <- filter(sim.delta.occSummary,EraContrast==paste(era2,"-",era1,sep=""))
  delta.occ.to.plot <- mutate(delta.occ.to.plot, true.delta.occ=true.delta.occ) %>% arrange(Mean) #Add column with true delta occ values then sort by ascending mean
  delta.occ.to.plot$Species <- factor(delta.occ.to.plot$Species, levels=unique(delta.occ.to.plot$Species),ordered=T) #Do this to prevent ggplot from automatically sorting the species on the x axis
  assign(paste("sim.delta.occSummary.",era1,".",era2,sep=""), delta.occ.to.plot, envir=.GlobalEnv) #save this dataframe with the true delta occ values to the global environment
  #Create caterpillar plot using ggplot2
  ggplot(delta.occ.to.plot, aes(x = Species, y = Mean)) +
    geom_errorbar(aes(ymin = Lower95, ymax = Upper95), width = 0, color="grey", position = position_dodge(width = 0), linewidth=2) +
    geom_vline(xintercept = which(delta.occ.to.plot$Species=="cross.sp.mean"), linetype = "dashed", color = "grey") +
    geom_point(position = position_dodge(width = 0), size=1.5, color="red") +
    geom_point(data=delta.occ.to.plot, aes(x=Species,y=true.delta.occ), position = position_dodge(width = 0), size=1.5, color="blue") +
    scale_y_continuous(labels = percent_format(scale = 100)) +
    #ylim(-0.2,0.4) +
    #labs(title = paste("Estimated vs. True Delta Occupancy between eras ",era1," and ",era2,sep=""), subtitle= "(blue dot = true delta occupancy, red dot = mean estimated delta occupancy,\n  bars = 95% credible interval)", x = paste("Species",sep=""), y = "Delta Occupancy") + #Full title, in case plotting only one contrast
    labs(x = "Species Identity", y = paste("Change in occupancy between 1900\u20131930 and 1991\u20132020 (%)")) + #Removes titles in order to plot multiple graphs in a grid with the function found below
    theme_minimal() +
    theme(panel.grid.major.x=element_blank(), plot.title = element_text(face = "bold"), axis.text.x= element_blank(), axis.title = element_text(size=12), axis.text.y = element_text(size=10) ) #If you want to bring back x axis labels, use axis.text.x = element_text(size= 3,angle = 0, vjust = 0.5, hjust=1)
}
caterpillar_plot_sim.delta.occ(1,4)

#What percent of true values lie within the estimated 95% CI of historic occ change?
outside.95 <- 0
for (i in 1:nrow(sim.delta.occSummary.1.4)){
  if(sim.delta.occSummary.1.4$true.delta.occ[i] < sim.delta.occSummary.1.4$Lower95[i] | sim.delta.occSummary.1.4$true.delta.occ[i] > sim.delta.occSummary.1.4$Upper95[i]){
    outside.95 <- outside.95 + 1}}
(nrow(sim.delta.occSummary.1.4)-outside.95)/nrow(sim.delta.occSummary.1.4)


#Plot all delta.occ
library(cowplot)
library(scales)
plot.all_caterpillar_plot_sim.delta.occ <- function(){
  title <- ggplot() + 
    labs(title = paste("Estimated vs. True Delta Occupancy between all Eras",sep=""), subtitle= "(blue dot = true delta occupancy, red dot = mean estimated delta occupancy,\n  bars = 95% credible interval)") +
    theme(plot.title = element_text(face = "bold"))
  a <- caterpillar_plot_sim.delta.occ(4,3)
  b <- caterpillar_plot_sim.delta.occ(4,2)
  c <- caterpillar_plot_sim.delta.occ(4,1)
  d <- caterpillar_plot_sim.delta.occ(3,2)
  e <- caterpillar_plot_sim.delta.occ(3,1)
  f <- caterpillar_plot_sim.delta.occ(2,1)
  plots <- plot_grid(a,b,c,d,e,f, nrow=2,labels=c("Eras 4-3","Eras 4-2","Eras 4-1","Eras 3-2","Eras 3-1","Eras 2-1"), label_x = 0, label_y = 1, label_size = 12)
  plot_grid(title,plots,ncol=1, rel_heights=c(1,12))
}
plot.all_caterpillar_plot_sim.delta.occ()










############################################################
# Appendix 4F: Simulation results from simulations 2 and 3 #
############################################################

caterpillar_plot_sim.occ.appendix <- function(){
  occ.to.plot <- sim.OccupancySummary %>% arrange(Mean) #Sort by ascending Mean
  occ.to.plot$Label <- paste(occ.to.plot$Species,"-",occ.to.plot$Era,sep="")
  occ.to.plot$Label <- factor(occ.to.plot$Label, levels=unique(occ.to.plot$Label),ordered=T) #Do this to prevent ggplot from automatically sorting the species on the x axis alphabetically
  #Create caterpillar plot using ggplot2
  plot <- ggplot(occ.to.plot, aes(x = Label, y = Mean)) +
    geom_errorbar(aes(ymin = Lower95, ymax = Upper95), width = 0, color="grey", position = position_dodge(width = 0), linewidth=0.8) +
    geom_point(position = position_dodge(width = 0), size=0.6, color="red") +
    geom_point(data=occ.to.plot, aes(x=Label,y=TrueOccupancy), position = position_dodge(width = 0), size=0.6, color="blue") +
    #labs(title = paste("Estimated vs. True Occupancy",sep=""), subtitle= "(blue dot= true occupancy, red dot = mean estimated occupancy, bars = 95% credible interval)", x = paste("Species & era",sep=""), y = "Occupancy" ) +
    scale_y_continuous(limits=c(0.1,1)) +
    labs(x = ifelse(SimulationNumber==2,"Species & Era","")) +
    labs(y = ifelse(SimulationNumber==1, "Occupancy","")) +
    theme_minimal() +
    if(SimulationNumber==1){
      theme(panel.grid.major.x=element_blank(), axis.text.x = element_blank(), axis.title = element_text(size=12), axis.text.y = element_text(size=10))}
    else{
      theme(panel.grid.major.x=element_blank(), axis.text.x = element_blank(), axis.title = element_text(size=12), axis.text.y = element_blank())}
  return(plot)
  }

caterpillar_plot_sim.delta.occ.appendix <- function(era1,era2){
  if(era1 > era2){tmp <- era1; era1 <- era2; era2 <- tmp}
  true.delta.occ <- as.vector(sim.compute.true.delta.occ(era1,era2))
  delta.occ.to.plot <- filter(sim.delta.occSummary,EraContrast==paste(era2,"-",era1,sep=""))
  delta.occ.to.plot <- mutate(delta.occ.to.plot, true.delta.occ=true.delta.occ) %>% arrange(Mean) #Add column with true delta occ values then sort by ascending mean
  delta.occ.to.plot$Species <- factor(delta.occ.to.plot$Species, levels=unique(delta.occ.to.plot$Species),ordered=T) #Do this to prevent ggplot from automatically sorting the species on the x axis
  assign(paste("sim.delta.occSummary.",era1,".",era2,sep=""), delta.occ.to.plot, envir=.GlobalEnv) #save this dataframe with the true delta occ values to the global environment
  #Create caterpillar plot using ggplot2
  if(era2-era1==1){
    limit <- c(-1,3)
  }else{
    limit <- c(-0.7,1)
  }
  plot <- ggplot(delta.occ.to.plot, aes(x = Species, y = Mean)) +
    geom_errorbar(aes(ymin = Lower95, ymax = Upper95), width = 0, color="grey", position = position_dodge(width = 0), linewidth=3) +
    geom_point(position = position_dodge(width = 0), size=1.5, color="red") +
    geom_point(data=delta.occ.to.plot, aes(x=Species,y=true.delta.occ), position = position_dodge(width = 0), size=1.5, color="blue") +
    scale_y_continuous(limits= limit, labels = percent_format()) +
    #labs(title = paste("Estimated vs. True Delta Occupancy between eras ",era1," and ",era2,sep=""), subtitle= "(blue dot = true delta occupancy, red dot = mean estimated delta occupancy,\n  bars = 95% credible interval)", x = paste("Species",sep=""), y = "Delta Occupancy") + #Full title, in case plotting only one contrast
    labs(x = ifelse(SimulationNumber==2,"Species Identity","")) +
    labs(y = ifelse(SimulationNumber==1, ifelse(era2-era1==1, "Recent Change in Occupancy (%)","Historic Change in Occupancy (%)"),"" ) ) + #Removes titles in order to plot multiple graphs in a grid with the function found below
    theme_minimal() +
    if(SimulationNumber==1){
      theme(panel.grid.major.x=element_blank(), plot.title = element_text(face = "bold"), axis.text.x= element_blank(), axis.title = element_text(size=12), axis.text.y = element_text(size=10))}
    else{
      theme(panel.grid.major.x=element_blank(), plot.title = element_text(face = "bold"), axis.text.x= element_blank(), axis.title = element_text(size=12), axis.text.y = element_blank())}
  return(plot)}


occ1 <- caterpillar_plot_sim.occ.appendix()
recent1 <- caterpillar_plot_sim.delta.occ.appendix(3,4)
historic1 <- caterpillar_plot_sim.delta.occ.appendix(1,4)

occ2 <- caterpillar_plot_sim.occ.appendix()
recent2 <- caterpillar_plot_sim.delta.occ.appendix(3,4)
historic2 <- caterpillar_plot_sim.delta.occ.appendix(1,4)

occ3 <- caterpillar_plot_sim.occ.appendix()
recent3 <- caterpillar_plot_sim.delta.occ.appendix(3,4)
historic3 <- caterpillar_plot_sim.delta.occ.appendix(1,4)

spacer <- nullGrob()
grid.arrange(occ1,occ2,occ3,spacer,spacer,spacer,recent1,recent2,recent3,spacer,spacer,spacer,historic1,historic2,historic3,ncol=3,heights=c(1,0.1,1,0.1,1),widths=c(1,0.9,0.9))


















