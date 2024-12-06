# Assessing Historical Occupancy Trends in North American Flower Flies (Diptera: Syrphidae)

**Overview**

Given the declines occurring in many insect populations around the world, there is an urgent need to assess large-scale trends in important groups like pollinators. Syrphidae is one of the largest families of flies and one of the most common pollinators of major crops after bees. Museum data and biological collections present unique opportunities to assess population trends over long historical timescales. Here, I use a hierarchical Bayesian multi-season, multi-species occupancy model on digital syrphid records from seven collections across North America to estimate both species-specific and overall syrphid occupancy change between the current time period (1991–2020) and two baseline time periods: recent (1961–1990) and historical (1900–1930). I also validate the accuracy and sensitivity of the model by using Prior Predictive Checks, Simulations, Posterior Predictive Checks, and a Sensitivity Analysis. Finally, I compare species traits (range size, body size, subfamily, larval feeding guild, and rot-hole specialization) between the declining and increasing species, to assess potential factors contributing to occupancy change.
 
 
**Description of Contents**

*Datasets* —— This folder (which must be unzipped) contains digital syrphid specimen records gathered from 7 collections across North America. It consists entirely of publicly accessible records.

*occupancy_syrphidae_main.R* —— This is the main script. It prepares the data, runs the model, outputs results, creates figures, etc.

*SyrphidModel.txt* —— This text file contains the occupancy model, written in BUGS language to be used by JAGS. To run the model, the R scripts use Rjags to interpret this file.

*FinalModel* —— This folder contains the outputs of the main script, including species tables and occupancy estimates.

*PriorPredictiveCheck.R* —— This script handles the analysis and plotting of the Prior Predictive Check outputs, which must first be obtained from syrphidae_occupancy_main.R.

*PriorPredCheck* —— This folder contains the outputs of the Prior Predictive Check, produced by syrphidae_occupancy_main.R to be assessed/plotted in PriorPredictiveCheck.R.

*Simulation.R* —— This script simulates data with known occupancy values, then runs the model on that data to assess the accuracy/reliability of the model.

*Simulations* —— This folder contains the outputs of the simulation tests produced by Simulation.R.

*PosteriorPredictiveChecks.RData* —— This R object is produced by syrphidae_occupancy_main.R to hold the outputs of the Posterior Predictive Checks, which simulate data using model estimates to evaluate how well the model estimates can explain the real observations.

*SensitivityModels.R* —— This script is used to run the occupancy model under different grid cell and time period sizes, to assess the sensitivity of my conclusions to the chosen spatial/temporal scale.

*SensitivityAnalysis* —— This folder contains the model outputs for the model variants run in SensitivityModels.R.

<br />&emsp;

NOTE: Due to file size constraints, the raw model outputs ("model.out.RData") could not be uploaded. However, this is not much of an issue because the raw outputs were only necessary for computing occupancy and occupancy change, which are summarized in the uploaded objects OccupancySummary.RData and delta.occSummary.RData, respectively.
