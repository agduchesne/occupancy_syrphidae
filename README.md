# Assessing Historical Occupancy Trends in North American Flower Flies (Diptera: Syrphidae)

**Overview**

This code is an accompaniment to the paper:
<br>
Duchesne, A.G., Des Marteaux, L., Skevington, J., Dawson, J., & Martin, A.E. Evidence of recent declines in flower flies (Diptera: Syrphidae) across North America.

Given the declines occurring in many insect populations around the world, there is an urgent need to assess large-scale trends in important groups like pollinators. Syrphidae is one of the largest families of flies and one of the most common pollinators of major crops after bees. Museum data and biological collections present unique opportunities to assess population trends over long historical timescales. Here, I use a hierarchical Bayesian multi-season, multi-species occupancy model on digital syrphid records from seven collections across North America to estimate both species-specific and overall syrphid occupancy change between the current time period (1991–2020) and two baseline time periods: recent (1961–1990) and historical (1900–1930). I also validate the accuracy and sensitivity of the model by using Prior Predictive Checks, Simulations, Posterior Predictive Checks, and a Sensitivity Analysis.
 
**Description of Contents**

*Datasets* —— This folder contains digital syrphid specimen records gathered from 7 collections across North America. It consists entirely of publicly accessible records.

*occupancy_syrphidae_main.R* —— This is the main R script. It prepares the data, runs the model, outputs results, creates figures, etc.

*SyrphidModel.txt* —— This text file contains the occupancy model, written in BUGS language to be used by JAGS. To run the model, the R scripts use the package 'Rjags' to interpret this file.

*FinalModel* —— This folder contains the outputs of the main script, including species tables and occupancy estimates.

*PriorPredictiveCheck.R* —— This R script analyzes and plots the Prior Predictive Check output, which is created by syrphidae_occupancy_main.R and stored in the PriorPredCheck folder.

*PriorPredCheck* —— This folder contains the outputs of the Prior Predictive Check, produced by syrphidae_occupancy_main.R to be assessed/plotted in PriorPredictiveCheck.R.

*Simulation.R* —— This R script simulates data using known occupancy values, then runs the model on that data to assess the accuracy/reliability of the model.

*Simulations* —— This folder contains the outputs of the Simulation tests produced by Simulation.R.

*PosteriorPredCheck* —— This folder contains the outputs of the Posterior Predictive Checks performed by syrphidae_occupancy_main.R. The Posterior Predictive Check simulates data using model estimates to evaluate how well the model estimates can explain the real observations.

*SensitivityModels.R* —— This R script was used to run the occupancy model under different grid cell and time period sizes, to assess the sensitivity of my conclusions to the chosen spatial/temporal scales.

*SensitivityAnalysis* —— This folder contains the model outputs for each model variant of the Sensitivity Analysis created in SensitivityModels.R.

<br />&emsp;

NOTE: Due to file size constraints, the raw model outputs (referenced in the R scripts as "model.out.RData") could not be uploaded. These raw outputs contained the 6000 posterior samples generated for each parameter in the model. The absence of this file should not be an issue because the raw outputs were only necessary for computing occupancy and occupancy change, which are summarized in the uploaded objects OccupancySummary.RData and delta.occSummary.RData, respectively.
