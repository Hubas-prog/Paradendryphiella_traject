# Paradendriphyella

Supplementary materials for paper entitled: Identity and sequence: The effect of multiple stressors on microphytobenthos assemblages by James E V Rimmer, Cédric Hubas, Adam Wyness, Bruno Jesus, Morgan Hartley, Andrew J Blight, Antoine Prins, David M Paterson

xlsx files contain the data used in the analysis, tailored for ease of manipulation with R. 'Previous results to guide priors.xlsx' contains data from previous experiments which were used to guide the priors used in the Bayesian analysis. Rda files contain R processed objects (e.g. dataframes); their production can be recreated using code in the scripts. R files are the scripts used for the key analysis.

Many functions used in the modelling are unique to the rstanarm package (https://mc-stan.org/rstanarm/). Stan priors can be specified manually or be allowed to scale to the data automatically (autoscale = TRUE) - the scaled priors used in this analysis are specified in the code to future-proof against changes to the package. See ?stan_lmer or ?stan_glm for more details.

Script correspondance:

'Diatom community analysis.R' requires the xlsx file 'Diatom assemblage.xlsx'.

'Data object creation.R' and 'Control models.R' requires the xlsx files 'Chlorophyll.xlsx', 'Critical erosion thresholds.xlsx', 'PAM time 1.xlsx', 'PAM time 2.xlsx', 'PAM time 3.xlsx'.

'Disturbance and sequence models.R' relies on the R objects created by 'Data object creation.R', or the provided .rda files.
