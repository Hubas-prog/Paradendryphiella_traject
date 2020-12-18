# Paradendriphyella

Supplementary materials for paper entitled: Laminariales Host Does Impact Lipid Temperature Trajectories of the Fungal Endophyte Paradendryphiella salina (Sutherland.) by Vallet_et_al. published in MDPI, Marine Drugs.

MDPI and ACS Style reference: Vallet, M.; Meziane, T.; Thiney, N.; Prado, S.; Hubas, C. Laminariales Host Does Impact Lipid Temperature Trajectories of the Fungal Endophyte Paradendryphiella salina (Sutherland.). Mar. Drugs 2020, 18, 379.

DOI: https://doi.org/10.3390/md18080379

txt files contain raw GC-FID data used in the analysis, tailored for ease of manipulation with R. 'R_Script_Vallet_et_al_mar_drugs.R' file contains the in-house R script used for data processing, univariate and multivariate statistics as well as figures. 'fill.C23copy.txt' file contains the amount of internal standard (C23) and the weight (in mg) of each sample. 

List of packages :

library(reshape) => to perform a pivot table from the list of txt files.
library(ggplot2)=> for data visualization.
library(rstatix)=> for simple and intuitive pipe-friendly framework (eeded to perform Welch ANOVA).
library(ade4) => for multivariate statistics.
library(factoextra)=> for improved data visualization of ade4 outputs.
library(cowplot)=> for publication-quality figures with 'ggplot2'.

Script correspondance:

'R_Script_Vallet_et_al_mar_drugs.R' generate stats and figures. It requires the txt files to generate a list of table and the fatty acid semiquantitative table (in %). It also requires, in addition, 'fill.C23copy.txt' to generate the concentration table.
