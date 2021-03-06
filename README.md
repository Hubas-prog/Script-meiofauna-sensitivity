R scripts of the statistical analysis for the article entitled:

# First assessment of the benthic meiofauna sensitivity to low human-impacted mangroves in French Guiana

Citation : Michelet C, Zeppilli D, Hubas C, Baldrighi E, Cuny P, Dirberg G, Militon C, Walcker R, Lamy D, Jézéquel R, Receveur J, Gilbert F, Houssainy AE, Dufour A, Heimbürger-Boavida L-E, Bihannic I, Sylvi L, Vivier B, Michaud E. First Assessment of the Benthic Meiofauna Sensitivity to Low Human-Impacted Mangroves in French Guiana. Forests. 2021; 12(3):338.

Published in Forest (MDPI) : https://doi.org/10.3390/f12030338

The current deposit corresponds only to the statistical analysis. It gives details about the statistical procedure that was developped for this article and called Between Class Multiple Factor Analysis (BC-MFA). The raw data are available at the following address: https://doi.org/10.5281/zenodo.4592299

The study provides measurments of various parameters on mangrove sediment core at 3 sampling sites in French Guiana. Statistical Individuals are coded according to the sanpling station S1,S2,S3 and sediment core depth layers L1,L2,L3 (increasing number with sediment depth).

'SCRIPT ENVIRONNEMENT.R' file contains the in-house R script used for data processing, univariate and multivariate statistics as well as figures corresponding to the analysis of environmental data. The script call the dataset hosted in https://doi.org/10.5281/zenodo.4592299. Contaminants.csv (sep=";",dec=",") contains trace metals and organic pollutants. pigments.csv contains lipophilic pigments measured by HPLC (sep=";",dec=","). CHN.csv contains Carbone and Nitrogen data measured by Element Analyser (sep=";",dec=","), Communaute-Bacterienne.csv (sep=";",dec=",") contains data about DNA quantity, bacteria and archea biomass. dataENV.csv (sep=";",dec=",") contains environmental data such as redox potential or grain size analysis

'SCRIPT NEMATODES.R' contains the in-house R script used for data processing, univariate and multivariate statistics as well as figures corresponding to the analysis of nematode communities. The script call the dataset hosted in https://doi.org/10.5281/zenodo.4592299. Nematodes_pourcent_bacterivores.csv (sep=";",dec=","), Nematodes_pourcent_detritivores.csv (sep=";",dec=","), Nematodes_pourcent_brouteurs.csv (sep=";",dec=",")and Nematodes_pourcent_predateurs.csv (sep=";",dec=",") contains relative abundances (in %) of bacterivorous, detritivorous, grazers and Omnivorous/predators nematodes.

'SCRIPT NEMATODES BIOMASS - Supp_Mat.R' contains the in-house R script used for data processing, univariate and multivariate statistics as well as figures corresponding to the analysis of nematode communities. The script call the dataset hosted in https://doi.org/10.5281/zenodo.4592299. Nematodes_biomasscm3_pourcent_bacterivores.csv (sep=";",dec=","), Nematodes_biomasscm3_pourcent_detritivores.csv (sep=";",dec=","), Nematodes_biomasscm3_pourcent_brouteurs.csv (sep=";",dec=",")and Nematodes_biomasscm3_pourcent_predateurs.csv (sep=";",dec=",") contains relative biomass (in %) of bacterivorous, detritivorous, grazers and Omnivorous/predators nematodes. They contain relative biomasses of the above mentioned free living nematodes trophic groups. 

List of R packages :
library(ade4) => for multivariate statistics and creattion of duality diagrams.
library(vegan)=> for improved data visualization of ade4 outputs.
library(scales)=> for controlling opacity of plot features.
