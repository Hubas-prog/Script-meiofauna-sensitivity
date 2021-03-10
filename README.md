# Script-meiofauna-sensitivity

R scripts of the statistical analysis for the article entitled:

First assessment of the benthic meiofauna sensitivity to low human-impacted mangroves in French Guiana

by  Claire Michelet, Daniela Zeppilli, Cédric Hubas, Elisa Baldrighi, Philippe Cuny, Guillaume Dirberg, Cécile Militon, Romain Walcker, Dominique Lamy, Ronan Jézéquel, Justine Receveur, Franck Gilbert, Amonda El Houssainy, Aurélie Dufour, Lars-Eric Heimbürger-Boavida, Isabelle Bihannic, Léa Sylvi, Baptiste Vivier, Emma Michaud

To be published in Forest (MDPI)

The current deposit corresponds only to the statistical analysis. The raw data are available at the following address: https://doi.org/10.5281/zenodo.4592299

The study provides measurments of various parameters on mangrove sediment core at 3 sampling sites in French Guiana. Statistical Individuals are coded according to the sanpling station S1,S2,S3 and sediment core depth layers L1,L2,L3 (increasing number with sediment depth).

'SCRIPT ENVIRONNEMENT.R' file contains the in-house R script used for data processing, univariate and multivariate statistics as well as figures corresponding to the analysis of environmental data. The script call the dataset hosted by the above mentioned address. Contaminants.csv (sep=";",dec=",") contains trace metals and organic pollutants. pigments.csv contains lipophilic pigments measured by HPLC (sep=";",dec=","). CHN.csv contains Carbone and Nitrogen data measured by Element Analyser (sep=";",dec=","), Communaute-Bacterienne.csv (sep=";",dec=",") contains data about DNA quantity, bacteria and archea biomass. dataENV.csv (sep=";",dec=",") contains environmental data such as redox potential or grain size analysis

'SCRIPT NEMATODES.R' contains the in-house R script used for data processing, univariate and multivariate statistics as well as figures corresponding to the analysis of nematode communities.

List of packages :
library(ade4) => for multivariate statistics.
library(factoextra)=> for improved data visualization of ade4 outputs.
library(cowplot)=> for publication-quality figures with 'ggplot2'.
