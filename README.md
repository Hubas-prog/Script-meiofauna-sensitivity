# Script-meiofauna-sensitivity

R scripts of the statistical analysis for the article entitled:

First assessment of the benthic meiofauna sensitivity to low human-impacted mangroves in French Guiana

by  Claire Michelet, Daniela Zeppilli, Cédric Hubas, Elisa Baldrighi, Philippe Cuny, Guillaume Dirberg, Cécile Militon, Romain Walcker, Dominique Lamy, Ronan Jézéquel, Justine Receveur, Franck Gilbert, Amonda El Houssainy, Aurélie Dufour, Lars-Eric Heimbürger-Boavida, Isabelle Bihannic, Léa Sylvi, Baptiste Vivier, Emma Michaud

To be published in Forest (MDPI)

The current deposit corresponds only to the statistical analysis. The raw data are available at the following address: https://doi.org/10.5281/zenodo.4592299

'SCRIPT ENVIRONNEMENT.R' file contains the in-house R script used for data processing, univariate and multivariate statistics as well as figures corresponding to the analysis of environmental data. The script call the dataset hosted by the above mentioned address. # Contaminants
# Contaminants are called by the command 
conta<-read.csv("Contaminants.csv",h=T,sep=";",dec=",")
# Pigments
pigments<-read.csv("pigments.csv",h=T,sep=";",dec=",")
# C/N ratio
CHN<-read.csv("CHN.csv",h=T,sep=";",dec=",")
# Bacterial community
bact<-read.csv("Communaute-Bacterienne.csv",h=T,sep=";",dec=",")
# Environemental data
enviro<-read.csv("dataENV.csv",h=T,sep=";",dec=",")



'SCRIPT NEMATODES.R' contains the in-house R script used for data processing, univariate and multivariate statistics as well as figures corresponding to the analysis of nematode communities.

List of packages :

library(ade4) => for multivariate statistics.
library(factoextra)=> for improved data visualization of ade4 outputs.
library(cowplot)=> for publication-quality figures with 'ggplot2'.

Script correspondance:

'SCRIPT ENVIRONNEMENT.R' generate stats and figures. It requires the .csv files available at the above mentioned address.
