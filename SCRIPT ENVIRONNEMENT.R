########################### SCRIPT ################################
# First assessment of the benthic meiofauna sensitivity to low human-impacted mangroves in French Guiana
# by : Claire Michelet, Daniela Zeppilli, Cédric Hubas, Elisa Baldrighi, Philippe Cuny, Guillaume Dirberg, Cécile Militon, Romain Walcker, Dominique Lamy, Ronan Jézéquel, Justine Receveur, Franck Gilbert, Amonda El Houssainy, Aurélie Dufour, Lars-Eric Heimbürger-Boavida, Isabelle Bihannic, Léa Sylvi, Baptiste Vivier, Emma Michaud
# For publication in Forest (MDPI)
# *Corresponding author : emma.michaud@univ-brest.fr
########################### SCRIPT ################################

#####################
# PACKAGES
#####################
library(ade4)
library(vegan)
library(scales)

#####################
# AESTHETICS
#####################
my.palette<-colorRampPalette(c("purple","black","red1","green3","blue"))

#####################
# DATA upload
#####################

# Contaminants
conta<-read.csv("Contaminants.csv",h=T,sep=";",dec=",")
# Pigments
pigments<-read.csv("pigments.csv",h=T,sep=";",dec=",")
# C/N ratio
CHN<-read.csv("CHN.csv",h=T,sep=";",dec=",")
# Bacterial community
bact<-read.csv("Communaute-Bacterienne.csv",h=T,sep=";",dec=",")
# Environemental data
enviro<-read.csv("dataENV.csv",h=T,sep=";",dec=",")

# verification that all statistical individuals are paired
data.frame(conta[,1],pigments[,1],enviro[,1],bact[,1],CHN[,1])

# final dataset
DATA<-cbind(conta[,2:32],pigments[,2:19],bact[,3:4],enviro[,3:7],CHN[,2:3])
DATA$sites<-substr(conta$groupe,1,2)
DATA$layer<-paste("L",substr(conta$groupe,4,4),sep="")
DATA$core<-substr(conta$groupe,3,3)
names(DATA)

bloc<-c(dim(conta[,2:32])[2],
        dim(pigments[,2:19])[2],
        dim(bact[,3:4])[2],
        dim(enviro[,3:7])[2],
        dim(CHN[,2:3])[2])

#####################
# Multiple FActor Analysis (unsupervised)
#####################

# Extraction of eigen values
eig.conta<-dudi.pca(conta[,2:32],scannf=F,nf=2)$eig 
eig.pig<-dudi.pca(pigments[,2:19],scannf=F,nf=2)$eig 
eig.bact<-dudi.pca(bact[,3:4],scannf=F,nf=2)$eig 
eig.enviro<-dudi.pca(enviro[,3:7],scannf=F,nf=2)$eig 
eig.CHN<-dudi.pca(CHN[,2:3],scannf=F,nf=2)$eig

# MFA
MFA<-dudi.pca(DATA[,1:58],col.w=rep(c(1/eig.conta[1],
                               1/eig.pig[1],
                               1/eig.bact[1],
                               1/eig.enviro[1],
                               1/eig.CHN[1]
                               ),
                             bloc),
              scannf=F,
              nf=2)

varexp1<-MFA$eig*100/sum(MFA$eig)
fact<-factor(DATA$sites)
fact2<-factor(paste(DATA$sites,DATA$layer))

#####################
# Between class MFA (supervised MFA analysis with factor = sites)
#####################

BCA<-bca(MFA,fact,scannf=F,nf=2)
BCA$ratio
varexp2<-BCA$eig*100/sum(BCA$eig)

randtest(BCA)

#####################
# Plots
#####################

par(layout(matrix(c(1,2,3,3),nrow=2)),mar=c(4,4,4,4))

plot(MFA$li[,2]~MFA$li[,1],
     pch=21,cex=2,col="white",
     bg=rainbow(length(levels(fact)))[fact],
     xlab=paste("Axis 1 : ",round(varexp1[1],2),"%"),
     ylab=paste("Axis 2 : ",round(varexp1[2],2),"%"),
     main="MFA scores")
ordihull(MFA$li,fact2,lab=T)

plot(BCA$ls[,2]~BCA$ls[,1],
     pch=21,cex=2,col="white",
     bg=rainbow(length(levels(fact)))[fact],
     xlab=paste("Axis 1 : ",round(varexp2[1],2),"%"),
     ylab=paste("Axis 2 : ",round(varexp2[2],2),"%"),
     main="BCA scores")

legend("topright",c("instrumental variable = Sites",
					"Total inertia explained = 39.8%",
					"p-value=0.001"),cex=0.8,box.lty=0)

ordihull(BCA$ls,fact,lab=T)

plot(BCA$co[,2]~BCA$co[,1],type="n",
     xlab=paste("Axis 1 : ",round(varexp2[1],2),"%"),
     ylab=paste("Axis 2 : ",round(varexp2[2],2),"%"),
     ylim=c(-1,1),xlim=c(-1,1),main="BCA loadings")

arrows(x0=0,y0=0,x1=BCA$co[,1],y1=BCA$co[,2],col="lightgrey",length=0.1)

cos2 <- as.matrix(BCA$co[,1:2])*as.matrix(BCA$co[,1:2]) 

var.group<-factor(rep(c("conta","pig","bact","enviro","CHN"),bloc),
                  levels=c("enviro","conta","bact","CHN","pig"))
BCA$co$col<-my.palette(length(bloc))[var.group]

newco <- BCA$co[cos2[,1]>0.4 | cos2[,2]>0.4,]
oldco <- BCA$co[cos2[,1]<0.4 & cos2[,2]<0.4,]

text(newco[,1], newco[,2],rownames(newco),col=newco$col,cex=1.2)
text(oldco[,1], oldco[,2],rownames(oldco),col=alpha(oldco$col,0.3),cex=0.8)
abline(v=0,lty="dashed")
abline(h=0,lty="dashed")

legend("topleft",
       c("Physicochemical parameters",
         "Contaminants",
         "Prokaryotic biomass",
         "Organic matter",
         "Pigments"),
       text.col=my.palette(length(bloc)),cex=0.8,box.lty=0)
