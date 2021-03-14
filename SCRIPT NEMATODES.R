########################### SCRIPT ################################
# Supplementary Material of the following article :
# First assessment of the benthic meiofauna sensitivity to low human-impacted mangroves in French Guiana
# Claire Michelet, Daniela Zeppilli, Cédric Hubas, Elisa Baldrighi, Philippe Cuny, Guillaume Dirberg, Cécile Militon, Romain Walcker, Dominique Lamy, Ronan Jézéquel, Justine Receveur, Franck Gilbert, Amonda El Houssainy, Aurélie Dufour, Lars-Eric Heimbürger-Boavida, Isabelle Bihannic, Léa Sylvi, Baptiste Vivier, Emma Michaud
# Forests. 2021; 12(3):338.
# https://doi.org/10.3390/f12030338
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

# Bacterivorous nematodes - density (in %)
nemad.pourcent.bact<-read.csv("Nematodes_pourcent_bacterivores.csv",
                              h=T,sep=";",dec=",")
# Detritivorous nematodes - density (in %)
nemad.pourcent.detri<-read.csv("Nematodes_pourcent_detritivores.csv",
                               h=T,sep=";",dec=",")
# Grazers nematodes - density (in %)
nemad.pourcent.brou<-read.csv("Nematodes_pourcent_brouteurs.csv",
                              h=T,sep=";",dec=",")
# Omnivorous/Predators nematodes - density (in %)
nemad.pourcent.pred<-read.csv("Nematodes_pourcent_predateurs.csv",
                              h=T,sep=";",dec=",")

# verification that all statistical individuals are paired
names(nemad.pourcent.bact)
names(nemad.pourcent.detri)
names(nemad.pourcent.brou)
names(nemad.pourcent.pred)

data.frame(nemad.pourcent.bact[,1],
           nemad.pourcent.detri[,1],
           nemad.pourcent.brou[,1],
           nemad.pourcent.pred[,1])

# final datasett
DATA<-cbind(nemad.pourcent.bact[,2:17],
            nemad.pourcent.detri[,2:27],
            nemad.pourcent.brou[,2:27],
            nemad.pourcent.pred[,2:8])

DATA$sites<-substr(nemad.pourcent.bact$groupe,1,2)
DATA$layer<-paste("L",substr(nemad.pourcent.bact$groupe,4,4),sep="")
DATA$core<-substr(nemad.pourcent.bact$groupe,3,3)
names(DATA)

bloc<-c(dim(nemad.pourcent.bact[,2:17])[2],
       dim(nemad.pourcent.detri[,2:27])[2],
       dim(nemad.pourcent.brou[,2:27])[2],
       dim(nemad.pourcent.pred[,2:8])[2])

#####################
# Multiple FActor Analysis (unsupervised)
#####################

# Extraction of eigen values
eig.nemad.pourcent.bact<-dudi.pca(nemad.pourcent.bact[,2:17],
                                  scannf=F,nf=2)$eig 
eig.nemad.pourcent.detri<-dudi.pca(nemad.pourcent.detri[,2:27],
                                   scannf=F,nf=2)$eig
eig.nemad.pourcent.brou<-dudi.pca(nemad.pourcent.brou[,2:27],
                                  scannf=F,nf=2)$eig 
eig.nemad.pourcent.pred<-dudi.pca(nemad.pourcent.pred[,2:8],
                                  scannf=F,nf=2)$eig

# MFA
names(DATA)
MFA<-dudi.pca(DATA[,1:75],col.w=rep(c(1/eig.nemad.pourcent.bact[1],
                                      1/eig.nemad.pourcent.detri[1],
                                      1/eig.nemad.pourcent.brou[1],
                                      1/eig.nemad.pourcent.pred[1]),
                                    bloc),
              scannf=F,
              nf=2)

varexp<-MFA$eig*100/sum(MFA$eig)
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
     xlab=paste("Axis 1 : ",round(varexp[1],2),"%"),
     ylab=paste("Axis 2 : ",round(varexp[2],2),"%"),
     main="MFA scores")

ordihull(MFA$li,fact2,lab=T)

plot(BCA$ls[,2]~BCA$ls[,1],
     pch=21,cex=2,col="white",
     bg=rainbow(length(levels(fact)))[fact],
     xlab=paste("Axis 1 : ",round(varexp2[1],2),"%"),
     ylab=paste("Axis 2 : ",round(varexp2[2],2),"%"),
     main="BCA scores")

legend("topleft",
       c("instrumental variable = Sites",
         "Total inertia explained = 14.3%",
        "p=0.001"),cex=0.8,box.lty=0)

ordihull(BCA$ls,fact,lab=T)

plot(BCA$co[,2]~BCA$co[,1],type="n",
     xlab=paste("Axis 1 : ",round(varexp2[1],2),"%"),
     ylab=paste("Axis 2 : ",round(varexp2[2],2),"%"),
     ylim=c(-1,1),xlim=c(-1,1),
     main="BCA loadings")

arrows(x0=0,y0=0,x1=BCA$co[,1],y1=BCA$co[,2],col="lightgrey",length=0.1)

cos2 <- as.matrix(BCA$co[,1:2])*as.matrix(BCA$co[,1:2])

var.group<-factor(rep(c("bac","detri","brou","pred"),bloc),
                  levels=c("bac","detri","brou","pred"))
BCA$co$col<-my.palette(length(bloc))[var.group]

newco <- BCA$co[cos2[,1]>0.2 | cos2[,2]>0.2,] 
oldco <- BCA$co[cos2[,1]<0.2 & cos2[,2]<0.2,]

text(newco[,1], newco[,2],rownames(newco),col=newco$col,cex=1.2)
text(oldco[,1], oldco[,2],rownames(oldco),col=alpha(oldco$col,0.3),cex=0.8) # la valeur alpha contrôle la transparence des variables < au seuil de cos2
summary(BCA)
abline(h=0,lty="dashed")
abline(v=0,lty="dashed")

legend("topleft",
       c("Bacterivorous",
         "Detritivores",
         "Grazers",
         "Omnivorous-Predators"),
       text.col=my.palette(length(bloc)),cex=0.8,box.lty=0)
