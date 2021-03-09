#ingédients#
library(ggplot2)#
library(ade4)#
library(vegan)#
library(RVAideMemoire)#
library(emmeans)#
library(pls)#
library(plotrix)#
library(wordcloud)

# crÃ©ation d'une palette de couleurs personalisÃ©e#
my.palette=colorRampPalette(c("black","red3","#009E73","blue"))

# Importation des tableaux :#
###DONNEES EN POURCENTAGES##

##POURCENT NEMAD BACT##
nemad.pourcent.bact=read.csv("C:/Users/clair/OneDrive/Bureau/papier-meio/Nematodes_pourcent_bacterivores.csv",h=T,sep=";",dec=",")
#nemad.pourcent.bact=read.csv("Nematodes_pourcent_bacterivores.csv",h=T,sep=";",dec=",")
names(nemad.pourcent.bact)
#POURCENT NEMAD DETRITIVORES#
nemad.pourcent.detri=read.csv("C:/Users/clair/OneDrive/Bureau/papier-meio/Nematodes_pourcent_detritivores.csv",h=T,sep=";",dec=",")
#nemad.pourcent.detri=read.csv("Nematodes_pourcent_detritivores.csv",h=T,sep=";",dec=",")
names(nemad.pourcent.detri)
#POURCENT NEMAD BROUTEURS#
nemad.pourcent.brou=read.csv("C:/Users/clair/OneDrive/Bureau/papier-meio/Nematodes_pourcent_brouteurs.csv",h=T,sep=";",dec=",")
#nemad.pourcent.brou=read.csv("Nematodes_pourcent_brouteurs.csv",h=T,sep=";",dec=",")
names(nemad.pourcent.brou)
#POURCENT NEMAD PREDATEURS#
nemad.pourcent.pred=read.csv("C:/Users/clair/OneDrive/Bureau/papier-meio/Nematodes_pourcent_predateurs.csv",h=T,sep=";",dec=",")
#nemad.pourcent.pred=read.csv("Nematodes_pourcent_predateurs.csv",h=T,sep=";",dec=",")
names(nemad.pourcent.pred)



##POURCENT NEMAB BACT##
nemab.pourcent.bact=read.csv("C:/Users/clair/OneDrive/Bureau/papier-meio/Nematodes_biomasscm3_pourcent_bacterivores.csv",h=T,sep=";",dec=",")
#nemab.pourcent.bact=read.csv("Nematodes_biomasscm3_pourcent_bacterivores.csv",h=T,sep=";",dec=",")
names(nemab.pourcent.bact)
#POURCENT NEMAB DETRITIVORES#
nemab.pourcent.detri=read.csv("C:/Users/clair/OneDrive/Bureau/papier-meio/Nematodes_biomasscm3_pourcent_detritivores.csv",h=T,sep=";",dec=",")
#nemab.pourcent.detri=read.csv("Nematodes_biomasscm3_pourcent_detritivores.csv",h=T,sep=";",dec=",")
names(nemab.pourcent.detri)
#POURCENT NEMAB BROUTEURS#
nemab.pourcent.brou=read.csv("C:/Users/clair/OneDrive/Bureau/papier-meio/Nematodes_biomasscm3_pourcent_brouteurs.csv",h=T,sep=";",dec=",")
#nemab.pourcent.brou=read.csv("Nematodes_biomasscm3_pourcent_brouteurs.csv",h=T,sep=";",dec=",")
names(nemab.pourcent.brou)
#POURCENT NEMAB PREDATEURS#
nemab.pourcent.pred=read.csv("C:/Users/clair/OneDrive/Bureau/papier-meio/Nematodes_biomasscm3_pourcent_predateurs.csv",h=T,sep=";",dec=",")
#nemab.pourcent.pred=read.csv("Nematodes_biomasscm3_pourcent_predateurs.csv",h=T,sep=";",dec=",")
names(nemab.pourcent.pred)



# verification des individus statistiques (doivent être appariés)#
data.frame(enz[,1],AG2[,1],conta[,1],pigments[,1],enviro[,1],bact[,1],gconta[,1],AG.g[,1])


# Standardisation#
vp.nemad.bact=dudi.pca(nemad.bact[,2:17],scannf=F,nf=2)$eig # extraction des valeurs propres#
nemad.bact.std=nemad.bact[,2:17]/vp.nemad.bact[1] # pondération du tableau contaminants par la première vp#
vp.nemad.detri=dudi.pca(nemad.detri[,2:27],scannf=F,nf=2)$eig # extraction des valeurs propres#
nemad.detri.std=nemad.detri[,2:27]/vp.nemad.detri[1] # pondération du tableau contaminants par la première vp#
vp.nemad.brou=dudi.pca(nemad.brou[,2:27],scannf=F,nf=2)$eig # extraction des valeurs propres#
nemad.brou.std=nemad.brou[,2:27]/vp.nemad.brou[1] # pondération du tableau contaminants par la première vp#
vp.nemad.pred=dudi.pca(nemad.pred[,2:8],scannf=F,nf=2)$eig # extraction des valeurs propres#
nemad.pred.std=nemad.pred[,2:8]/vp.nemad.pred[1] # pondération du tableau contaminants par la première vp#

vp.nemab.bact=dudi.pca(nemab.bact[,2:17],scannf=F,nf=2)$eig # extraction des valeurs propres#
nemab.bact.std=nemab.bact[,2:17]/vp.nemab.bact[1] # pondération du tableau contaminants par la première vp#
vp.nemab.detri=dudi.pca(nemab.detri[,2:27],scannf=F,nf=2)$eig # extraction des valeurs propres#
nemab.detri.std=nemab.detri[,2:27]/vp.nemab.detri[1] # pondération du tableau contaminants par la première vp#
vp.nemab.brou=dudi.pca(nemab.brou[,2:26],scannf=F,nf=2)$eig # extraction des valeurs propres#
nemab.brou.std=nemab.brou[,2:26]/vp.nemab.brou[1] # pondération du tableau contaminants par la première vp#
vp.nemab.pred=dudi.pca(nemab.pred[,2:8],scannf=F,nf=2)$eig # extraction des valeurs propres#
nemab.pred.std=nemab.pred[,2:8]/vp.nemab.pred[1] # pondération du tableau contaminants par la première vp#

##NEMA POURCENTAGE DENSITY GT##
vp.nemad.pourcent.bact=dudi.pca(nemad.pourcent.bact[,2:17],scannf=F,nf=2)$eig # extraction des valeurs propres#
nemad.pourcent.bact.std=nemad.pourcent.bact[,2:17]/vp.nemad.pourcent.bact[1] # pondération du tableau contaminants par la première vp#
vp.nemad.pourcent.detri=dudi.pca(nemad.pourcent.detri[,2:27],scannf=F,nf=2)$eig # extraction des valeurs propres#
nemad.pourcent.detri.std=nemad.pourcent.detri[,2:27]/vp.nemad.pourcent.detri[1] # pondération du tableau contaminants par la première vp#
vp.nemad.pourcent.brou=dudi.pca(nemad.pourcent.brou[,2:27],scannf=F,nf=2)$eig # extraction des valeurs propres#
nemad.pourcent.brou.std=nemad.pourcent.brou[,2:27]/vp.nemad.pourcent.brou[1] # pondération du tableau contaminants par la première vp#
vp.nemad.pourcent.pred=dudi.pca(nemad.pourcent.pred[,2:8],scannf=F,nf=2)$eig # extraction des valeurs propres#
nemad.pourcent.pred.std=nemad.pourcent.pred[,2:8]/vp.nemad.pourcent.pred[1] # pondération du tableau contaminants par la première vp#

##NEMA POURCENTAGE BIOMASS GT##
vp.nemab.pourcent.bact=dudi.pca(nemab.pourcent.bact[,2:17],scannf=F,nf=2)$eig # extraction des valeurs propres#
nemab.pourcent.bact.std=nemab.pourcent.bact[,2:17]/vp.nemab.pourcent.bact[1] # pondération du tableau contaminants par la première vp#
vp.nemab.pourcent.detri=dudi.pca(nemab.pourcent.detri[,2:27],scannf=F,nf=2)$eig # extraction des valeurs propres#
nemab.pourcent.detri.std=nemab.pourcent.detri[,2:27]/vp.nemab.pourcent.detri[1] # pondération du tableau contaminants par la première vp#
vp.nemab.pourcent.brou=dudi.pca(nemab.pourcent.brou[,2:26],scannf=F,nf=2)$eig # extraction des valeurs propres#
nemab.pourcent.brou.std=nemab.pourcent.brou[,2:26]/vp.nemab.pourcent.brou[1] # pondération du tableau contaminants par la première vp#
vp.nemab.pourcent.pred=dudi.pca(nemab.pourcent.pred[,2:8],scannf=F,nf=2)$eig # extraction des valeurs propres#
nemab.pourcent.pred.std=nemab.pourcent.pred[,2:8]/vp.nemab.pourcent.pred[1] # pondération du tableau contaminants par la première vp#


# final dataset#
DATA=cbind(nemab.pourcent.bact.std,nemab.pourcent.detri.std,nemab.pourcent.brou.std,nemab.pourcent.pred.std)#
DATA$sites=substr(nemad.pourcent.bact$groupe,1,2)#
DATA$horizon=substr(nemad.pourcent.bact$groupe,4,4)#
DATA$core=substr(nemad.pourcent.bact$groupe,3,3)#
names(DATA)#
bloc=c(dim(nemab.pourcent.bact.std)[2],dim(nemab.pourcent.detri.std)[2],dim(nemab.pourcent.brou.std)[2],dim(nemab.pourcent.pred.std)[2])

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# pca dataset
DATA2=cbind(nemab.pourcent.bact[,2:17],nemab.pourcent.detri[,2:27],nemab.pourcent.brou[,2:26],nemab.pourcent.pred[,2:8])#
DATA2$sites=substr(nemad.pourcent.bact$groupe,1,2)#
DATA2$horizon=substr(nemad.pourcent.bact$groupe,4,4)#
DATA2$core=substr(nemad.pourcent.bact$groupe,3,3)#
names(DATA2)#
bloc=c(dim(nemab.pourcent.bact[,2:17])[2],dim(nemab.pourcent.detri[,2:27])[2],dim(nemab.pourcent.brou[,2:26])[2],dim(nemab.pourcent.pred[,2:8])[2])
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

# MFA new script
names(DATA2)
par(layout(matrix(c(1,2,3,3),nrow=2)),mar=c(4,4,4,4))
MFA2=dudi.pca(DATA2[,1:74],col.w=rep(c(1/vp.nemad.pourcent.bact[1],1/vp.nemad.pourcent.detri[1],1/vp.nemad.pourcent.brou[1],1/vp.nemad.pourcent.pred[1]),bloc),scannf=F,nf=2)
varexp1=MFA2$eig*100/sum(MFA2$eig)#
facteur1=factor(paste(DATA2$site))#
facteur2=factor(paste(DATA2$site,DATA2$horizon))#
facteur3=factor(paste(DATA2$horizon))#
plot(MFA2$li[,2]~MFA2$li[,1],pch=21,cex=2,col="white",bg=rainbow(length(levels(facteur1)))[facteur1],xlab=paste("Axis 1 : ",round(varexp1[1],2),"%"),ylab=paste("Axis 2 : ",round(varexp1[2],2),"%"),main="MFA scores")#
ordihull(MFA2$li,facteur2,lab=T)#

# Between class analysis (analyse supervisée - FACTEUR 1 SITES#
BCA2=bca(MFA2,facteur1,scannf=F,nf=2)#
BCA2$ratio#
varexp2=BCA2$eig*100/sum(BCA2$eig)#
plot(BCA2$ls[,2]~BCA2$ls[,1],pch=21,cex=2,col="white",bg=rainbow(length(levels(facteur1)))[facteur1],xlab=paste("Axis 1 : ",round(varexp2[1],2),"%"),ylab=paste("Axis 2 : ",round(varexp2[2],2),"%"),main="BCA scores")#
legend("topleft",c("instrumental variable = Sites","Total inertia explained = 12.5%"),cex=0.8,box.lty=0)
ordihull(BCA2$ls,facteur1,lab=T)#
plot(BCA2$co[,2]~BCA2$co[,1],type="n",xlab=paste("Axis 1 : ",round(varexp2[1],2),"%"),ylab=paste("Axis 2 : ",round(varexp2[2],2),"%"),ylim=c(-1,1),xlim=c(-1,1),main="BCA loadings")#
BCA2$co$col=rep(my.palette(4),bloc)#
newco=BCA2$co[abs(BCA2$co[,1])>0 | abs(BCA2$co[,2])>0,]#
arrows(x0=0,y0=0,x1=BCA2$co[,1],y1=BCA2$co[,2],col="lightgrey",length=0.1)#
cos2 <- as.matrix(BCA2$co[,1:2])*as.matrix(BCA2$co[,1:2]) # calcul des cos2 pour filtrer les variables
newco <- BCA2$co[cos2[,1]>0.2 | cos2[,2]>0.2,] # utiliser une valeur de cos2 pour filtrer (ici 0.5)
oldco <- BCA2$co[cos2[,1]<0.2 & cos2[,2]<0.2,]
text(newco[,1], newco[,2],rownames(newco),col=newco$col,cex=1.7)
text(oldco[,1], oldco[,2],rownames(oldco),col=alpha(oldco$col,0.3),cex=1.5) # la valeur alpha contrÃ´le la transparence des variables < au seuil de cos2
summary(BCA2) ### AJOUT
#draw.circle(0,0,1)#
#draw.circle(0,0,0.3,lty="dashed")#
#text(newco[,1], newco[,2],rownames(newco),col=newco[,1],cex=0.6) # 4 = nombre de martrices dans le tableau DATA#
abline(h=0,lty="dashed")
abline(v=0,lty="dashed")
 # 4 = nombre de martrices dans le tableau DATA#
legend("topleft",c("Bacterivorous", "Detritivores", "Grazers","Omnivorous-Predators"),text.col=my.palette(length(bloc)),cex=1.2,box.lty=0)
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
