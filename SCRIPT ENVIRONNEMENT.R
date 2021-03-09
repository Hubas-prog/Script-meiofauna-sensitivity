
# ingÃ©dients#
library(ggplot2)#
library(ade4)#
library(vegan)#
library(RVAideMemoire)#
library(emmeans)#
library(pls)#
library(plotrix)#
library(factoextra)
#


# crÃ©ation d'une palette de couleurs personalisÃ©e#
my.palette=colorRampPalette(c("purple","black","red1","green3","blue"))


# Importation des tableaux :#


# Contaminants
conta=read.csv("C:/Users/clair/OneDrive/Bureau/papier-meio/Contaminants.csv",h=T,sep=";",dec=",")
#conta=read.csv("Contaminants.csv",h=T,sep=";",dec=",")

names(conta)

# Pigments
pigments=read.csv("C:/Users/clair/OneDrive/Bureau/papier-meio/pigments.csv",h=T,sep=";",dec=",")
#pigments=read.csv("pigments.csv",h=T,sep=";",dec=",")
names(pigments)

#Rapport C/N
CHN=read.csv("C:/Users/clair/OneDrive/Bureau/papier-meio/CHN.csv",h=T,sep=";",dec=",")
#CHN=read.csv("CHN.csv",h=T,sep=";",dec=",")
names(CHN)

#COMMUNAUTE MICROBIENNE
bact=read.csv("C:/Users/clair/OneDrive/Bureau/papier-meio/Communaute-Bacterienne.csv",h=T,sep=";",dec=",")
#bact=read.csv("Communaute-Bacterienne.csv",h=T,sep=";",dec=",")
names(bact)

#ENVIRONNEMENT (GRANULO,PH,REDOX)##
enviro=read.csv("C:/Users/clair/OneDrive/Bureau/papier-meio/enviro.csv",h=T,sep=";",dec=",")
#enviro=read.csv("enviro.csv",h=T,sep=";",dec=",")
names(enviro)

# verification des individus statistiques (doivent Ãªtre appariÃ©s)#
data.frame(conta[,1],pigments[,1],enviro[,1],bact[,1],CHN[,1])


# Standardisation#
vp.conta=dudi.pca(conta[,2:32],scannf=F,nf=2)$eig # extraction des valeurs propres
conta.std=conta[,2:32]/vp.conta[1] # pondÃ©ration du tableau contaminants par la premiÃ¨re vp
vp.pig=dudi.pca(pigments[,2:19],scannf=F,nf=2)$eig # extraction des valeurs propres
pig.std=pigments[,2:19]/vp.pig[1] # pondÃ©ration du tableau contaminants par la premiÃ¨re vp
vp.bact=dudi.pca(bact[,3:4],scannf=F,nf=2)$eig # extraction des valeurs propres
bact.std=bact[,3:4]/vp.bact[1] # pondÃ©ration du tableau nematodes par la premiÃ¨re vp
vp.enviro=dudi.pca(enviro[,2:6],scannf=F,nf=2)$eig # extraction des valeurs propres#
enviro.std=enviro[,2:6]/vp.enviro[1] # pondÃ©ration du tableau contaminants par la premiÃ¨re vp#
vp.CHN=dudi.pca(CHN[,2:3],scannf=F,nf=2)$eig # extraction des valeurs propres#
CHN.std=CHN[,2:3]/vp.CHN[1] # pondÃ©ration du tableau contaminants par la premiÃ¨re vp#

# final dataset#
DATA=cbind(enviro.std,conta.std,bact.std,CHN.std,pig.std)#
DATA$sites=substr(enviro$groupe,1,2)#
DATA$horizon=substr(enviro$groupe,4,4)#
DATA$core=substr(enviro$groupe,3,3)#
names(DATA)#
bloc=c(dim(enviro.std)[2],dim(conta.std)[2],dim(bact.std)[2],dim(CHN.std)[2],dim(pig.std)[2])



###### All data#
# Multiple factor analysis (non supervisÃ©e)#
par(layout(matrix(c(1,2,3,3),nrow=2)),mar=c(4,4,4,4))
MFA=dudi.pca(DATA[,1:58],scannf=F,nf=2)#
varexp1=MFA$eig*100/sum(MFA$eig)#
facteur1=factor(paste(DATA$site))#
facteur2=factor(paste(DATA$site,DATA$horizon))#
facteur3=factor(paste(DATA$horizon))#
plot(MFA$li[,2]~MFA$li[,1],pch=21,cex=2,col="white",bg=rainbow(length(levels(facteur1)))[facteur1],xlab=paste("Axis 1 : ",round(varexp1[1],2),"%"),ylab=paste("Axis 2 : ",round(varexp1[2],2),"%"),main="MFA scores")#
ordihull(MFA$li,facteur2,lab=T)#
#plot(MFA$co[,2]~MFA$co[,1],type="n",xlab=paste("Axis 1 : ",round(varexp1[1],2),"%"),ylab=paste("Axis 2 : ",round(varexp1[2],2),"%"),xlim=c(-1,1),ylim=c(-1,1))#
#MFA$co$col=rep(my.palette(5),bloc)#
#newco.MFA=MFA$co[abs(MFA$co[,1])>0 | abs(MFA$co[,2])>0,]#
#arrows(x0=0,y0=0,x1=newco.MFA[,1],y1=newco.MFA[,2],col="lightgrey",length=0.1)#
#draw.circle(0,0,1)#
#draw.circle(0,0,0.5,lty="dashed")#
#text(newco.MFA[,1], newco.MFA[,2],rownames(newco.MFA),col=rep(my.palette(length(bloc)),bloc),cex=0.6) # 4 = nombre de martrices dans le tableau DATA#



# Between class analysis (analyse supervisÃ©e - FACTEUR 1 SITES#
BCA=bca(MFA,facteur1,scannf=F,nf=2)#
BCA$ratio#
varexp2=BCA$eig*100/sum(BCA$eig)#
plot(BCA$ls[,2]~BCA$ls[,1],pch=21,cex=2,col="white",bg=rainbow(length(levels(facteur1)))[facteur1],xlab=paste("Axis 1 : ",round(varexp2[1],2),"%"),ylab=paste("Axis 2 : ",round(varexp2[2],2),"%"),main="BCA scores")#
legend("topright",c("instrumental variable = Sites","Total inertia explained = 34.5%"),cex=0.8,box.lty=0)
ordihull(BCA$ls,facteur1,lab=T)#
plot(BCA$co[,2]~BCA$co[,1],type="n",xlab=paste("Axis 1 : ",round(varexp2[1],2),"%"),ylab=paste("Axis 2 : ",round(varexp2[2],2),"%"),ylim=c(-1,1),xlim=c(-1,1),main="BCA loadings")#
BCA$co$col=rep(my.palette(5),bloc)
arrows(x0=0,y0=0,x1=BCA$co[,1],y1=BCA$co[,2],col="lightgrey",length=0.1)#
cos2 <- as.matrix(BCA$co[,1:2])*as.matrix(BCA$co[,1:2]) # calcul des cos2 pour filtrer les variables
newco <- BCA$co[cos2[,1]>0.4 | cos2[,2]>0.4,] # utiliser une valeur de cos2 pour filtrer (ici 0.5)
oldco <- BCA$co[cos2[,1]<0.4 & cos2[,2]<0.4,]
text(newco[,1], newco[,2],rownames(newco),col=newco$col,cex=1.2)
text(oldco[,1], oldco[,2],rownames(oldco),col=alpha(oldco$col,0.3),cex=1) # la valeur alpha contrÃ´le la transparence des variables < au seuil de cos2
summary(BCA) ### AJOUT
#draw.circle(0,0,1)#
#draw.circle(0,0,0.3,lty="dashed")#
abline(h=0,lty="dashed")
abline(v=0,lty="dashed")
 # 4 = nombre de martrices dans le tableau DATA#
legend("topleft",c("Physicochemical parameters", "Contaminants", "Prokaryotic biomass","Organic Matter", "Pigments"),text.col=my.palette(length(bloc)),cex=1.1,box.lty=0)



# Between class analysis (analyse supervisÃ©e - FACTEUR 2 SITES ET HORIZONS#
BCA=bca(MFA,facteur2,scannf=F,nf=2)#
BCA$ratio#
varexp2=BCA$eig*100/sum(BCA$eig)#
plot(BCA$ls[,2]~BCA$ls[,1],pch=21,cex=2,col="white",bg=rainbow(length(levels(facteur1)))[facteur1],xlab=paste("Axis 1 : ",round(varexp2[1],2),"%"),ylab=paste("Axis 2 : ",round(varexp2[2],2),"%"),main="BCA scores")#
legend("topleft",c("Interaction sites:sediment depth","Total inertia explained = 60%"),cex=0.8,box.lty=0)
ordihull(BCA$ls,facteur2,lab=T)#
plot(BCA$co[,2]~BCA$co[,1],type="n",xlab=paste("Axis 1 : ",round(varexp2[1],2),"%"),ylab=paste("Axis 2 : ",round(varexp2[2],2),"%"),ylim=c(-1,1),xlim=c(-1,1),main="BCA loadings")#
BCA$co$col=rep(my.palette(5),bloc)#
newco=BCA$co[abs(BCA$co[,1])>0 | abs(BCA$co[,2])>0,]#
arrows(x0=0,y0=0,x1=newco[,1],y1=newco[,2],col="lightgrey",length=0.1)#
text(newco[,1], newco[,2],rownames(newco),col= rep(my.palette(5),bloc),cex=0.8)
summary(BCA) ### AJOUT
#draw.circle(0,0,1)#
#draw.circle(0,0,0.3,lty="dashed")#
abline(h=0,lty="dashed")
abline(v=0,lty="dashed")
# 4 = nombre de martrices dans le tableau DATA#
legend("topleft",c("Environmental variables", "Pollutants", "Procaryotes","CHN", "Pigments"),text.col=my.palette(length(bloc)),cex=0.8,box.lty=0)


# Between class analysis (analyse supervisÃ©e - FACTEUR 3 HORIZON#
BCA=bca(MFA,facteur3,scannf=F,nf=2)#
BCA$ratio#
varexp2=BCA$eig*100/sum(BCA$eig)#
par(mfrow=c(1,2))#
plot(BCA$ls[,2]~BCA$ls[,1],pch=21,cex=2,col="white",bg=rainbow(length(levels(facteur3)))[facteur3],xlab=paste("Axis 1 : ",round(varexp2[1],2),"%"),ylab=paste("Axis 2 : ",round(varexp2[2],2),"%"))#
ordihull(BCA$ls,facteur3,lab=T)#
plot(BCA$co[,2]~BCA$co[,1],type="n",xlab=paste("Axis 1 : ",round(varexp2[1],2),"%"),ylab=paste("Axis 2 : ",round(varexp2[2],2),"%"),ylim=c(-1,1),xlim=c(-1,1))#
BCA$co$col=rep(my.palette(2),bloc)#
newco=BCA$co[abs(BCA$co[,1])>0 | abs(BCA$co[,2])>0,]#
text(newco[,1], newco[,2],rownames(newco),col= rep(my.palette(2),bloc))
arrows(x0=0,y0=0,x1=newco[,1],y1=newco[,2],col="lightgrey",length=0.1)#
text(newco[,1], newco[,2],rownames(newco),col= rep(my.palette(2),bloc))## AJOUT
summary(BCA) ### AJOUT
draw.circle(0,0,1)#
draw.circle(0,0,0.3,lty="dashed")#
text(newco[,1], newco[,2],rownames(newco),col=newco[,2],cex=0.6) # 4 = nombre de martrices dans le tableau DATA#
#


############ MFA as weighted PCA

par(layout(matrix(c(1,2,3,3),nrow=2)),mar=c(4,4,4,4))
MFA2=dudi.pca(DATA[,1:58],scannf=F,nf=2)
varexp1=MFA2$eig*100/sum(MFA2$eig)#
facteur1=factor(paste(DATA$site))#
facteur2=factor(paste(DATA$site,DATA$horizon))#
facteur3=factor(paste(DATA$horizon))#
plot(MFA2$li[,2]~MFA2$li[,1],pch=21,cex=2,col="white",bg=rainbow(length(levels(facteur1)))[facteur1],xlab=paste("Axis 1 : ",round(varexp1[1],2),"%"),ylab=paste("Axis 2 : ",round(varexp1[2],2),"%"),main="MFA scores")#
ordihull(MFA2$li,facteur2,lab=T)#

# Between class analysis (analyse supervisée - FACTEUR 1 SITES#
BCA2=bca(MFA2,facteur1,scannf=F,nf=2)#
BCA2$ratio#
varexp2=BCA$eig*100/sum(BCA$eig)#
plot(BCA2$ls[,2]~BCA2$ls[,1],pch=21,cex=2,col="white",bg=rainbow(length(levels(facteur1)))[facteur1],xlab=paste("Axis 1 : ",round(varexp2[1],2),"%"),ylab=paste("Axis 2 : ",round(varexp2[2],2),"%"),main="BCA scores")#
legend("topright",c("instrumental variable = Sites","Total inertia explained = 38%"),cex=0.8,box.lty=0)
ordihull(BCA2$ls,facteur1,lab=T)#
plot(BCA2$co[,2]~BCA2$co[,1],type="n",xlab=paste("Axis 1 : ",round(varexp2[1],2),"%"),ylab=paste("Axis 2 : ",round(varexp2[2],2),"%"),ylim=c(-1,1),xlim=c(-1,1),main="BCA loadings")#
BCA2$co$col=rep(rainbow(5),bloc)
arrows(x0=0,y0=0,x1=BCA2$co[,1],y1=BCA2$co[,2],col="lightgrey",length=0.1)#
cos2 <- as.matrix(BCA2$co[,1:2])*as.matrix(BCA2$co[,1:2]) # calcul des cos2 pour filtrer les variables
newco <- BCA2$co[cos2[,1]>0.5 | cos2[,2]>0.5,] # utiliser une valeur de cos2 pour filtrer (ici 0.5)
oldco <- BCA2$co[cos2[,1]<0.5 & cos2[,2]<0.5,]
text(newco[,1], newco[,2],rownames(newco),col=newco$col,cex=1.2)
text(oldco[,1], oldco[,2],rownames(oldco),col=alpha(oldco$col,0.3),cex=1) # la valeur alpha contrôle la transparence des variables < au seuil de cos2
summary(BCA2) ### AJOUT
#draw.circle(0,0,1)#
#draw.circle(0,0,0.3,lty="dashed")#~/Desktop/plot-old-script.pdf
abline(h=0,lty="dashed")
abline(v=0,lty="dashed")
 # 4 = nombre de martrices dans le tableau DATA#
legend("topleft",c("Physicochemical parameters", "Contaminants", "Prokaryotic biomass","Organic Matter", "Pigments"),text.col=my.palette(length(bloc)),cex=1.1,box.lty=0)


# Multiple factor analysis (non supervisée)#
par(layout(matrix(c(1,2,3,3),nrow=2)),mar=c(4,4,4,4))
MFA=dudi.pca(DATA[,1:58],scannf=F,nf=2)#
varexp1=MFA$eig*100/sum(MFA$eig)#
facteur1=factor(paste(DATA$site))#
facteur2=factor(paste(DATA$site,DATA$horizon))#
facteur3=factor(paste(DATA$horizon))#
plot(MFA$li[,2]~MFA$li[,1],pch=21,cex=2,col="white",bg=rainbow(length(levels(facteur1)))[facteur1],xlab=paste("Axis 1 : ",round(varexp1[1],2),"%"),ylab=paste("Axis 2 : ",round(varexp1[2],2),"%"),main="MFA scores")#
ordihull(MFA$li,facteur2,lab=T)#
#plot(MFA$co[,2]~MFA$co[,1],type="n",xlab=paste("Axis 1 : ",round(varexp1[1],2),"%"),ylab=paste("Axis 2 : ",round(varexp1[2],2),"%"),xlim=c(-1,1),ylim=c(-1,1))#
#MFA$co$col=rep(my.palette(5),bloc)#
#newco.MFA=MFA$co[abs(MFA$co[,1])>0 | abs(MFA$co[,2])>0,]#
#arrows(x0=0,y0=0,x1=newco.MFA[,1],y1=newco.MFA[,2],col="lightgrey",length=0.1)#
#draw.circle(0,0,1)#
#draw.circle(0,0,0.5,lty="dashed")#
#text(newco.MFA[,1], newco.MFA[,2],rownames(newco.MFA),col=rep(my.palette(length(bloc)),bloc),cex=0.6) # 4 = nombre de martrices dans le tableau DATA#

# Between class analysis (analyse supervisée - FACTEUR 1 SITES#
BCA=bca(MFA,facteur1,scannf=F,nf=2)#
BCA$ratio#
varexp2=BCA$eig*100/sum(BCA$eig)#
plot(BCA$ls[,2]~BCA$ls[,1],pch=21,cex=2,col="white",bg=rainbow(length(levels(facteur1)))[facteur1],xlab=paste("Axis 1 : ",round(varexp2[1],2),"%"),ylab=paste("Axis 2 : ",round(varexp2[2],2),"%"),main="BCA scores")#
legend("topright",c("instrumental variable = Sites","Total inertia explained = 36%"),cex=0.8,box.lty=0)
ordihull(BCA$ls,facteur1,lab=T)#
plot(BCA$co[,2]~BCA$co[,1],type="n",xlab=paste("Axis 1 : ",round(varexp2[1],2),"%"),ylab=paste("Axis 2 : ",round(varexp2[2],2),"%"),ylim=c(-1,1),xlim=c(-1,1),main="BCA loadings")#
BCA$co$col=rep(rainbow(5),bloc)
arrows(x0=0,y0=0,x1=BCA$co[,1],y1=BCA$co[,2],col="lightgrey",length=0.1)#
cos2 <- as.matrix(BCA$co[,1:2])*as.matrix(BCA$co[,1:2]) # calcul des cos2 pour filtrer les variables
newco <- BCA$co[cos2[,1]>0.5 | cos2[,2]>0.5,] # utiliser une valeur de cos2 pour filtrer (ici 0.5)
oldco <- BCA$co[cos2[,1]<0.5 & cos2[,2]<0.5,]
text(newco[,1], newco[,2],rownames(newco),col=newco$col,cex=1.2)
text(oldco[,1], oldco[,2],rownames(oldco),col=alpha(oldco$col,0.3),cex=1) # la valeur alpha contrôle la transparence des variables < au seuil de cos2
summary(BCA) ### AJOUT
#draw.circle(0,0,1)#
#draw.circle(0,0,0.3,lty="dashed")#
abline(h=0,lty="dashed")
abline(v=0,lty="dashed")
 # 4 = nombre de martrices dans le tableau DATA#
legend("topleft",c("Environmental variables", "Pollutants", "Procaryotes","CHN", "Pigments"),text.col=rainbow(length(bloc)),cex=0.8,box.lty=0)

################
# Test BCA
randtest(BCA)
###############


# Multiple factor analysis (non supervisée)#
par(layout(matrix(c(1,2,3,3),nrow=2)),mar=c(4,4,4,4))
MFA=dudi.pca(DATA[,1:58],scannf=F,nf=2)#
varexp1=MFA$eig*100/sum(MFA$eig)#
facteur1=factor(paste(DATA$site))#
facteur2=factor(paste(DATA$site,DATA$horizon))#
facteur3=factor(paste(DATA$horizon))#
plot(MFA$li[,2]~MFA$li[,1],pch=21,cex=2,col="white",bg=rainbow(length(levels(facteur1)))[facteur1],xlab=paste("Axis 1 : ",round(varexp1[1],2),"%"),ylab=paste("Axis 2 : ",round(varexp1[2],2),"%"),main="MFA scores")#
ordihull(MFA$li,facteur2,lab=T)#
#plot(MFA$co[,2]~MFA$co[,1],type="n",xlab=paste("Axis 1 : ",round(varexp1[1],2),"%"),ylab=paste("Axis 2 : ",round(varexp1[2],2),"%"),xlim=c(-1,1),ylim=c(-1,1))#
#MFA$co$col=rep(my.palette(5),bloc)#
#newco.MFA=MFA$co[abs(MFA$co[,1])>0 | abs(MFA$co[,2])>0,]#
#arrows(x0=0,y0=0,x1=newco.MFA[,1],y1=newco.MFA[,2],col="lightgrey",length=0.1)#
#draw.circle(0,0,1)#
#draw.circle(0,0,0.5,lty="dashed")#
#text(newco.MFA[,1], newco.MFA[,2],rownames(newco.MFA),col=rep(my.palette(length(bloc)),bloc),cex=0.6) # 4 = nombre de martrices dans le tableau DATA#

# Between class analysis (analyse supervisée - FACTEUR 1 SITES#
BCA=bca(MFA,facteur1,scannf=F,nf=2)#
BCA$ratio#
varexp2=BCA$eig*100/sum(BCA$eig)#
plot(BCA$ls[,2]~BCA$ls[,1],pch=21,cex=2,col="white",bg=rainbow(length(levels(facteur1)))[facteur1],xlab=paste("Axis 1 : ",round(varexp2[1],2),"%"),ylab=paste("Axis 2 : ",round(varexp2[2],2),"%"),main="BCA scores")#
legend("topright",c("instrumental variable = Sites","Total inertia explained = 36%"),cex=0.8,box.lty=0)
ordihull(BCA$ls,facteur1,lab=T)#
plot(BCA$co[,2]~BCA$co[,1],type="n",xlab=paste("Axis 1 : ",round(varexp2[1],2),"%"),ylab=paste("Axis 2 : ",round(varexp2[2],2),"%"),ylim=c(-1,1),xlim=c(-1,1),main="BCA loadings")#
BCA$co$col=rep(rainbow(5),bloc)
arrows(x0=0,y0=0,x1=BCA$co[,1],y1=BCA$co[,2],col="lightgrey",length=0.1)#
cos2 <- as.matrix(BCA$co[,1:2])*as.matrix(BCA$co[,1:2]) # calcul des cos2 pour filtrer les variables
newco <- BCA$co[cos2[,1]>0.5 | cos2[,2]>0.5,] # utiliser une valeur de cos2 pour filtrer (ici 0.5)
oldco <- BCA$co[cos2[,1]<0.5 & cos2[,2]<0.5,]
text(newco[,1], newco[,2],rownames(newco),col=newco$col,cex=1.2)
text(oldco[,1], oldco[,2],rownames(oldco),col=alpha(oldco$col,0.3),cex=1) # la valeur alpha contrôle la transparence des variables < au seuil de cos2
summary(BCA) ### AJOUT
#draw.circle(0,0,1)#
#draw.circle(0,0,0.3,lty="dashed")#
abline(h=0,lty="dashed")
abline(v=0,lty="dashed")
 # 4 = nombre de martrices dans le tableau DATA#
legend("topleft",c("Physicochemical parameters", "Contaminants", "Prokaryotic biomass","Organic Matter", "Pigments"),text.col=my.palette(length(bloc)),cex=1.1,box.lty=0)


-----------------------------------------------------------------------------------------------------------------------------

# ingédients#
library(ggplot2)#
library(ade4)#
library(vegan)#
library(RVAideMemoire)#
library(emmeans)#
library(pls)#
library(plotrix)#
library(factoextra)
#


# création d'une palette de couleurs personalisée#
my.palette=colorRampPalette(c("purple","black","red1","green3","blue"))


# Importation des tableaux :#

# Contaminants
conta=read.csv("C:/Users/clair/OneDrive/Bureau/papier-meio/Contaminants.csv",h=T,sep=";",dec=",")
#conta=read.csv("Contaminants.csv",h=T,sep=";",dec=",")

names(conta)

# Pigments
pigments=read.csv("C:/Users/clair/OneDrive/Bureau/papier-meio/pigments.csv",h=T,sep=";",dec=",")
#pigments=read.csv("pigments.csv",h=T,sep=";",dec=",")
names(pigments)

#Rapport C/N
CHN=read.csv("C:/Users/clair/OneDrive/Bureau/papier-meio/CHN.csv",h=T,sep=";",dec=",")
#CHN=read.csv("CHN.csv",h=T,sep=";",dec=",")
names(CHN)

#COMMUNAUTE MICROBIENNE
bact=read.csv("C:/Users/clair/OneDrive/Bureau/papier-meio/Communaute-Bacterienne.csv",h=T,sep=";",dec=",")
#bact=read.csv("Communaute-Bacterienne.csv",h=T,sep=";",dec=",")
names(bact)

#ENVIRONNEMENT (GRANULO,PH,REDOX)##
enviro=read.csv("C:/Users/clair/OneDrive/Bureau/papier-meio/enviro.csv",h=T,sep=";",dec=",")
#enviro=read.csv("enviro.csv",h=T,sep=";",dec=",")
names(enviro)

# verification des individus statistiques (doivent Ãªtre appariÃ©s)#
data.frame(conta[,1],pigments[,1],enviro[,1],bact[,1],CHN[,1])


# Standardisation#
vp.conta=dudi.pca(conta[,2:32],scannf=F,nf=2)$eig # extraction des valeurs propres
conta.std=conta[,2:32]/vp.conta[1] # pondération du tableau contaminants par la première vp
vp.pig=dudi.pca(pigments[,2:19],scannf=F,nf=2)$eig # extraction des valeurs propres
pig.std=pigments[,2:19]/vp.pig[1] # pondération du tableau contaminants par la première vp
vp.bact=dudi.pca(bact[,3:4],scannf=F,nf=2)$eig # extraction des valeurs propres
bact.std=bact[,3:4]/vp.bact[1] # pondération du tableau nematodes par la première vp
vp.enviro=dudi.pca(enviro[,2:6],scannf=F,nf=2)$eig # extraction des valeurs propres#
enviro.std=enviro[,2:6]/vp.enviro[1] # pondération du tableau contaminants par la première vp#
vp.CHN=dudi.pca(CHN[,2:3],scannf=F,nf=2)$eig # extraction des valeurs propres#
CHN.std=CHN[,2:3]/vp.CHN[1] # pondération du tableau contaminants par la première vp#

# final dataset#
DATA=cbind(enviro.std,conta.std,bact.std,CHN.std,pig.std)#
DATA$sites=substr(enviro$groupe,1,2)#
DATA$horizon=substr(enviro$groupe,4,4)#
DATA$core=substr(enviro$groupe,3,3)#
names(DATA)#
bloc=c(dim(enviro.std)[2],dim(conta.std)[2],dim(bact.std)[2],dim(CHN.std)[2],dim(pig.std)[2])

# PCA dataset
pca.data <- cbind(enviro[,2:6],conta[,2:32],bact[,3:4],CHN[,2:3],pigments[,2:19])
dim(pca.data)
sum(bloc)

###### All data#
# Multiple factor analysis (non supervisée)#
par(layout(matrix(c(1,2,3,3),nrow=2)),mar=c(4,4,4,4))
MFA=dudi.pca(DATA[,1:58],scannf=F,nf=2)#
varexp1=MFA$eig*100/sum(MFA$eig)#
facteur1=factor(paste(DATA$site))#
facteur2=factor(paste(DATA$site,DATA$horizon))#
facteur3=factor(paste(DATA$horizon))#
plot(MFA$li[,2]~MFA$li[,1],pch=21,cex=2,col="white",bg=rainbow(length(levels(facteur1)))[facteur1],xlab=paste("Axis 1 : ",round(varexp1[1],2),"%"),ylab=paste("Axis 2 : ",round(varexp1[2],2),"%"),main="MFA scores")#
ordihull(MFA$li,facteur2,lab=T)#
#plot(MFA$co[,2]~MFA$co[,1],type="n",xlab=paste("Axis 1 : ",round(varexp1[1],2),"%"),ylab=paste("Axis 2 : ",round(varexp1[2],2),"%"),xlim=c(-1,1),ylim=c(-1,1))#
#MFA$co$col=rep(my.palette(5),bloc)#
#newco.MFA=MFA$co[abs(MFA$co[,1])>0 | abs(MFA$co[,2])>0,]#
#arrows(x0=0,y0=0,x1=newco.MFA[,1],y1=newco.MFA[,2],col="lightgrey",length=0.1)#
#draw.circle(0,0,1)#
#draw.circle(0,0,0.5,lty="dashed")#
#text(newco.MFA[,1], newco.MFA[,2],rownames(newco.MFA),col=rep(my.palette(length(bloc)),bloc),cex=0.6) # 4 = nombre de martrices dans le tableau DATA#

# Between class analysis (analyse supervisée - FACTEUR 1 SITES#
BCA=bca(MFA,facteur1,scannf=F,nf=2)#
BCA$ratio#
varexp2=BCA$eig*100/sum(BCA$eig)#
plot(BCA$ls[,2]~BCA$ls[,1],pch=21,cex=2,col="white",bg=rainbow(length(levels(facteur1)))[facteur1],xlab=paste("Axis 1 : ",round(varexp2[1],2),"%"),ylab=paste("Axis 2 : ",round(varexp2[2],2),"%"),main="BCA scores")#
legend("topright",c("instrumental variable = Sites","Total inertia explained = 36%"),cex=0.8,box.lty=0)
ordihull(BCA$ls,facteur1,lab=T)#
plot(BCA$co[,2]~BCA$co[,1],type="n",xlab=paste("Axis 1 : ",round(varexp2[1],2),"%"),ylab=paste("Axis 2 : ",round(varexp2[2],2),"%"),ylim=c(-1,1),xlim=c(-1,1),main="BCA loadings")#
BCA$co$col=rep(rainbow(5),bloc)
arrows(x0=0,y0=0,x1=BCA$co[,1],y1=BCA$co[,2],col="lightgrey",length=0.1)#
cos2 <- as.matrix(BCA$co[,1:2])*as.matrix(BCA$co[,1:2]) # calcul des cos2 pour filtrer les variables
newco <- BCA$co[cos2[,1]>0.5 | cos2[,2]>0.5,] # utiliser une valeur de cos2 pour filtrer (ici 0.5)
oldco <- BCA$co[cos2[,1]<0.5 & cos2[,2]<0.5,]
text(newco[,1], newco[,2],rownames(newco),col=newco$col,cex=1.7)
text(oldco[,1], oldco[,2],rownames(oldco),col=alpha(oldco$col,0.3),cex=1.5) # la valeur alpha contrôle la transparence des variables < au seuil de cos2
summary(BCA) ### AJOUT
#draw.circle(0,0,1)#
#draw.circle(0,0,0.3,lty="dashed")#
abline(h=0,lty="dashed")
abline(v=0,lty="dashed")
 # 4 = nombre de martrices dans le tableau DATA#
legend("topleft",c("Physicochemical parameters", "Contaminants", "Prokaryotic biomass","Organic matter", "Pigments"),text.col=my.palette(length(bloc)),cex=0.8,box.lty=0)

################
# Test BCA
randtest(BCA)
###############


############ MFA as weighted PCA

par(layout(matrix(c(1,2,3,3),nrow=2)),mar=c(4,4,4,4))
MFA2=dudi.pca(pca.data,col.w=rep(c(1/vp.enviro[1],1/vp.conta[1],1/vp.bact[1],1/vp.CHN[1],1/vp.pig[1]),bloc),scannf=F,nf=2)
varexp1=MFA2$eig*100/sum(MFA2$eig)#
facteur1=factor(paste(DATA$site))#
facteur2=factor(paste(DATA$site,DATA$horizon))#
facteur3=factor(paste(DATA$horizon))#
plot(MFA2$li[,2]~MFA2$li[,1],pch=21,cex=3,col="white",bg=rainbow(length(levels(facteur1)))[facteur1],xlab=paste("Axis 1 : ",round(varexp1[1],2),"%"),ylab=paste("Axis 2 : ",round(varexp1[2],2),"%"),main="MFA scores")#
ordihull(MFA2$li,facteur2,lab=T)#

# Between class analysis (analyse supervisée - FACTEUR 1 SITES#
BCA2=bca(MFA2,facteur1,scannf=F,nf=2)#
BCA2$ratio#
varexp2=BCA$eig*100/sum(BCA$eig)#
plot(BCA2$ls[,2]~BCA2$ls[,1],pch=21,cex=3,col="white",bg=rainbow(length(levels(facteur1)))[facteur1],xlab=paste("Axis 1 : ",round(varexp2[1],2),"%"),ylab=paste("Axis 2 : ",round(varexp2[2],2),"%"),main="BCA scores")#
legend("topleft",c("instrumental variable = Sites","Total inertia explained = 39.8%"),cex=1.2,box.lty=0)
ordihull(BCA2$ls,facteur1,lab=T)#
plot(BCA2$co[,2]~BCA2$co[,1],type="n",xlab=paste("Axis 1 : ",round(varexp2[1],2),"%"),ylab=paste("Axis 2 : ",round(varexp2[2],2),"%"),ylim=c(-1,1),xlim=c(-1,1),main="BCA loadings")#
BCA2$co$col=rep(my.palette(5),bloc)
arrows(x0=0,y0=0,x1=BCA2$co[,1],y1=BCA2$co[,2],col="lightgrey",length=0.1)#
cos2 <- as.matrix(BCA2$co[,1:2])*as.matrix(BCA2$co[,1:2]) # calcul des cos2 pour filtrer les variables
newco <- BCA2$co[cos2[,1]>0.4 | cos2[,2]>0.4,] # utiliser une valeur de cos2 pour filtrer (ici 0.5)
oldco <- BCA2$co[cos2[,1]<0.4 & cos2[,2]<0.4,]
text(newco[,1], newco[,2],rownames(newco),col=newco$col,cex=1.7)
text(oldco[,1], oldco[,2],rownames(oldco),col=alpha(oldco$col,0.3),cex=1.5) # la valeur alpha contrôle la transparence des variables < au seuil de cos2
summary(BCA2) ### AJOUT
#draw.circle(0,0,1)#
#draw.circle(0,0,0.3,lty="dashed")#~/Desktop/plot-old-script.pdf
abline(h=0,lty="dashed")
abline(v=0,lty="dashed")
 # 4 = nombre de martrices dans le tableau DATA#
legend("topleft",c("Physicochemical parameters", "Contaminants", "Prokaryotic biomass","Organic matter", "Pigments"),text.col=my.palette(length(bloc)),cex=1.4,box.lty=0)






# Between class analysis (analyse supervisée - FACTEUR 2 SITES ET HORIZONS#
BCA=bca(MFA,facteur2,scannf=F,nf=2)#
BCA$ratio#
varexp2=BCA$eig*100/sum(BCA$eig)#
plot(BCA$ls[,2]~BCA$ls[,1],pch=21,cex=2,col="white",bg=rainbow(length(levels(facteur1)))[facteur1],xlab=paste("Axis 1 : ",round(varexp2[1],2),"%"),ylab=paste("Axis 2 : ",round(varexp2[2],2),"%"),main="BCA scores")#
legend("topleft",c("Interaction sites:sediment depth","Total inertia explained = 60%"),cex=0.8,box.lty=0)
ordihull(BCA$ls,facteur2,lab=T)#
plot(BCA$co[,2]~BCA$co[,1],type="n",xlab=paste("Axis 1 : ",round(varexp2[1],2),"%"),ylab=paste("Axis 2 : ",round(varexp2[2],2),"%"),ylim=c(-1,1),xlim=c(-1,1),main="BCA loadings")#
BCA$co$col=rep(my.palette(5),bloc)#
newco=BCA$co[abs(BCA$co[,1])>0 | abs(BCA$co[,2])>0,]#
arrows(x0=0,y0=0,x1=newco[,1],y1=newco[,2],col="lightgrey",length=0.1)#
text(newco[,1], newco[,2],rownames(newco),col= rep(my.palette(5),bloc),cex=0.8)
summary(BCA) ### AJOUT
#draw.circle(0,0,1)#
#draw.circle(0,0,0.3,lty="dashed")#
abline(h=0,lty="dashed")
abline(v=0,lty="dashed")
# 4 = nombre de martrices dans le tableau DATA#
legend("topleft",c("Physicochemical parameters", "Contaminants", "Prokaryotic biomass","Organic matter", "Pigments"),text.col=my.palette(length(bloc)),cex=0.8,box.lty=0)


# Between class analysis (analyse supervisée - FACTEUR 3 HORIZON#
BCA=bca(MFA,facteur3,scannf=F,nf=2)#
BCA$ratio#
varexp2=BCA$eig*100/sum(BCA$eig)#
par(mfrow=c(1,2))#
plot(BCA$ls[,2]~BCA$ls[,1],pch=21,cex=2,col="white",bg=rainbow(length(levels(facteur3)))[facteur3],xlab=paste("Axis 1 : ",round(varexp2[1],2),"%"),ylab=paste("Axis 2 : ",round(varexp2[2],2),"%"))#
ordihull(BCA$ls,facteur3,lab=T)#
plot(BCA$co[,2]~BCA$co[,1],type="n",xlab=paste("Axis 1 : ",round(varexp2[1],2),"%"),ylab=paste("Axis 2 : ",round(varexp2[2],2),"%"),ylim=c(-1,1),xlim=c(-1,1))#
BCA$co$col=rep(my.palette(2),bloc)#
newco=BCA$co[abs(BCA$co[,1])>0 | abs(BCA$co[,2])>0,]#
text(newco[,1], newco[,2],rownames(newco),col= rep(my.palette(2),bloc))
arrows(x0=0,y0=0,x1=newco[,1],y1=newco[,2],col="lightgrey",length=0.1)#
text(newco[,1], newco[,2],rownames(newco),col= rep(my.palette(5),bloc))## AJOUT
summary(BCA) ### AJOUT
draw.circle(0,0,1)#
draw.circle(0,0,0.3,lty="dashed")#
legend("topleft",c("Physicochemical parameters", "Contaminants", "Prokaryotic biomass","Organic matter", "Pigments"),text.col=my.palette(length(bloc)),cex=0.8,box.lty=0)





