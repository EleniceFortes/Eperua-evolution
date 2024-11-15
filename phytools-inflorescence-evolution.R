#######################################################################################
### Discrete character evolution with SIMMAP   #######################################
######################################################################################

#load libraries
library(ape)
library(phytools)
library(treeio)
library(ggtree)
library(ggplot2)
library(BioGeoBEARS)

##################################################

trfn = np("Beast_concat_30exons_27tips_renamed_MCCtree.nwk")
moref(trfn)

eperua.tree = read.tree(trfn)
eperua.tree
dev.off() #in case it is not showing
plot(eperua.tree)
title("Eperua phylogeny")
axisPhylo() # plots timescale

#read data
eperua.trait<-read.csv("matrix_traits.csv", row.names = 1, stringsAsFactors = T)
eperua.trait

#checking names tree x data 
library(ape)
library(phangorn)
library(phytools)
library(geiger)
name.check(eperua.trait,eperua.tree)

##extract Eperua traits as a vector
floral.type<-setNames(eperua.trait[,1], rownames(eperua.trait))
floral.type

##set colors for plotting
cols<-setNames(c("grey","purple"),levels(floral.type))
cols
##plot tree & data
plotTree.datamatrix(eperua.tree,as.data.frame(floral.type),
                    colors=list(cols),header=F,fsize=0.45)
##add legend
legend("topleft",legend=levels(floral.type), pch=22,
       pt.cex = 1.5,pt.bg = cols,bty = "n",cex = 0.8)

#testing models
##fit ER model
fitER<-fitMk(eperua.tree, floral.type, model="ER")
##fit ARD model
fitARD<-fitMk(eperua.tree, floral.type, model="ARD")
##fit  ARD model
fitSYM<-fitMk(eperua.tree, floral.type, model="SYM")
##extract AIC values for each model
aic<-c(AIC(fitER), AIC(fitARD), AIC(fitSYM))
##print summary table
data.frame(model=c("ER","ARD","SYM"),logL=c( 
      logLik(fitER), logLik(fitARD), logLik(fitSYM)), AIC=aic,delta.AIC=aic-min(aic))

############################################################################################
#generate 1,00 stochastic maps in which the transiction rate rate 
#is sampled from its posterior distribuition

mtrees<-make.simmap(eperua.tree,floral.type,model="ER",nsim=1000,
                    Q="mcmc",vQ=0.01,prior=list(use.empirical=T),samplefreq=10)
mtrees


#create a 10x10 grid of plot cells
par(mfrow=c(10,10))
##graph 100 stochastic map trees, sampled evenly from our set of 1000
null<-sapply(mtrees[seq(10,1000,by=10)],plot,colors=cols,lwd=1,ftype="off")

#compute posterior probabilities at nodes
pd<-summary(mtrees)
pd
##create a plot showing PP at all nodes of the tree
#tiff device
tiff("plot_inflorescence_cladogram_ER.tiff", width = 10, height = 10, units = 'in', res = 600)
#plotting
par(mfrow=c(1,1))
plot(pd,colors=cols,fsize=1.2,ftype="i",lwd=2,
     offset=0.4,ylim=c(-1,Ntip(eperua.tree)),cex=c(0.7,0.45))
legend("topleft",legend=c("Long and pendulous", "Short and erect"),pch=22,pt.cex=4.5,
       pt.bg=cols,bty="n",cex=1.5)
axisPhylo(1,pos=0.5)
dev.off()

