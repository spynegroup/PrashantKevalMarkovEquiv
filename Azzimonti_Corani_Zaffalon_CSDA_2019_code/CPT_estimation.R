###############################################################################
###    IPG - IDSIA
###    Authors: L. Azzimonti, G. Corani
###    Reference: "Hierarchical estimation of parameters in Bayesian networks",
###    Date: May 16, 2019
###
###    This code is for demostration purpose only
###    for licensing contact laura@idsia.ch
###    COPYRIGHT (C) 2019 IDSIA
###############################################################################
source("src/genData_given_alpha.R")
source("src/VB_hierMD.R")

#library(entropy)
library(MCMCpack)

I=2
#|X|, i.e., number of states of the child variable
J=4
#|Y|, i.e.,number of states of the child variable
n=80
#number of observations

alpha <- rdirichlet(1,rep(1,I))
#true underlying parameter vector alpha


data_all=genData_given_alpha(n, I, J,alpha)
counts=data_all$counts
data=data_all	
data$xValues=data_all$xValues
data$yValues=data_all$yValues
dataValues=data.frame(xValues = factor(data_all$xValues,levels=1:I), yValues = factor(data_all$yValues,levels=1:J))
theta_true=as.numeric(data_all$thetaXgivenY)
theta_freq=as.vector(t(t(counts)/colSums(counts)))
theta_bayes1=as.vector(t(t(counts+1/(I))/(colSums(counts)+1)))
theta_bayes10=as.vector(t(t(counts+10/(I))/(colSums(counts)+10)))
theta_bayesI=as.vector(t(t(counts+I/(I))/(colSums(counts)+I)))

#variational inference
dirichlet_prior=rep(1,I) #alpha0
s=I #equivalent sample size

nu=matrix(1,I,J)/I
kappa=rep(1,I)/I
tau=100
VB_ris=VB_hierMD(data$counts,kappa, tau, nu, I, J, dirichlet_prior, s)

alpha_VB=VB_ris$kappa
theta_VB=as.numeric(t(t(VB_ris$nu)/colSums(VB_ris$nu)))

err_theta=cbind(theta_freq-theta_true,theta_bayes1-theta_true,theta_bayes10-theta_true,theta_bayesI-theta_true,theta_VB-theta_true)

MSE_theta=c(colMeans(err_theta^2), n, I, J)
names(MSE_theta)=c("ML","BAYES1","BAYES10","BAYESI", "VB", "N", "I", "J")
print(MSE_theta)
