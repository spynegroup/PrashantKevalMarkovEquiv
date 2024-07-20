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

VB_stan_hierMD<- function(xValues, yValues,xStates=2,yStates=2) {

  import.package("rstan")
  import.package("MCMCpack")

  trainingSamples=length(xValues)
  #we aim at estimating the probability of X|Y
  #yStates is the number of states of the conditioning variable Y

  rstan_options(auto_write = TRUE)
  options(mc.cores = parallel::detectCores())

  yCounts <- as.numeric(table(factor(yValues, levels = 1:yStates)))

  xyCounts <- matrix (nrow = xStates, ncol = yStates)
  for (i in 1:yStates){
    xyCounts[,i] <- as.numeric(table(factor(xValues[yValues==i], levels = 1:xStates)))
  }

  #prepare data for Stan
  dataList = list(
    yStates = yStates,
    xStates = xStates,
    Ntotal = trainingSamples,
    yValues = yValues,
    xValues = xValues,
    xyCounts = xyCounts,
    alpha0 = rep(1,xStates),
    s =  xStates
  )


  #shrinkage estimation, we need two additional parameters:
  #alpha0, the vector of pars of the highest-level prior
  #s, the amount of smoothing performed locally
  stanModel <-  stan_model(model_code="data {
      int<lower=2> yStates;
      int<lower=2> xStates;
      int Ntotal; //number of instances
      vector[yStates] yCounts;
      matrix[xStates,yStates] xyCounts;
      vector[xStates] alpha0; //the highest-level Dir coefficients
      real<lower=0> s; //amount of local smoothing. No prior on s for the moment, to be added later however.
}


parameters {
simplex[yStates] thetaY;
simplex[xStates] thetaX[yStates];
simplex[xStates] alpha;//the prior  for the local states
}

model {

//sample the thetaY from the posterior  Dirichlet
thetaY ~ dirichlet (1.0 + yCounts);

//we treat alpha0 as a fixed vector of ones
alpha   ~ dirichlet (alpha0);

for (y in 1:yStates){
//sample the thetaXY from a Dir(1,1,1,...1)
thetaX[y] ~ dirichlet (s*alpha + col(xyCounts,y));
}
}")
    stanVb <-  vb(stanModel, data =  dataList)
    stanResultsShrinkage<- extract(stanVb, permuted = TRUE)

  postSamples <- dim(stanResultsShrinkage$thetaY)[1]


  #pre-allocate
  currentThetaY <- vector(length = yStates)
  currentThetaXgivenY <- array(dim=c(xStates,yStates))

  for (currentSample in 1:postSamples){
    currentThetaY <- stanResultsShrinkage$thetaY[currentSample,]
    for (currentYstate in 1:yStates){
      currentThetaXgivenY[,currentYstate] <- stanResultsShrinkage$thetaX[currentSample,currentYstate,]
    }
  }

   postThetaXgivenY=NULL
   for(k in 1:dim(stanResultsShrinkage$thetaX)[2]){
   postThetaXgivenY=cbind(postThetaXgivenY,stanResultsShrinkage$thetaX[,k,])
   }

  return(list(theta=postThetaXgivenY, thetaY=stanResultsShrinkage$thetaY,alpha=stanResultsShrinkage$alpha))
}

