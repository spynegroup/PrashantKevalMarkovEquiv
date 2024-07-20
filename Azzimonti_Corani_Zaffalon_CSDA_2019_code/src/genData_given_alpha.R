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

genData_given_alpha <- function(trainingSamples=10,xStates=4,yStates=2, alphaXY) {
  library(MCMCpack)  
  
  xyStates=xStates*yStates
  thetaXY <- rdirichlet(1,rep(alphaXY*xyStates,yStates))
  thetaY=colSums(matrix(thetaXY,xStates,yStates))
  thetaXgivenY=t(t(matrix(thetaXY,xStates,yStates))/thetaY)

  xyValues = sample (1:xyStates,trainingSamples,replace = TRUE, prob = thetaXY)
  yValues=((xyValues-1) %/% xStates) +1
  xValues=xyValues- ((xyValues-1) %/% xStates) *xStates 
  xyCounts <- matrix(as.numeric(table(factor(xyValues, levels = 1:xyStates))),xStates,yStates)
  thetaXY_mat=matrix(thetaXY,xStates,yStates)
  return(list(counts=xyCounts, xyValues=xyValues,xValues=xValues,yValues=yValues,thetaXgivenY=thetaXgivenY,thetaY=thetaY,thetaXY=thetaXY_mat))
}

