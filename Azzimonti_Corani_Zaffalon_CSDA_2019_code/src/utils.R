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

log_sum=function(log_a, log_b){
  if (log_a < log_b){
    v = log_b+log(1 + exp(log_a-log_b));
  }else{
    v = log_a+log(1 + exp(log_b-log_a));
  }
  return(v);
}

digamma=function(x){
  x=x+6;
  p=1/(x*x);
  p=(((0.004166666666667*p-0.003968253986254)*p+
        0.008333333333333)*p-0.083333333333333)*p;
  p=p+log(x)-0.5/x-1/(x-1)-1/(x-2)-1/(x-3)-1/(x-4)-1/(x-5)-1/(x-6);
  return(p);
}

trigamma=function(x){
  x=x+6;
  p=1/(x*x);
  p=(((((0.075757575757576*p-0.033333333333333)*p+0.0238095238095238)
       *p-0.033333333333333)*p+0.166666666666667)*p+1)/x+0.5*p;
  for (i in 1:6){
    x=x-1;
    p=1/(x*x)+p;
  }
  return(p);
}


tetragamma=function(x){
  x=x+6;
  p=1/(x*x);
  p=(((((0.3 - 0.833333333333333 * p) * p - 0.166666666666666) * p + 0.166666666666666) * p - 0.5) * p - 1/x - 1) * p;
  for (i in 1:6){
    x=x-1;
    p = p - 2 / (x*x*x);
  }
  return(p);
}

import.package = function(pkg) {
  if (!requireNamespace(pkg))
    stop("this function requires the ", pkg, " package.")
}
