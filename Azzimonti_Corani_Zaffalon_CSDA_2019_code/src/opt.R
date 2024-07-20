###############################################################################
###    IPG - IDSIA
###    Authors: L. Azzimonti, G. Corani, S. Toniolo
###    Reference: "Hierarchical estimation of parameters in Bayesian networks",
###    Date: May 16, 2019
###
###    This code is for demostration purpose only
###    for licensing contact laura@idsia.ch
###    COPYRIGHT (C) 2019 IDSIA
###############################################################################
source("src/utils.R")

opt_nu=function(counts, kappa, x_dim, y_dim, alpha0, s ){
  nu=matrix(0, x_dim, y_dim)
  likelihood=0

  for (i in 1:x_dim) {
    for(j in 1:y_dim){
      nu[i,j] = s*kappa[i]+counts[i,j];
    }
  }
  sum_x_nu=colSums(nu)
  common=0
  common2=0
  common3=0
  common4=0
  common5=0
  fact=0
  for(j in 1:y_dim){
    for(i in 1:x_dim){
      common=common+lgamma(nu[i,j])# terzo termine della loglike
      common=common-(nu[i,j]-1)*(digamma(nu[i,j])-digamma(sum_x_nu[j]))#quarto termine della loglike - cambiato segno
      common3=common3- digamma(nu[i,j])# parte del secondo termine della loglike
      common2=common2+counts[i,j]*(digamma(nu[i,j])-digamma(sum_x_nu[j]))#primo termine loglike
      fact=fact+lfactorial(counts[i,j])
    }
    common4=common4 +lgamma(sum_x_nu[j])#quinto termine della loglike
    common5=common5+ (s-x_dim)*digamma(sum_x_nu[j])# seconda parte del secondo termine della loglike (manca ancora un pezzo)
  }
  likelihood=likelihood-common4+common+common3+common2+common5 # cambiato segno di common5
  #likelihood=likelihood+y_dim*lfactorial(sum(counts))-fact # non so da dove derivi - eliminata questa riga


  return(list(nu=nu, likelihood=likelihood))
}



opt_kappa_tau=function(nu, kappa, x_dim, y_dim, alpha0, s, KAPPA_TAU_MAX_ITER, KAPPA_TAU_THRESH, KAPPA_MAX_ITER, KAPPA_NEWTON_THRESH, TAU_MAX_ITER, TAU_NEWTON_THRESH){
  dep_likelihood = 0;
  indep_likelihood = 0;

  kappa_tau_likelihood = 0;
  kappa_tau_likelihood_old = 0;
  kappa_tau_converged = 1;
  kappa_tau_loop = 0;

  while((kappa_tau_loop < 2) | ((abs(kappa_tau_converged) > KAPPA_TAU_THRESH) & (kappa_tau_loop < KAPPA_TAU_MAX_ITER))) {
    kappa_tau_loop =kappa_tau_loop + 1;

    kappa_tau_likelihood = 0;
    OPT_TAU=opt_tau(kappa, x_dim, y_dim, alpha0, s, TAU_MAX_ITER, TAU_NEWTON_THRESH)
    tau=OPT_TAU$tau
    tau_like=OPT_TAU$likelihood

    OPT_KAPPA=opt_kappa(nu, tau, x_dim, y_dim, alpha0, s, KAPPA_MAX_ITER, KAPPA_NEWTON_THRESH)
    kappa=OPT_KAPPA$kappa
    kappa_like=OPT_KAPPA$likelihood

    kappa_tau_likelihood = kappa_tau_likelihood+ kappa_like

    for (i in 1:x_dim) {
      kappa_tau_likelihood = kappa_tau_likelihood +digamma(tau) * (-alpha0[i] + tau * kappa[i] + (s * kappa[i] - 1) * y_dim);
    }
    kappa_tau_likelihood =kappa_tau_likelihood - lgamma(tau);
    kappa_tau_likelihood =kappa_tau_likelihood - s * y_dim * (x_dim - 1) / tau;

    kappa_tau_converged = (kappa_tau_likelihood - kappa_tau_likelihood_old) / abs(kappa_tau_likelihood_old);
    kappa_tau_likelihood_old = kappa_tau_likelihood;
  }

  kappa
  tau

  #alpha
  precompute = 0;
  alpha_likelihood = 0;
  precompute = y_dim * (digamma(tau) - (x_dim - 1) / tau);
  digamma_sum_over_children=rep(0,x_dim)

  for (i in 1:x_dim) {
    alphakappai = s * kappa[i];
    for(j in 1:y_dim){
      digamma_sum_over_children[i]=digamma_sum_over_children[i]+ digamma(nu[i,j])
    }
    alpha_likelihood =alpha_likelihood - lgamma(alphakappai);
    precompute = precompute +y_dim * kappa[i] * (log(kappa[i]) - digamma(tau * kappa[i]));
    precompute = precompute +kappa[i] * digamma_sum_over_children[i];
  }

  alpha_likelihood = y_dim * (alpha_likelihood + lgamma(s)) + s * precompute;
  indep_likelihood = indep_likelihood + alpha_likelihood;

  digamma_tau = digamma(tau);

  indep_likelihood = indep_likelihood - y_dim * x_dim * digamma_tau;

  for (i in 1:x_dim) {
    kappai = kappa[i];
    taukappai = tau * kappai;
    digammataukappai = digamma(taukappai);
    common = (digammataukappai - digamma_tau);

    indep_likelihood = indep_likelihood -y_dim * (log(kappai) - digammataukappai);
    indep_likelihood = indep_likelihood +lgamma(taukappai);
    indep_likelihood = indep_likelihood - taukappai * common;
    dep_likelihood = dep_likelihood +alpha0[i] * common;
  }

  indep_likelihood =indep_likelihood - lgamma(tau);
  likelihood= indep_likelihood+dep_likelihood
  return(list(tau=tau,kappa=kappa, likelihood=likelihood, indep_likelihood=indep_likelihood, dep_likelihood=dep_likelihood))
}


opt_kappa=function(nu, tau, x_dim, y_dim, alpha0, s, KAPPA_MAX_ITER, KAPPA_NEWTON_THRESH){
  dep_new_likelihood=0
  indep_new_likelihood=0

  kappa=rep(1,x_dim)/x_dim #beta
  digamma_sum_over_children=rep(0,x_dim)
  for (i in 1:x_dim) {
    taukappai = tau * kappa[i];
    alphakappai = s * kappa[i];
    common = alpha0[i] + y_dim * (1 - alphakappai) - taukappai;
    digammataukappai = digamma(taukappai);
    logkappai = log(kappa[i]);
    for(j in 1:y_dim){
      digamma_sum_over_children[i]=digamma_sum_over_children[i]+ digamma(nu[i,j])
    }
    dep_new_likelihood = dep_new_likelihood + digammataukappai * common;
    indep_new_likelihood = indep_new_likelihood -y_dim * (lgamma(alphakappai) + (1 - alphakappai) * logkappai);
    indep_new_likelihood = indep_new_likelihood +alphakappai * digamma_sum_over_children[i];
    dep_new_likelihood = dep_new_likelihood +lgamma(taukappai);
  }
  new_likelihood = indep_new_likelihood + dep_new_likelihood;
  initial_likelihood = new_likelihood;

  sqr_newton_decrement=10* KAPPA_NEWTON_THRESH
  step_size=1
  iter=0

  #qua inizia il ciclo while
  while(iter < KAPPA_MAX_ITER & sqr_newton_decrement > KAPPA_NEWTON_THRESH * 2 & step_size > 1e-8 ){
    iter=iter+1
    #print(iter)
    g=rep(0,x_dim)
    h=rep(0,x_dim)
    invhsum=0
    goverhsum=0
    coefficient = 0

    for(i in 1:x_dim){
      taukappai = tau * kappa[i];
      alphakappai = s * kappa[i];
      common = alpha0[i] + y_dim * (1 - alphakappai) - taukappai;
      digammataukappai = digamma(taukappai);
      trigammataukappai = trigamma(taukappai);
      logkappai = log(kappa[i]);
      #digamma_sum_over_children=y_dim * digamma(nu[i])

      g[i] = tau * trigammataukappai * common +
        - y_dim * s * (digamma(alphakappai) - logkappai + digammataukappai - 1) +
        - y_dim / kappa[i]+ s * digamma_sum_over_children[i];

      h[i] = tau * tau * tetragamma(taukappai) * common +
        - tau * trigammataukappai * (tau + 2 * s * y_dim) +
        - s * s * trigamma(alphakappai) * y_dim +
        + s * y_dim / kappa[i] +
        + y_dim / (kappa[i] * kappa[i]);

      invhsum = invhsum + 1 / h[i];
      goverhsum = goverhsum + g[i] / h[i];
    }
    old_likelihood = new_likelihood;

    coefficient = goverhsum / invhsum;
    sqr_newton_decrement = 0;
    expected_increase = 0;
    step_size = 1;

    delta_kappa=rep(0, x_dim)

    for (i in 1:x_dim) {
      delta_kappa[i] = (coefficient - g[i]) / h[i];
      sqr_newton_decrement = sqr_newton_decrement - h[i] * delta_kappa[i] * delta_kappa[i]; # this one is maximization
      expected_increase =expected_increase + g[i] * delta_kappa[i];
      if (delta_kappa[i] < 0) {
        limit = (kappa[i] - 1e-10) / -(delta_kappa[i]);
        if (step_size > limit) {
          step_size = limit;
        }
      }
    }

    new_kappa=kappa
    while(new_likelihood < old_likelihood + 0.4 * step_size * expected_increase & step_size > 1e-8){
      indep_new_likelihood = 0.0;
      dep_new_likelihood = 0.0;
      for (i in 1:x_dim) {
        new_kappa[i] = kappa[i] + step_size * delta_kappa[i];

        taukappai = tau * new_kappa[i];
        alphakappai = s * new_kappa[i];
        common = alpha0[i] + y_dim * (1 - alphakappai) - taukappai;# DA VERIFICARE
        logkappai = log(new_kappa[i]);
        #digamma_sum_over_children=y_dim * digamma(nu[i])
        dep_new_likelihood = dep_new_likelihood +digamma(taukappai) * common;
        indep_new_likelihood = indep_new_likelihood -y_dim * (lgamma(alphakappai) + (1 - alphakappai) * logkappai);
        indep_new_likelihood = indep_new_likelihood +alphakappai * digamma_sum_over_children[i];
        dep_new_likelihood = dep_new_likelihood +lgamma(taukappai);
      }

      new_likelihood = indep_new_likelihood + dep_new_likelihood;
      step_size = step_size * 0.9;
    }
    kappa=new_kappa
  }

  return(list(kappa=kappa, likelihood=new_likelihood))
}


opt_tau=function(kappa, x_dim, y_dim, alpha0, s, TAU_MAX_ITER, TAU_NEWTON_THRESH){

  iter = 0;
  d1 = 1;
  d2 = 0;

  init_tau = 100;
  likelihood = 0;
  old_likelihood = 0;

  tau = init_tau;
  log_tau = log(tau);

  while (iter < TAU_MAX_ITER & abs(d1) > TAU_NEWTON_THRESH){
    iter=iter+1
    #print(iter)

    d1 = 0;
    d2 = 0;
    likelihood = 0;

    common2 = y_dim * s * (x_dim - 1) / tau;
    trigammatau = trigamma(tau);
    digammatau = digamma(tau);
    for (i in 1:x_dim) {
      taukappai = tau * kappa[i];
      trigammataukappai = trigamma(taukappai);
      common = alpha0[i] - taukappai + y_dim * (1 - s * kappa[i]);

      d1 =d1 + (trigammataukappai * kappa[i] - trigammatau) * common;
      d2 =d2 + kappa[i] * kappa[i] * ( tetragamma(taukappai) * common - trigammataukappai);
      d2 =d2 - tetragamma(tau) * common;
      likelihood = likelihood +(digamma(taukappai) - digammatau) * common + lgamma(taukappai);
    }
    d1 =d1 + common2 / tau;
    d2 =d2 + trigammatau - 2 * common2 / tau / tau;

    likelihood =likelihood - lgamma(tau);
    likelihood =likelihood - common2;

    old_likelihood = likelihood;

    log_tau = log_tau - d1 / (d2 * tau + d1);
    tau = exp(log_tau);
    if (is.nan(tau) | tau < 1e-10) {
      init_tau = init_tau * 10;
      tau = init_tau;
      log_tau = log(tau);
      old_likelihood = 0;
    }
  }

  return (list(tau=tau,likelihood=likelihood ));
}

