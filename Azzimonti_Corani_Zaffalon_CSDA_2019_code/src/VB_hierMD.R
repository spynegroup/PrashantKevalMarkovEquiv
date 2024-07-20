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
source("src/opt.R")

VB_hierMD=function(counts, kappa, tau, nu, x_dim, y_dim, alpha0, s, param_VI = NULL, print_step = FALSE){

    if(is.null(param_VI)){
        KAPPA_NEWTON_THRESH = 1e-6
        KAPPA_MAX_ITER = 5000
        LIKELIHOOD_DECREASE_ALLOWANCE = 1e-5
        TAU_NEWTON_THRESH = 1e-6
        TAU_MAX_ITER = 5000
        KAPPA_TAU_THRESH = 1e-6
        KAPPA_TAU_MAX_ITER = 5000
        EM_ITER_MAX = 1000
        EM_THRESH = 1e-6
    } else {
        KAPPA_NEWTON_THRESH = param_VI$KAPPA_NEWTON_THRESH
        KAPPA_MAX_ITER = param_VI$KAPPA_MAX_ITER
        LIKELIHOOD_DECREASE_ALLOWANCE = param_VI$LIKELIHOOD_DECREASE_ALLOWANCE
        TAU_NEWTON_THRESH = param_VI$TAU_NEWTON_THRESH
        TAU_MAX_ITER = param_VI$TAU_MAX_ITER
        KAPPA_TAU_THRESH = param_VI$KAPPA_TAU_THRESH
        KAPPA_TAU_MAX_ITER = param_VI$KAPPA_TAU_MAX_ITER
        EM_ITER_MAX = param_VI$EM_ITER_MAX
        EM_THRESH = param_VI$EM_THRESH
    }

    iter_em = 0
    whole_converged = 1
    whole_likelihood_old = 0
    whole_likelihood = 0
    kappa_tau_like = 0
    nu_like = 0
    whole_likelihood_vec=NULL

    while(iter_em < EM_ITER_MAX & abs(whole_converged) > EM_THRESH){
        iter_em = iter_em + 1

        OPT_KAPPA_TAU = opt_kappa_tau(nu, kappa, x_dim, y_dim, alpha0, s,
                                      KAPPA_TAU_MAX_ITER, KAPPA_TAU_THRESH, KAPPA_MAX_ITER, KAPPA_NEWTON_THRESH,
                                      TAU_MAX_ITER, TAU_NEWTON_THRESH)
        tau = OPT_KAPPA_TAU$tau
        kappa = OPT_KAPPA_TAU$kappa
        kappa_tau_like = OPT_KAPPA_TAU$likelihood

        OPT_NU = opt_nu(counts,kappa, x_dim, y_dim, alpha0, s)
        nu = OPT_NU$nu
        nu_like = OPT_NU$likelihood

        common=0
        for(i in 1:x_dim){
          common=common-lgamma(alpha0[i])
        }
        whole_likelihood = kappa_tau_like + nu_like - common + lgamma(sum(alpha0))
        whole_likelihood_vec=c(whole_likelihood_vec,whole_likelihood)

        whole_converged = (whole_likelihood - whole_likelihood_old) / abs(whole_likelihood_old);
        whole_likelihood_old = whole_likelihood

        if (print_step) {
            print(iter_em)
            print(paste("converged: ", whole_converged))
            print(paste("value: ", whole_likelihood))
        }
    }

    return(list(nu = nu, tau = tau, kappa = kappa,
                whole_likelihood = whole_likelihood, whole_converged = whole_converged,
                kappa_tau_like = kappa_tau_like, nu_like = nu_like, whole_likelihood_vec = whole_likelihood_vec))
}



