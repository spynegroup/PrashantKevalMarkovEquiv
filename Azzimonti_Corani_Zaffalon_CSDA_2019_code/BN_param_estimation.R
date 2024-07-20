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
source("src/bn.fit.hier.R")

library(bnlearn)
data=bnlearn::coronary
features=names(data)
class="Family"
features_evidence=names(data)[features!=class]

#Single domain hierachical parameter estimation
dag = tree.bayes(data, class, features_evidence)
plot(dag)
fit_hier=bn.fit.hier(dag, data)

#fit_hier=bn.fit.hier(dag, data)

group="Smoking"
features_evidence_g=features_evidence[features_evidence!=group]

# Multi domain hierarchical parameter estimation
dag_multi = tree.bayes(data, class, features_evidence_g)
plot(dag_multi)
fit_hier=bn.fit.hier(dag_multi, data,method="multi-domain", group=group)
