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
source("src/VB_hierMD.R")
source("src/convert.R")
source("src/predict.R")


bn.fit.hier = function(dag, data, method="CPT", group=NULL, alpha=1,laplace=TRUE) {
  stan=FALSE
  n = nrow(data)


  levels_old=list()
  levels_new=list()
  for(node in names(data)){
    levels_old[[node]]=levels(data[,node])
    levels_split=lapply(levels(data[,node]),strsplit,split="/")
    levels_new[[node]]=unlist(lapply(lapply(levels_split,unlist), paste0,collapse=""))
    levels(data[,node])=levels_new[[node]]
  }

  # Closed form solution for parameter estimation of Multinomial- Dirichlet model
  compute_theta_MD_mean=function(N,eta0){
    # Parameter learning by means of traditional Multinomial-Dirichlet model.
    # N = vector of counts
    # eta0 = dirichlet prior

    Ydim=length(N)
    theta_mean=(N+eta0/Ydim)/(sum(N)+eta0)
    return(theta_mean)
  }

  hier.bayes=function(dag,group){
    # Creates a new DAG with group as parent variable for every node (only for method="multi-domain")

    # dag = dag structure in bnlearn format
    # group = name of the variable associated to the domain (only for method="multi-domain")

    for(node in names(dag$nodes)){
      dag$nodes[[node]]$parents=c(dag$nodes[[node]]$parents,group)
      dag$nodes[[node]]$mb=c(dag$nodes[[node]]$mb,group)
      dag$nodes[[node]]$nbr=c(dag$nodes[[node]]$nbr,group)
    }
    nodes=names(dag$nodes)
    dag$nodes[[group]]$mb=nodes
    dag$nodes[[group]]$nbr=nodes
    dag$nodes[[group]]$parents=NULL
    dag$nodes[[group]]$childern=nodes

    dag$arcs=rbind(dag$arcs,cbind(rep("m",length(nodes)),nodes))
    return(dag)
  }

  # Define the fitting function.
  # Fitting function for CPT estimation by means of hierarchical Multinomial-Dirichlet model.
  # The hierarchial model induces a shrinkage effect between columns of the CPT.
  if(method=="CPT"){
    fit = function(node,alpha) {
      # Store the labels of the parents and the children
      parents = dag$nodes[[node]]$parents
      children = dag$nodes[[node]]$children

      if(length(dag$nodes[[node]]$parents)>0){
        # Hierarchical estimation of CPT when the current node has non empty parent set
        # Store in data all the useful information
        summary_data=NULL
        # Rows of the CPT correspond to the states of the current node
        summary_data$xValues=as.numeric(data[,node])
        levelsx=levels(data[,node])
        summary_data$xValues=as.numeric(data[,node])

        # Columns of the CPT correspond to the joint states of all the parents
        # levelsy contains all the possible joint states of the parents
        newstates=data[,dag$nodes[[node]]$parents[length(dag$nodes[[node]]$parents)]]
        levelsy=levels(data[,dag$nodes[[node]]$parents[length(dag$nodes[[node]]$parents)]])
        if(length(dag$nodes[[node]]$parents)>1){
          for(j in (length(dag$nodes[[node]]$parents)-1):1){
            newstates=paste0(newstates, "/", data[,dag$nodes[[node]]$parents[j]])
            levels_new=levels(data[,dag$nodes[[node]]$parents[j]])
            levelsy=paste0(rep(levelsy,each=length(levels_new)), "/",rep(levels_new,length(levelsy)))
          }
        }
        meaning_yStates=levelsy
        summary_data$yValues=as.numeric(factor(newstates,levels=levelsy))

        # Contingency table
        summary_data$counts=table(data[,node],factor(newstates,levels=levelsy) )

        # Dimension of the contingency table
        Xdim=length(levelsx)
        Ydim=length(levelsy)

        # Set parameters of the hierarchical model: alpha0
        dirichlet_prior=rep(1,Xdim)

        # Set the starting point for variational inference
        nu=matrix(1,Xdim,Ydim)/Xdim
        kappa=rep(1,Xdim)/Xdim
        tau=100

        if(laplace){
          s=alpha*Xdim
        }else{
          s=alpha/Ydim
        }

        # CPT estimation by means of variational inference
        #variational inference algorithm for hierarchical Multinomial Dirichlet model
        if(!stan){
          VB=VB_hierMD(summary_data$counts,kappa, tau, nu, Xdim, Ydim, dirichlet_prior, s)
          theta_VB=t(t(VB$nu)/colSums(VB$nu))
        }else{
          VB_stan=VB_stan_hierMD(summary_data$xValues,summary_data$yValues,Xdim, Ydim)
          theta_VB=matrix(colMeans(VB_stan$theta),Xdim,Ydim)
        }
        theta=theta_VB

        # Convert the matrix theta_VB into an array with dimension equal to the numer of parents + 1
        tab=convert2array(theta_VB,meaning_yStates,data,node,dag,Xdim,Ydim)
        # Set parameters of tab
        dim(tab)=dim(tab)
        dimnames(tab)[[1]]=levels_old[[node]]
        for(j in 1:(length(dag$nodes[[node]]$parents))){
          dimnames(tab)[[j+1]]=levels_old[[dag$nodes[[node]]$parents[j]]]
        }
        names(dimnames(tab))=c(node,dag$nodes[[node]]$parents)
      }else{
        # Traditional Bayesian estimation of CPT when the current node has empty parent set
        # Store in data all the useful information
        summary_data=NULL

        # Rows of the CPT correspond to the states of the current node
        summary_data$xValues=as.numeric(data[,node])
        levelsx=levels(data[,node])
        Xdim=length(levelsx)
        Ydim=1
        summary_data$xCounts=as.numeric(table(factor( summary_data$xValues,levels=1:Xdim)))

        if(laplace){
          s=alpha*Xdim
        }else{
          s=alpha/Ydim
        }
        # Closed form solution for parameter estimation of Multinomial- Dirichlet model
        theta=compute_theta_MD_mean( summary_data$xCounts,s)
        tab=theta

        # Set parameters of tab
        names(tab)=levels_old[[node]]
        dim(tab)=length(levelsx)
        dimnames(tab)[[1]]=levels_old[[node]]
        names(dimnames(tab))=c(node)
      }

      # Set class
      class(tab)="table"
      class = ifelse(is(data[, node], "ordered"), "bn.fit.onode", "bn.fit.dnode")

      structure(list(node = node, parents = parents, children = children,
                     prob = tab), class = class)
    }
  # End fit for CPT by means of hierarchical Multinomial Dirichlet model
  }else if(method=="multi-domain"){
    fit = function(node,alpha) {
      #print(node)
      # Store the labels of the parents and the children
      parents = dag$nodes[[node]]$parents
      children = dag$nodes[[node]]$children

      if(length(dag$nodes[[node]]$parents)>=0){
        # Hierarchical estimation of CPT when the current node has non empty parent set
        # Store in data all the useful information
        summary_data=NULL
        # Rows of the CPT correspond to the states of the current node
        #summary_data$xValues=as.numeric(data[,node])
        levelsx=levels(data[,node])
        #summary_data$xValues=as.numeric(data[,node])

        # Columns of the CPT correspond to the joint states of all the parents
        # levelsx contains all the possible joint states of the parents
        newstates=data[,node]
        if(length(dag$nodes[[node]]$parents)>1){
          for(j in (length(dag$nodes[[node]]$parents)-1):1){
            levels_new=levels(data[,dag$nodes[[node]]$parents[j]])
            levelsx=paste0(rep(levelsx,each=length(levels_new)), "/",rep(levels_new,length(levelsx)))
            newstates=paste0(newstates, "/", data[,dag$nodes[[node]]$parents[j]])

          }
        }
        meaning_xStates=levelsx
        summary_data$xValues=as.numeric(factor(newstates,levels=levelsx))
        summary_data$yValues=as.numeric(as.factor(data[,group]))

        # Contingency table
        summary_data$counts=table(factor(newstates,levels=levelsx),summary_data$yValues)

        # Dimension of the contingency table
        Xdim=length(levelsx)
        Ydim=length(levels(as.factor(data[,group])))

        # Set parameters of the hierarchical model: alpha0
        dirichlet_prior=rep(1,Xdim)

        # Set the starting point for variational inference
        nu=matrix(1,Xdim,Ydim)/Xdim
        kappa=rep(1,Xdim)/Xdim
        tau=100

        if(laplace){
          s=alpha*Xdim
        }else{
          s=alpha/Ydim
        }

        # CPT estimation by means of variational inference
        #variational inference algorithm for hierarchical Multinomial Dirichlet model
        if(!stan){
          VB=VB_hierMD(summary_data$counts,kappa, tau, nu, Xdim, Ydim, dirichlet_prior, s)
          theta_VB=t(t(VB$nu)/colSums(VB$nu))
        }else{
          VB_stan=VB_stan_hierMD(summary_data$xValues,summary_data$yValues,Xdim, Ydim)
          theta_VB=matrix(colMeans(VB_stan$theta),Xdim,Ydim)
        }
        theta=theta_VB

        # Convert the matrix theta_VB into an array with dimension equal to the numer of parents + 1
        tab=convert2CPT(theta_VB,meaning_xStates,data,node,dag,Xdim,Ydim)
        # Set parameters of tab
        dim(tab)=dim(tab)
        dimnames(tab)[[1]]=levels_old[[node]]
        for(j in 1:(length(dag$nodes[[node]]$parents))){
          dimnames(tab)[[j+1]]=levels_old[[dag$nodes[[node]]$parents[j]]]
        }
        names(dimnames(tab))=c(node,dag$nodes[[node]]$parents)
      }

      # Set class
      class(tab)="table"
      class = ifelse(is(data[, node], "ordered"), "bn.fit.onode", "bn.fit.dnode")

      structure(list(node = node, parents = parents, children = children,
                     prob = tab, joint=theta), class = class)
    }

  }

  #in multi-domain CPT estimation the variable group is added to the dag as parent variable for all the variables
  if(method=="multi-domain"){
    dag_hier=hier.bayes(dag, group)
    dag_orig=dag
    dag=dag_hier
  }else{
    dag_orig=dag
  }

  # Fit the parameters of each node
  fitted = sapply(names(dag_orig$nodes), fit, alpha=alpha, simplify = FALSE)

  # Set class: any additional class of the original bn object is preserved
  orig.class = class(dag)
  class = c(orig.class[orig.class != "bn"], "bn.fit", "bn.fit.dnet")
  # The training node label from Bayesian network classifiers is preserved
  fitted = structure(fitted, class = class, training = dag$learning$args$training)

  return(fitted)
}



