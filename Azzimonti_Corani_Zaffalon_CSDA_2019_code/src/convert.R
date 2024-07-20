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

convert2array=function(theta,meaning_zStates,dati,node,dag,I,J){
  #theta=matrix(colMeans(sample_theta),I,J)

  dimvec=length(levels(as.factor(dati[,node])))
  dimnames = list(levels(as.factor(dati[,node])))

  for(parent in dag$nodes[[node]]$parents){
    dimvec=c(dimvec,length(levels(as.factor(dati[,parent]))))
    dimnames[[length(dimnames)+1]]=levels(as.factor(dati[,parent]))
  }
  MAT=array(0,dim=dimvec,dimnames=dimnames)


  meaning_mat=t(data.frame(sapply(meaning_zStates,strsplit,split="/")))
  colnames(meaning_mat)=rev(dag$nodes[[node]]$parents)
  rownames(meaning_mat)=NULL


  for(id_row in 1:dim(meaning_mat)[1]){
    meaning_mat[id_row,]
    stringa="MAT["
    for(d in 2:length(dim(MAT))){
      id=which(dimnames(MAT)[[d]]==meaning_mat[id_row,length(dim(MAT))-d+1])
      stringa=paste(stringa, id, sep=",")
    }
    stringa=paste(stringa,"] = theta[,id_row]",sep="")
    eval(parse(text=stringa))
  }

  return(MAT)
}




convert2CPT=function(theta,meaning_xStates,dati,node,dag,I,J){
  #theta=matrix(colMeans(sample_theta),I,J)

  dimvec=length(levels(as.factor(dati[,node])))
  dimnames = list(levels(as.factor(dati[,node])))

  for(parent in dag$nodes[[node]]$parents){
    dimvec=c(dimvec,length(levels(as.factor(dati[,parent]))))
    dimnames[[length(dimnames)+1]]=levels(as.factor(dati[,parent]))
  }
  MAT=array(0,dim=dimvec,dimnames=dimnames)



  meaning_mat=t(data.frame(sapply(meaning_xStates,strsplit,split="/")))
  meaning_mat=cbind(meaning_mat[,-1],meaning_mat[,1])
  colnames(meaning_mat)=c(rev(dag$nodes[[node]]$parents)[-1],node)
  rownames(meaning_mat)=NULL

  for(id_row in 1:dim(meaning_mat)[1]){
    meaning_mat[id_row,]
    stringa="MAT["
    for(d in 1:(length(dim(MAT))-1)){
      id=which(dimnames(MAT)[[d]]==meaning_mat[id_row,length(dim(MAT))-d])
      stringa=paste(stringa, id, sep="")
      stringa=paste(stringa,",", sep="")
    }
    stringa=paste(stringa,"] = theta[id_row,]",sep="")
    eval(parse(text=stringa))
  }


  MAT_cond=colSums(MAT)

  MAT_temp=MAT
  for(id_row in 1:dim(MAT)[1]){
    stringa1="MAT_temp[id_row"
    stringa2="MAT[id_row"
    for(d in 1:(length(dim(MAT))-1)){
      stringa1=paste(stringa1,",", sep="")
      stringa2=paste(stringa2,",", sep="")
    }
    stringa1=paste(stringa1,"]", sep="")
    stringa2=paste(stringa2,"]", sep="")
    stringa=paste(stringa1,"=",stringa2,"/MAT_cond",sep="")
    eval(parse(text=stringa))
    #MAT_temp[id_row,,]=MAT[id_row,,]/MAT_cond
  }


  return(MAT_temp)
}


