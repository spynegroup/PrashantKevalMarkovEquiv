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
library(gRain)

predict.bn=function(dati_test,fit_bn,class=NULL){
  P=NULL
  Ptrue=NULL
  if(!is.null(class)){
    id_class=which(names(dati_test)==class)
    names_evidence=names(dati_test)[-id_class]
    dati_test_evidence=dati_test[,-id_class]
  }else{
    names_evidence=names(dati_test)
    dati_test_evidence=dati_test
  }
  junction <- compile(as.grain(fit_bn))

  if(!is.null(class)){
  for(Drow in 1:(dim(dati_test)[1])){
    class_true=dati_test[Drow,class]
    set_junction <- setEvidence(junction, nodes = names_evidence, states = sapply(dati_test_evidence[Drow,],as.character))
    Ptemp=querygrain(set_junction, nodes = class)[[class]]
    P=rbind(P,Ptemp)
    Ptrue=c(Ptrue,Ptemp[as.numeric(dati_test[Drow,class])])
  }
  lev=levels(dati_test[,class])
  p_most_prob=apply(P,1,max)
  if(sum(is.nan(p_most_prob))==0){
    most_prob=apply(P,1,which.max)
  }else{
    most_prob=rep(NaN,dim(P)[1])
    most_prob[!is.nan(p_most_prob)]=apply(P[!is.nan(p_most_prob),],1,which.max)
  }
  dati_pred=data.frame(true=dati_test[,class], pred=lev[most_prob],p_true=Ptrue, p_pred=p_most_prob, P ,row.names = NULL)
  dati_pred[,2]=factor(dati_pred[,2],levels=levels(dati_pred[,1]),ordered=is.ordered(dati_pred[,1]))
  acc=mean(dati_pred[,1]==dati_pred[,2])

  names(dati_pred)[5:length(names(dati_pred))]=lev
  cont=table(dati_pred[,1],dati_pred[,2])

  log_loss=sum(log(dati_pred[,3]))


  ROC_vec=NULL
  for(l in lev){
    if(rowSums(cont)[l]>0){
      ROC_vec=c(ROC_vec,wilcox.test(dati_pred[dati_pred$true==l,l],dati_pred[dati_pred$true!=l,l])$statistic/(rowSums(cont)[l]*(sum(cont)-rowSums(cont)[l])))
    }else{
      ROC_vec=c(ROC_vec,0)
    }
  }
  ROC=sum(ROC_vec*rowSums(cont)/sum(cont))

  return(list(dati_pred=dati_pred,acc=acc, cont=cont, log_loss=log_loss,ROC=as.numeric(ROC)))
  }else{
    log_loss=0
    for(Drow in 1:(dim(dati_test)[1])){
      set_junction <- setEvidence(junction, nodes = names_evidence, states = sapply(dati_test_evidence[Drow,],as.character))
      log_loss=log_loss+log(pEvidence(set_junction))
    }
    return(log_loss)
  }
}

predict.RF_all=function(dati_test,RF,class){
  lev=levels(dati_test[,class])
  prob_RF=predict.RF(RF,dati_test,type="prob")
  class_pred_RF=predict.RF(RF,dati_test)
  p_most_prob=apply(prob_RF,1,max)
  prob_RF_extended=cbind(prob_RF,matrix(0,dim(prob_RF)[1],length(lev)-dim(prob_RF)[2]))
  colnames(prob_RF_extended)=unique(c(colnames(prob_RF),lev))
  #prob_RF[1,as.character(dati_test[1,class])]
  Ptrue=NULL
  for(Drow in 1:(dim(dati_test)[1])){
    Ptrue=c(Ptrue,prob_RF_extended[Drow,as.character(dati_test[Drow,class])])
  }

  dati_pred=data.frame(true=dati_test[,class], pred=class_pred_RF,p_true=Ptrue, p_pred=p_most_prob,prob_RF_extended, row.names = NULL)
  dati_pred[,2]=factor(dati_pred[,2],levels=levels(dati_pred[,1]),ordered=is.ordered(dati_pred[,1]))
  acc=mean(dati_pred[,1]==dati_pred[,2])
  cont=table(dati_pred[,1],dati_pred[,2])
  names(dati_pred)[5:length(names(dati_pred))]=colnames(prob_RF_extended)

  log_loss=-sum(log(dati_pred[,3]))


  ROC_vec=NULL
  for(l in lev){
    if(rowSums(cont)[l]>0){
      ROC_vec=c(ROC_vec,wilcox.test(dati_pred[dati_pred$true==l,l],dati_pred[dati_pred$true!=l,l])$statistic/(rowSums(cont)[l]*(sum(cont)-rowSums(cont)[l])))
    }else{
      ROC_vec=c(ROC_vec,0)
    }
  }
  ROC=sum(ROC_vec*rowSums(cont)/sum(cont))


  precision=NULL
  recall=NULL
  for(l in lev){
    precision=c(precision,sum(dati_pred[dati_pred$true==l,"pred"]==l)/sum(dati_pred$pred==l))
    recall=c(recall, sum(dati_pred[dati_pred$true==l,"pred"]==l)/sum(dati_pred$true==l))
  }
  PR=rbind(precision,recall)
  colnames(PR)=lev
  return(list(dati_pred=dati_pred,acc=acc, cont=cont, log_loss=log_loss,ROC=as.numeric(ROC),prec_recall=PR))
}
