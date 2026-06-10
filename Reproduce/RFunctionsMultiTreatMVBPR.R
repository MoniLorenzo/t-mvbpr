#R Script Multi-tratment MultiView Bayesian Profile Regression
#author lorenzo moni


# Definitions of some usefull functions




#1. Compute the optimal posterior prediction for a categorical variables
# from a dosterior sample. Input: matrix n x nMCMC (n: numbr of units,  nMCMC number of mcmc sim), K: number of categories
OptPostCategorical=function(PostDraws, maxCat=NULL){
  if(is.null(maxCat)){maxCat=max(PostDraws)}

  if(maxCat<0 || maxCat< max(PostDraws)){stop("maxCat must be greather than 0 and >= max(PostDraws)")}


  nobs=dim(PostDraws)[1]
#  nMCMC=dim(PostDraws)[2]

 # ferqtable=apply(PostDraws, 1, function(ithunitPosteriordraws){

#    return(which.max(table( factor(ithunitPosteriordraws, levels = 0:maxCat) )/nMCMC ))} )

  Maxs=rep(-Inf, nobs)
  Result=rep(NA, nobs)
  for (k in 0:maxCat) {

    Countk=rowSums(PostDraws==k)
    #Chenge the max if count_i > max_i
    indexmax=Countk>Maxs
    Maxs[indexmax]=Countk[indexmax]
    Result[indexmax]=k

  }

  return(Result)
}




#2. this function build a list with the default hyper-parameters for the priors,
#in the specific format that is usable in the cpp function
genDefaultPrior=function(Nvie,Maxcategories,Ry, singlevalue=NULL, nprogn){
  l=list()
  cat("default hyperparameters for the discrete variables prior: \n")
  cat("1/Rd all entries or all <singlevalue> \n")
  D=length(Maxcategories)
  maxc=max(Maxcategories)
  M=matrix(-10,nrow = maxc,ncol = D)
  for (d in 1:D) {
    if(is.null(singlevalue)){
      M[1:Maxcategories[d], d]=rep(1/Maxcategories[d] ,Maxcategories[d])
    }else{
      M[1:Maxcategories[d], d]=rep(singlevalue ,Maxcategories[d])
    }
  }
  l[["a"]]=M
  cat("default hyperparameters view allocation variables: \n")
  l[["HParViewallocationlog"]]=log(matrix( rep(1/Nvie,Nvie*D), ncol = D, nrow = Nvie))

  cat("default hyper-hyperparameters hyperprior on aalpha_v:  \n")
  #prior alphaDP~Gamma(Shape, scale)
  l[["Hyp2AB"]]=c(2,4) #a=2, b=4

  cat("hyperparameters theta_k:  \n")
  lthata=list()
  lthata[["mean"]]=rep(0,Ry-1)
  lthata[["scale"]]=rep(2.5,Ry-1)
  lthata[["df"]]=rep(7,Ry-1)

  l[["HParthetak"]]=lthata #a=2, b=4



  cat("hyperparameters betas:  \n")
  lbetas=list()
  lbetas[["mean"]]=matrix(0,ncol=Ry-1,nrow =nprogn)
  lbetas[["scale"]]=matrix(2.5,ncol=Ry-1,nrow =nprogn)
  lbetas[["df"]]=matrix(7,ncol=Ry-1,nrow =nprogn)

  l[["HParbetas"]]=lbetas #a=2, b=4


  return(l)

}




Relabel=function(vect, gapscheck=FALSE){
  if(!is.numeric(vect)){stop()}
  if(any(vect < 0 )){stop("Negative values not allowed")}

  minvalue=min(vect)
  maxvalue=max(vect)

  #check for gaps
  gaps=setdiff(seq(min(vect), max(vect)), sort(unique(vect)))
  if(gapscheck){
   if(length(gaps) > 0){cat("vector has gaps: ",gaps,  " \n");stop()}
  }
  if(minvalue>0){
    cat("relabeled")
    return(vect-minvalue)
 #   return(vect)
  }else{
    return(vect)
  }
}



#3. MAIN FUNCTION: envelop the cpp MVBPRmultitreat and compute all the structure needed to call it
MultiTreatmentMVBPR=function(Dataset,DatasetPostPred=NULL,
                             nameTreat,
                             nameY,
                             namesPred=NULL,
                             namesProg=NULL){

  #if namesPred is NULL, the we consider all the cols of Dataset with "d" in the name
  if(is.null(namesPred)){
    namesPred=colnames(Dataset)[grep("^d", colnames(Dataset))]
  }
  #if namesProg is NULL, the we consider all the cols of Dataset with "w" in the name
  if(is.null(namesProg)){
    namesProg=colnames(Dataset)[grep("w", colnames(Dataset))]
  }

  #
  #count number of treatments, if only one treatment simple MultiViewBPR
  if(is.null(nameTreat)){stop("nameTreat not given")}

  TreatInceces=Dataset[,nameTreat]
  nTreat=length(unique(TreatInceces))

  #construct list containing the data
  ListData=list()

  ListData[["fitdata"]]=list()

  #DATA TO FIT THE MODEL
  #this save the cpp based index of the treatment, only used for debug and to match the
  #native treatment index on Dataset to the cpp treatment tracking
  cppit=0

  
  #relabeling if needed
  Datasetrel=Dataset[ ,namesPred, drop=FALSE]#apply(Dataset[ ,namesPred, drop=FALSE],2, Relabel)
  YRel=Dataset[ ,nameY]#Relabel(Dataset[ ,nameY])
  
  for (t in sort(unique(TreatInceces))) {
    listtreat_t=list()

    
    listtreat_t[["Y"]]=YRel[TreatInceces==t]

    #apply relabeling
    listtreat_t[["XDiscPred"]]=Datasetrel[TreatInceces==t, , drop=FALSE]#apply(Dataset[TreatInceces==t ,namesPred, drop=FALSE],2, Relabel)

    #figurative matrix: not used but needed [nx0]  since continuos predictive variables not timplemented in cpp
    listtreat_t[["XContPred"]]=Dataset[TreatInceces==t ,0,   drop=FALSE]

    listtreat_t[["WProgn"]]=Dataset[TreatInceces==t ,namesProg,  drop=FALSE]

    listtreat_t[["Info"]]=paste0("Treatment: ", t, " - cpp index of the treatment:", cppit)

    ListData[["fitdata"]][[t]]=listtreat_t
    cppit=cppit+1
    rm(listtreat_t)
  }

  #OUT OF SAMPLE DATA TO COMPUTE THE POSTERIOR PREDICTIVE
  ListData[["postpredictivedata"]]=list()

  if(is.null(DatasetPostPred)){
    ListData[["postpredictivedata"]][["XDiscPred"]]=Dataset[  ,namesPred, drop=FALSE]#apply(Dataset[  ,namesPred, drop=FALSE],2, Relabel)

    #figurative matrix: not used but needed [nx0]  since continuos predictive variables not timplemented in cpp
    ListData[["postpredictivedata"]][["XContPred"]]=Dataset[  ,0,   drop=FALSE]

    ListData[["postpredictivedata"]][["WProgn"]]=Dataset[  ,namesProg,  drop=FALSE]


  }else{



    #listtreat_pred[["Y"]]=Relabel(Dataset[TreatInceces==t ,nameY])

    #apply relabeling
    ListData[["postpredictivedata"]][["XDiscPred"]]=DatasetPostPred[,namesPred, drop=FALSE]#apply(DatasetPostPred[,namesPred, drop=FALSE],2, Relabel)

    #figurative matrix: not used but needed [nx0]  since continuos predictive variables not timplemented in cpp
    ListData[["postpredictivedata"]][["XContPred"]]=DatasetPostPred[,0,   drop=FALSE]

    ListData[["postpredictivedata"]][["WProgn"]]=DatasetPostPred[ ,namesProg,  drop=FALSE]


  }

  return(ListData)
}





#function to initialize the gammas

InitializerGammas=function(Listdata, p=0.05){
  ntreat=length(Listdata$fitdata)
  
  Pvalues=matrix(NA,ncol = ntreat, nrow = ncol(Listdata$fitdata[[1]]$XDiscPred))
  
  for (t in 1:ntreat) {
    Xt=Listdata$fitdata[[t]]$XDiscPred
    Wt=Listdata$fitdata[[t]]$WProgn
    Yt=Listdata$fitdata[[t]]$Y
    
    
    
    mod0=lm(Yt ~ Wt)
    
    
    for(d in 1:ncol(Xt)){
      
      Xtfact=try(as.factor(Xt[,d]))
      if(length(levels(Xtfact))!=1){
        mod_d=lm(Yt ~ Wt+Xtfact)
        Pvalues[d,t]=anova(mod0,mod_d)$`Pr(>F)`[2]
        
        
      }else{
        Pvalues[d,t]=1
      }
      
      
    }
    
  }
  
  return(Pvalues<p)
  
  
}










#4. returns the MAP for the posterior partition in each view
OptPostPart=function(ZarrPost){
  res=apply(ZarrPost,2,function(z_v){
    postsim=mcclust::comp.psm(t(z_v))
    mcclust::maxpear(postsim)$cl
  })

  colnames(res)=paste0("Zopt",1:ncol(res))
  res
}

#5. plot marginal posterior probs view allocation variables
PlotMPPViewAlloc=function(PostDraws, TypeofPlot=1,maxCat=NULL){#TypeofPlot=1 on plot # TypeofPlot=2 plot with facet_wrap
  require(ggplot2)
  if(is.null(maxCat)){maxCat=max(PostDraws)}

  if(maxCat<0 || maxCat< max(PostDraws)){stop("maxCat must be greather than 0 and >= max(PostDraws)")}
  nobs=dim(PostDraws)[1]
  nmcmc=dim(PostDraws)[2]

  Result=matrix(NA,nrow = nobs,ncol =maxCat+1 )
  colnames(Result)=paste0("View",0:maxCat)
  for (k in 0:maxCat) {

    Result[,k+1]=Countk=rowSums(PostDraws==k)
  }
  Result=Result/nmcmc
  Result=reshape2::melt(Result)

  if(TypeofPlot==1){
    ggplot(Result, aes(x = factor(Var1), y = value, fill = Var2)) +
      geom_bar(stat = "identity", width = 0.5) +  # narrower bars
      xlab("gamma_d") +
      ylab("Probability") +  scale_y_continuous(expand = c(0, 0)) +  theme_minimal()
  }


  if(TypeofPlot==1){

    ggplot(Result, aes(x = factor(Var1), y = value)) +
      geom_bar(stat = "identity", width = 0.2) +
      xlab("gamma_d") +
      ylab("Post marginal prob") +
      scale_y_continuous(expand = c(0, 0)) +
      facet_grid(Var2~.,) +
      #theme_minimal() +
      theme(
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank()
      )

  }

}


#6. plot betas and prior
PlotBeta=function(BetasArr, Listhyperparam=list(mu = 0, sigma = 2.5, nu = 7)){
  require(tidyverse)
  require(ggplot2)
  #names Beta in futuro mettere un altra funzione che da i corretti nom
  dimnames(BetasArr)= list(paste0("p",1:dim(BetasArr)[1]),
                           paste0("r",1:dim(BetasArr)[2]))


  if(is.null(dimnames(BetasArr))){stop("no names for betas")}

  #plot
  BetasArr %>% reshape2::melt() %>% rename(Progn=Var1,Category=Var2,mcmcindex=Var3) %>%
    ggplot(.,aes(x=value))+
    geom_density()+
    facet_wrap(Progn~Category)+
    stat_function(fun=LaplacesDemon::dst, args = Listhyperparam,
                  color="black",linetype="dotted")



}


#7. Estimated utility
PosteriorUtility=function(PostPredictivetreatment_t, weights=c(0,40,100)){
  if(max(PostPredictivetreatment_t)+1!=length(weights)){stop("n category respone != than weights length")}
  Utility_hat=0
  for( k in 1:length(weights)){
    Utility_hat=Utility_hat+rowSums(PostPredictivetreatment_t==(k-1))*weights[k]
  }
return(Utility_hat/ncol(PostPredictivetreatment_t))
}


#8. EstimatedOptimalTreatment
EstimOptTreatment=function(PostPredictiveAlltrt, weights=c(0,40,100)){
  MAPytrts=NULL
  Utilitys=NULL
  for( t in names(PostPredictiveAlltrt)[grep("Treatment",names(PostPredictiveAlltrt))]){
    MAPytrts=cbind(MAPytrts, OptPostCategorical(PostPredictiveAlltrt[[t]]$PostPredictiveY))
    Utilitys=cbind(Utilitys, PosteriorUtility(PostPredictiveAlltrt[[t]]$PostPredictiveY,weights ))
  }

  Estimatedopttrt=apply(Utilitys,1,which.max)

  #Outcome: y1 if utility 1> unility 2 ,ow y2
  Outpred=MAPytrts[cbind(1:nrow(MAPytrts), Estimatedopttrt)]
  return(list(EstimetedAssTrt=Estimatedopttrt, PredOutcome=Outpred, Utility=Utilitys))
}

#9. Compute Empirical Summary Measure EMS  
ComuteESM=function(PredOptTrt, PredResponse, ActualAssTrt, TrueResp, responders=c("2")){
  
  yresp=(TrueResp%in%responders)*1
  
  #weights 
  pact_asstrt1=mean(ActualAssTrt==1)
  pact_asstrt2=mean(ActualAssTrt==2)

  #prop respondents under rnd all
  P_resp_rnd=mean(yresp) 

  #prop patient ass trt1
  precom_asstrt1=mean(PredOptTrt==1)
  precom_asstrt2=mean(PredOptTrt==2)
  
  #prop ass 1 recomanded 1
  p_rec_ass1=mean(yresp[(ActualAssTrt==1 & PredOptTrt==1) ])
  p_rec_ass2=mean(yresp[(ActualAssTrt==2 & PredOptTrt==2) ])
  
  ESM=(p_rec_ass1*pact_asstrt1+p_rec_ass2*pact_asstrt2)-P_resp_rnd
  
  
  }


##esm di pedone-----
#ESM
# non ho definito come respondent anche i partial responent
#myoutot <- as.numeric(matchRTComp[,9])#simdata$yord[[k]][131:158,]
#mytab <- cbind(myass = predAPT_all[,3], rndass = trtsgn, resp = as.numeric(myoutot>2))
#pred1 <- subset(mytab, mytab[,1]==1)
#table1 <- table(pred1[,3],pred1[,2])
#pred2 <- subset(mytab, mytab[,1]==2)
#table2 <- table(pred2[,3], pred2[,2])
#p1 <- sum(table1)/(sum(table1)+sum(table2))
#p2 <- sum(table2)/(sum(table1)+sum(table2))

#if(length(table1) == 4){
#  crt1 <- table1[2,1]/sum(table1[,1])
#}
#if(length(table1) < 4){
#  crt1 <- as.numeric(row.names(table1))
#}

#if(length(table2) == 4){
#  crt2 <- table2[2,2]/sum(table2[,2])
#}
#if(length(table2) < 4){
#  crt2 <- as.numeric(row.names(table2))
#}
#ESM <- c(crt1*p1 + crt2*p2 - sum(as.numeric(myoutot>2))/npat)
### summary meaures

#######-----





#10. Plot estimated utility
PlotTrueveEstdiffUtility=function(Dataset, namesutility=c("utility1","utility2"), PredUtility){
  require(tidyverse)
  diffU1U2true=Dataset[,namesutility[1]]-Dataset[,namesutility[2]]

  #
  order=order(diffU1U2true)

  diffU1U2pred=PredUtility[,1]-PredUtility[,2]


  tibble(DiffUtility=c(diffU1U2true[order],diffU1U2pred[order]), index=c(1:length(diffU1U2true),1:length(diffU1U2true)),
         Legend=c(rep("True", length(diffU1U2true)), rep("Predicted", length(diffU1U2true)))) %>%
    ggplot(., aes(x=index, y=DiffUtility,color=Legend))+
    geom_line()
}



#10. Plot estimated utility
PlotTrueveEstdiffUtilityMethods=function(Dataset, namesutility=c("utility1","utility2"), PredUtilitylist){
  require(tidyverse)
  diffU1U2true=Dataset[,namesutility[1]]-Dataset[,namesutility[2]]

  #
  order=order(diffU1U2true)


  diffU1U2pred=lapply(PredUtilitylist, function(l){l[,1]-l[,2]})

  #reordering and melt
  diffU1U2pred= reshape2::melt(lapply(diffU1U2pred, function(l){l[order]}))


  tibble(DiffUtility=c(diffU1U2true[order],diffU1U2pred[,1]), reorderedindex=c(1:length(diffU1U2true),
                                                                               rep(1:length(diffU1U2true),length(unique(diffU1U2pred[,2])))),
         Legend=c(rep("True", length(diffU1U2true)), diffU1U2pred[,2])) %>%
    ggplot(., aes(x=reorderedindex, y=DiffUtility,color=Legend))+
    geom_line()
}



##Compute posterior Y for each treatment STAN
postTrtinfStan=function(PostY=A$ypred, weights=c(0,40,100)){
  n=nrow(PostY)
  #map yt1 and yt2
  Yt1=OptPostCategorical(PostY[1:(n/2),])
  Yt2=OptPostCategorical(PostY[(n/2+1):n,])

  #estimated utility
  u1=PosteriorUtility(PostY[1:(n/2),])
  u2=PosteriorUtility(PostY[(n/2+1):n,])

  #opt trearment
  Otrt=apply(cbind(u1,u2), 1, which.max)


  #Outcome: y1 if utility 1> unility 2 ,ow y2
  predoutc=Yt1*(Otrt==1)+Yt2*(Otrt==2)

  return(list(y1pred=Yt1, y2pred=Yt2, PredOptTrt=Otrt, PredOutcome=predoutc, Utility=cbind(u1,u2)))
}





###




































































