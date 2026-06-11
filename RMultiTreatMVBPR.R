#R Script Multi-tratment MultiView Bayesian Profile Regression
#@author: Lorenzo Moni, Silvia Liverani, Alberto Cassese, Francesco Claudio Stingo  


require(Rcpp)
require(viridis)
require(mcclust)
require(mltools)

#install.packages("BiocManager"); BiocManager::install("ComplexHeatmap")
require(ComplexHeatmap)


#sourceCpp("cpprev2DEBUG.cpp")
sourceCpp("~/Downloads/drive-download-20260610T160654Z-3-001/cppCodetMVBPR.cpp")

# Definitions of some useful functions




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
  names(Result)=rownames(PostDraws)
  
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
genDefaultPrior=function(Nvie,Maxcategories,Ry, singlevalue=NULL, HyperparDP=NULL,
                         nprogn, verbose=FALSE){
  l=list()
  if(verbose){
    cat("default hyperparameters for the discrete variables prior: \n")
    cat("1/Rd all entries or all <singlevalue> \n")
    
  }
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
  if(verbose){
    cat("default hyperparameters view allocation variables: \n")}
  l[["HParViewallocationlog"]]=log(matrix( rep(1/Nvie,Nvie*D), ncol = D, nrow = Nvie))

  if(verbose){
    cat("default hyper-hyperparameters hyperprior on aalpha_v:  \n")}
  #prior alphaDP~Gamma(Shape, scale)
  
  if(is.null(HyperparDP)){
    l[["Hyp2AB"]]=c(2,1) #a=2, b=4
  }else{
    if(all(HyperparDP>0) & length(HyperparDP)==2){
      l[["Hyp2AB"]]=HyperparDP
    }else{
      stop("hyper-hyperparameters DPs must be 2 non-negative: c(Shape, scale) ")
    }
  }

  if(verbose){
    cat("hyperparameters theta_k:  \n")}
  lthata=list()
  lthata[["mean"]]=rep(0,Ry-1)
  lthata[["scale"]]=rep(2.5,Ry-1)
  lthata[["df"]]=rep(7,Ry-1)


  l[["HParthetak"]]=lthata  



  if(verbose){
    cat("hyperparameters betas:  \n")}
  lbetas=list()
  lbetas[["mean"]]=matrix(0,ncol=Ry-1,nrow =nprogn)
  lbetas[["scale"]]=matrix(2.5,ncol=Ry-1,nrow =nprogn)
  lbetas[["df"]]=matrix(7,ncol=Ry-1,nrow =nprogn)

  l[["HParbetas"]]=lbetas #a=2, b=4


  return(l)

}


#Convert factor in numeric 0:nlevel(x)-1
FcttoNum=function(Input, Mapping=NULL ){
  if(!is.factor(Input)){stop("Outcome must be factor")}
  
  levs=levels(Input)
  nlev=nlevels(Input)
  
  #default relabeling -- consider the order of levels(dataset[[nameOutcome]])
  if(is.null(Mapping)){
    Mapping = setNames(0:(nlev - 1), levs )
    
  }else{# for personalized mapping only check validity of mapping
    if(!setequal(names(Mapping), levs) || !all(sort(Mapping) == 0:(nlev - 1))){
      stop("Invalid mapping")
    }
  }
  #relabel
  Output= unname(Mapping[as.character(Input)])
  
  #check one-to-one mapping
  tab = table(Input, Output)
  #print(tab)
  
  if (!(all(apply(tab, 1, function(x) sum(x != 0) <= 1)) &&
        all(apply(tab, 2, function(x) sum(x != 0) == 1)))) {
    stop("Relabeling is not one-to-one.")
  }
  return( list(Vect=Output, Mappingused=Mapping))
}
#relabel outcome (convert fct in numeric from 0:nlev(outcome)-1)



 
 


#3. MAIN FUNCTION: envelop the cpp MVBPRmultitreat and compute all the structure needed to call it
##########à
tMVBPR=function(dataset, predictiveDataset=NULL, 
                Outcomemapping=NULL, 
                NumberofViews=3,
                ListPrior=NULL,
                ListInitialValues=list(Numbinitialclusters=10),
                nameTreatment="Treatment", 
                nameOutcome="Outcome",
                namesPredictive=NULL,
                namesPrognostic=NULL,
                MCMCfinalSampleSize=20,
                MCMCburnin=10,
                MCMCThinning=2,
                ListMCMCparam=list(PropVarTheta=0.1, InitVarBeta=2, RelViewNeal8Mparam=10),
                HowgammasInit=-2
){
  #####

  listRelabeling=list()
  if(class(dataset)!="data.frame"){stop("dataset must be a data.frame")}
  
  if(!all(c(nameOutcome,namesPrognostic, namesPredictive) %in% colnames(dataset))){
    stop("dataset lacks one or more required columns")}
  
  
  
  #relabel outcome (convert fct in numeric from 0:nlev(outcome)-1)
  Relabtemp=FcttoNum(dataset[[nameOutcome]],Mapping = Outcomemapping)
  dataset[[nameOutcome]]=Relabtemp$Vect
  listRelabeling[[nameOutcome]]=Relabtemp$Mappingused
  
  
  #set n of categories outcome
  ncategoriesOutcome=length(Relabtemp$Mappingused)
  
  #relabel predictive variables
  ncategoriesPred=setNames(rep(NA,length(namesPredictive)), namesPredictive)
  for (predvarnam in namesPredictive) {
    
    Relabtemp=FcttoNum(dataset[[predvarnam]])
    dataset[[predvarnam]]=Relabtemp$Vect
    listRelabeling[[predvarnam]]=Relabtemp$Mappingused
    
    ncategoriesPred[predvarnam]= length(Relabtemp$Mappingused)
    
    
  }
  
  
  #Build dummies for factor prognostivariables
  oldnamesPrognostic=namesPrognostic
  #prognostic variables
  contrastsarg = lapply(
    dataset[oldnamesPrognostic],
    function(x) if (is.factor(x)) contrasts(x)
  )
  contrastsarg=contrastsarg[ !unlist(lapply(contrastsarg, is.null)) ]
  
  
  ModmatProgn=model.matrix(~.,dataset[oldnamesPrognostic], 
                           contrasts.arg =contrastsarg ) 
  Prognvarrecoded = as.data.frame( ModmatProgn[, setdiff(colnames(ModmatProgn),
                                                         "(Intercept)"), 
                                               drop = FALSE]  )
  #reference levels
  listRelabeling[["RefPrognosticfct"]]=lapply(dataset[oldnamesPrognostic], function(x) if (is.factor(x)) levels(x)[1] else NULL)
  
  dataset[namesPrognostic]=NULL
  dataset[colnames(Prognvarrecoded)]=Prognvarrecoded
  
  #change the names of prognostic variables with those recoded
  namesPrognostic=colnames(Prognvarrecoded)
  
  
  
  #Predictive dataset
  if (!is.null(predictiveDataset)) {
    if(!all(c(oldnamesPrognostic, namesPredictive) %in% colnames(predictiveDataset))) {
      stop("predictiveDataset lacks required columns: namesPrognostic, namesPredictive must be the same in predictiveDataset and dataset")
    }
    
    #relabel predictive variables 
    for (predvarnam in namesPredictive) {
      mapping = listRelabeling[[predvarnam]]
      
     predictiveDataset[[predvarnam]] =
        as.integer(mapping[as.character(predictiveDataset[[predvarnam]]) ])
     
     #check if there are more levels in predictive dataset 
     if(any(is.na(predictiveDataset[[predvarnam]]))){
       stop(paste0("In predictiveDataset the variable ", predvarnam, "has more levels than in dataset" ))
       
     }
      
    }
    
    #prognostic variables
    ModmatPredProgn=model.matrix(~., predictiveDataset[oldnamesPrognostic],
                                  contrasts.arg = contrastsarg)
    
    PrognvarrecodedPred = as.data.frame( ModmatPredProgn[, setdiff(colnames(ModmatPredProgn),
                                                           "(Intercept)"), 
                                                 drop = FALSE]  )
    predictiveDataset[oldnamesPrognostic]=NULL
    predictiveDataset[colnames(PrognvarrecodedPred)]=PrognvarrecodedPred

    }
  
  
  
  
  ##
  rm(ModmatProgn,Prognvarrecoded, Relabtemp, oldnamesPrognostic, ModmatPredProgn,PrognvarrecodedPred )
  
  #create suitable datastructure for the cpp function
  Datalist=CreateDataStructureforCpp(Dataset = dataset, DatasetPostPred=predictiveDataset,
                                     nameTreat=nameTreatment,
                                     nameY=nameOutcome,
                                     namesPred =namesPredictive,
                                     namesProg =namesPrognostic )
  #extract n_treatment and number of prognostic variables from the Datalist
  Q=dim(Datalist$fitdata[[1]]$WProgn)[2]
  nTreatments=length(Datalist$fitdata)
  
  D=dim(Datalist$fitdata[[1]]$XDiscPred)[2]
  
  #create suitable datastructure for the prior hyperparameters
  if(is.null(ListPrior)){
    PriorList=genDefaultPrior(Nvie =NumberofViews,HyperparDP = c(5,5),
                              Maxcategories =ncategoriesPred ,singlevalue = 1,
                              Ry =ncategoriesOutcome,
                              nprogn = Q)
  }
  ##
  #PriorList$HParViewallocationlog=(PriorList$HParViewallocationlog*0)+log(c(0.90/3,0.1,0.90/3,0.90/3))
  
  
  #initialization
  ainit=rep(6,NumberofViews-1) #alphas DPs initialization
  Betainit=matrix((0), nrow =Q, ncol = ncategoriesOutcome-1)
  
  ##
  ginit=InitializerGammas(Datalist)*1#matrix(c(rep(0,15),rep(1,5)), ncol=2,nrow = 20)
  ginit[is.na(ginit)]=0
  
  tMVBPRcpp=MVBPRmultitreat(DataList =Datalist ,
                            ListPrior = PriorList,
                            numbtreatments =nTreatments,
                            numbprognosticvar =Q,
                            maximumcatDiscrete = ncategoriesPred,
                            NumberofView = NumberofViews,
                            TypeofXmodel =1, #DO NOT CHANGE 
                            Maxcatresponse =ncategoriesOutcome ,
                            M = ListMCMCparam$RelViewNeal8Mparam,
                            lambdtheta = ListMCMCparam$PropVarTheta,
                            lambdbeta = ListMCMCparam$InitVarBeta,
                            Ninitclusters = ListInitialValues$Numbinitialclusters,
                            InitGammas =HowgammasInit, #-2 
                            InitGammasmat = ginit,##
                            AlphaDPsinit =ainit,
                            BetaMatinit = Betainit,
                            MCMCFinalSamplesize =MCMCfinalSampleSize,
                            Burnin = MCMCburnin,
                            Thinning = MCMCThinning )
  
  
  
  
  tMVBPRcpp=namingresults(tMVBPRcpp)
  
  tMVBPRcpp[["OriginalVariablesCoding"]]=listRelabeling
  return(tMVBPRcpp)
}



#4.
#NEW name: CreateDataStructureforCpp OLDNAME MultiTreatmentMVBPR 
CreateDataStructureforCpp=function(Dataset,DatasetPostPred=NULL,
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
    listtreat_t[["XDiscPred"]]=as.matrix(Datasetrel[TreatInceces==t, , drop=FALSE])  #apply(Dataset[TreatInceces==t ,namesPred, drop=FALSE],2, Relabel)

    #figurative matrix: not used but needed [nx0]  since continuos predictive variables not timplemented in cpp
    listtreat_t[["XContPred"]]=as.matrix(Dataset[TreatInceces==t ,0,   drop=FALSE])

    listtreat_t[["WProgn"]]=as.matrix(Dataset[TreatInceces==t ,namesProg,  drop=FALSE])

    listtreat_t[["Info"]]=paste0("Treatment: ", t, " - cpp index of the treatment:", cppit)

    ListData[["fitdata"]][[t]]=listtreat_t
    cppit=cppit+1
    rm(listtreat_t)
  }

  #OUT OF SAMPLE DATA TO COMPUTE THE POSTERIOR PREDICTIVE
  ListData[["postpredictivedata"]]=list()

  if(is.null(DatasetPostPred)){
    ListData[["postpredictivedata"]][["XDiscPred"]]=as.matrix(Dataset[  ,namesPred, drop=FALSE])

    #figurative matrix: not used but needed [nx0]  since continuos predictive variables not timplemented in cpp
    ListData[["postpredictivedata"]][["XContPred"]]=as.matrix(Dataset[  ,0,   drop=FALSE])

    ListData[["postpredictivedata"]][["WProgn"]]=as.matrix(Dataset[  ,namesProg,  drop=FALSE])


  }else{



    #listtreat_pred[["Y"]]=Relabel(Dataset[TreatInceces==t ,nameY])

    #apply relabeling
    ListData[["postpredictivedata"]][["XDiscPred"]]=as.matrix(DatasetPostPred[,namesPred, drop=FALSE])#apply(DatasetPostPred[,namesPred, drop=FALSE],2, Relabel)

    #figurative matrix: not used but needed [nx0]  since continuos predictive variables not timplemented in cpp
    ListData[["postpredictivedata"]][["XContPred"]]=as.matrix(DatasetPostPred[,0,   drop=FALSE])

    ListData[["postpredictivedata"]][["WProgn"]]=as.matrix(DatasetPostPred[ ,namesProg,  drop=FALSE])


  }

  return(ListData)
}






#frequentest init view allocations 
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


#dataset=readRDS("/home/lorenzo/Documents/UNIFI/PhD/Peoject/tcga/MYData1") 


#Post=readRDS("~/Downloads/ResultLGG70real701")
#Utility function to name the result objects, ie object in returned list
namingresults=function(Post){
  S=Post$Info$MCMCinfo["MCMCFinalSamplesize"]

  namesProg=Post$Info$VariablesNames$XprogNames
  namesPred=Post$Info$VariablesNames$XdiscNames
  
  
  for (treatmentname in names(Post)[grepl( "Treatment",names(Post))]) {
    
    #names cluster allocations
    dimnames(Post[[treatmentname]]$ZZallviews)=list(paste0("i=",1:dim(Post[[treatmentname]]$ZZallviews)[1]),
                                                  paste0("view",1:dim(Post[[treatmentname]]$ZZallviews)[2]),
                                                  NULL)
    
    #names gammas
    dimnames(Post[[treatmentname]]$Gammas)=list(paste0("gamma", namesPred),
                                          NULL)
    
    #names precision parameter DPs
    dimnames(Post[[treatmentname]]$alpha)=list(paste0("alpha",1:dim(Post[[treatmentname]]$alpha)[1]),
                                                    NULL)
    
    #names Unit-specific DP parameters 
    dimnames(Post[[treatmentname]]$Thetaus)=list(paste0("r",1:dim(Post[[treatmentname]]$Thetaus)[1]),
                                           paste0("thetaus_i",1:dim(Post[[treatmentname]]$Thetaus)[2]),
                                           NULL)
    
    #names predictive units: posterior predictive outcome
    dimnames(Post[[treatmentname]]$PostPredictiveY)=list(paste0("itilde=",1:dim(Post[[treatmentname]]$PostPredictiveY)[1]),
                                                   NULL)
  }
  #names betas, ieprognostic variables effects
  dimnames(Post$Beta)=list(paste0("beta_",namesProg),
                           paste0("r",1:dim(Post$Beta)[2]),
                           NULL)
  
  return(Post)
}



#helper function to keep only specific indices of iterations
KeepIters=function(Post,iter){
  S=Post$Info$MCMCinfo["MCMCFinalSamplesize"]
  
  namesProg=Post$Info$VariablesNames$XprogNames
  namesPred=Post$Info$VariablesNames$XdiscNames
  
  
  for (treatmentname in names(Post)[grepl( "Treatment",names(Post))]) {
    
    #cluster allocations
    Post[[treatmentname]]$ZZallviews=Post[[treatmentname]]$ZZallviews[,,iter]
    
    #gammas
    Post[[treatmentname]]$Gammas=Post[[treatmentname]]$Gammas[,iter]
    
    #precision parameter DPs
    Post[[treatmentname]]$alpha=Post[[treatmentname]]$alpha[,iter]
    
    #Unit-specific DP parameters 
    Post[[treatmentname]]$Thetaus=Post[[treatmentname]]$Thetaus[,,iter]
    
    # predictive units: posterior predictive outcome
    Post[[treatmentname]]$PostPredictiveY=Post[[treatmentname]]$PostPredictiveY[,iter]
  }
  # betas, ieprognostic variables effects
  Post$Beta=Post$Beta[,,iter] 
  

  
  return(Post)
}



#Combine multiple MCMC chains into a single list
combineChains=function(PostMultipleChains, iter=NULL){
  nchains=length(PostMultipleChains)
  nameschains=names(PostMultipleChains)
  
  #check if n chains is >1
  if(length(PostMultipleChains)<1){
    stop("More that two chains are required")
    
  }else{
    if(any(names(PostMultipleChains)!=paste0("Chain", 1:length(PostMultipleChains)))){
      stop("Each chain must be named Chain1,Chain2,etc.")
      
    }
      
    
  }
  
  Output=list()
  
  #keep only iter draws
  if(!is.null(iter)){
    PostMultipleChains=lapply(PostMultipleChains, function(chain) KeepIters(chain,iter ))  
  }
  
  #check if OriginalVariablesCoding are coherent across chains
  if(!all(sapply(2:length(PostMultipleChains), function(i) 
    identical(PostMultipleChains[[i]]$OriginalVariablesCoding, 
              PostMultipleChains[[1]]$OriginalVariablesCoding)))){
    stop("OriginalVariablesCoding differ for chains: may indicate 2 different models")
  }
  
  
  
  #check if treatments are coherent across chains
  if(!all(sapply(2:length(PostMultipleChains), function(i) 
    identical(names(PostMultipleChains[[i]])[grep("Treatment",names(PostMultipleChains[[i]]))], 
              names(PostMultipleChains[[1]])[grep("Treatment",names(PostMultipleChains[[1]]))])))){
    stop("Treatments differ for chains: may indicate 2 different models")
  }else{
    Trtnames=names(PostMultipleChains[[1]])[grep("Treatment",names(PostMultipleChains[[1]]))]
  }

  
  
  
  ###combine betas and checks on the structure 
  MCMCiterbetas=lapply(PostMultipleChains,function(x){
    list(dimname=dimnames(x$Beta),
    dims=dim(x$Beta))})
  
  if(!all(sapply(2:length(MCMCiterbetas), function(i) 
    identical(MCMCiterbetas[[i]]$dimname[1:2], 
              MCMCiterbetas[[1]]$dimname[1:2])))){
    stop("Dimesions and/or dimnames of of beta structures  differ for chains: may indicate 2 different models")
  }
  
  #total iterations: save for susequent checks
  totMCMC=sum(sapply(MCMCiterbetas, function(x) x$dims[3]))
  
  #combine the chains: Beta
  
  Output[["Beta"]]=array(NA, dim = c(MCMCiterbetas[[1]]$dims[1:2],totMCMC), 
                         dimnames =MCMCiterbetas[[1]]$dimname[1:2] )
  
  istart=0
  for(c in nameschains){
    iend=istart+MCMCiterbetas[[c]]$dims[3]
    Output[["Beta"]][,, (istart+1):iend]=PostMultipleChains[[c]]$Beta
    istart=iend
  }
  rm(MCMCiterbetas, istart, iend,c)

  
  for (treatmentname in Trtnames) {
    
     ###combine cluster allocations and checks on the structure 
    MCMCiterZZ=lapply(PostMultipleChains,function(x){
      list(dimname=dimnames(x[[treatmentname]]$ZZallviews),
           dims=dim(x[[treatmentname]]$ZZallviews))})
    
    if(!all(sapply(2:length(MCMCiterZZ), function(i) 
      identical(MCMCiterZZ[[i]]$dimname[1:2], 
                MCMCiterZZ[[1]]$dimname[1:2])))){
      stop("Dimesions and/or dimnames of cluster allocations  differ for chains: may indicate 2 different models")
    }else{
      if(  totMCMC!=sum(sapply(MCMCiterZZ, function(x) x$dims[3]))){
        stop("totMCMC iterations differ for beta and ZZ: may indicate 2 different models")
        
      }
      
    }
    #combine the chains: ZZ
    
    Output[[treatmentname]][["ZZallviews"]]=array(NA, dim = c(MCMCiterZZ[[1]]$dims[1:2],totMCMC), 
                           dimnames =MCMCiterZZ[[1]]$dimname[1:2] )
    istart=0
    for(c in nameschains){
      iend=istart+MCMCiterZZ[[c]]$dims[3]
      Output[[treatmentname]][["ZZallviews"]][,, (istart+1):iend]=PostMultipleChains[[c]][[treatmentname]]$ZZallviews
      istart=iend
    }
    rm(MCMCiterZZ, istart, iend,c)
    

    ###combine view allocations and checks on the structure 
    MCMCiterGammas=lapply(PostMultipleChains,function(x){
      list(dimname=dimnames(x[[treatmentname]]$Gammas),
           dims=dim(x[[treatmentname]]$Gammas))})
    
    if(!all(sapply(2:length(MCMCiterGammas), function(i) 
      identical(MCMCiterGammas[[i]]$dimname[1], 
                MCMCiterGammas[[1]]$dimname[1])))){
      stop("Dimesions and/or dimnames of view allocations  differ for chains: may indicate 2 different models")
    }else{
      if(  totMCMC!=sum(sapply(MCMCiterGammas, function(x) x$dims[2]))){
        stop("totMCMC iterations differ for beta and gammas: may indicate 2 different models")
        
      }
      
    }
    
    #combine the chains: gammas
    Output[[treatmentname]][["Gammas"]]=matrix(NA, nrow = MCMCiterGammas[[1]]$dims[1],
                                                           ncol = totMCMC, 
                                                  dimnames =MCMCiterGammas[[1]]$dimname[1] )
    istart=0
    for(c in nameschains){
      iend=istart+MCMCiterGammas[[c]]$dims[2]
      Output[[treatmentname]][["Gammas"]][, (istart+1):iend]=PostMultipleChains[[c]][[treatmentname]]$Gammas
      istart=iend
    }
    rm(MCMCiterGammas, istart, iend,c)
    
    
    
    ###combine precision parameter DPs and checks on the structure 
    MCMCiteralpha=lapply(PostMultipleChains,function(x){
      list(dimname=dimnames(x[[treatmentname]]$alpha),
           dims=dim(x[[treatmentname]]$alpha))})
    
    if(!all(sapply(2:length(MCMCiteralpha), function(i) 
      identical(MCMCiteralpha[[i]]$dimname[1], 
                MCMCiteralpha[[1]]$dimname[1])))){
      stop("Dimesions and/or dimnames of DPs precision param  differ for chains: may indicate 2 different models")
    }else{
      if(  totMCMC!=sum(sapply(MCMCiteralpha, function(x) x$dims[2]))){
        stop("totMCMC iterations differ for beta and gammas: may indicate 2 different models")
        
      }
      
    }
    
    #combine the chains: DPs hyperparameters
    Output[[treatmentname]][["alpha"]]=matrix(NA, nrow = MCMCiteralpha[[1]]$dims[1],
                                               ncol = totMCMC, 
                                               dimnames =MCMCiteralpha[[1]]$dimname[1] )
    istart=0
    for(c in nameschains){
      iend=istart+MCMCiteralpha[[c]]$dims[2]
      Output[[treatmentname]][["alpha"]][, (istart+1):iend]=PostMultipleChains[[c]][[treatmentname]]$alpha
      istart=iend
    }
    rm(MCMCiteralpha, istart, iend,c)
     
    
    ###combine theta unit specific checks on the structure 
    MCMCiterThetaUS=lapply(PostMultipleChains,function(x){
      list(dimname=dimnames(x[[treatmentname]]$Thetaus),
           dims=dim(x[[treatmentname]]$Thetaus))})
    
    if(!all(sapply(2:length(MCMCiterThetaUS), function(i) 
      identical(MCMCiterThetaUS[[i]]$dimname[1:2], 
                MCMCiterThetaUS[[1]]$dimname[1:2])))){
      stop("Dimesions and/or dimnames of theta unit specific differ for chains: may indicate 2 different models")
    }else{
      if(  totMCMC!=sum(sapply(MCMCiterThetaUS, function(x) x$dims[3]))){
        stop("totMCMC iterations differ for beta and theta unit specific: may indicate 2 different models")
        
      }
      
    }
    #combine the chains: Thetaus
    Output[[treatmentname]][["Thetaus"]]=array(NA, dim = c(MCMCiterThetaUS[[1]]$dims[1:2],totMCMC), 
                                                  dimnames =MCMCiterThetaUS[[1]]$dimname[1:2] )
    istart=0
    for(c in nameschains){
      iend=istart+MCMCiterThetaUS[[c]]$dims[3]
      Output[[treatmentname]][["Thetaus"]][,, (istart+1):iend]=PostMultipleChains[[c]][[treatmentname]]$Thetaus
      istart=iend
    }
    rm(MCMCiterThetaUS, istart, iend,c)
    
    
    
    ###combine posterior predictive outcome and checks on the structure 
    MCMCiterPostPredictiveY=lapply(PostMultipleChains,function(x){
      list(dimname=dimnames(x[[treatmentname]]$PostPredictiveY),
           dims=dim(x[[treatmentname]]$PostPredictiveY))})
    
    if(!all(sapply(2:length(MCMCiterPostPredictiveY), function(i) 
      identical(MCMCiterPostPredictiveY[[i]]$dimname[1], 
                MCMCiterPostPredictiveY[[1]]$dimname[1])))){
      stop("Dimesions and/or dimnames of posterior predictive  differ for chains: may indicate 2 different models")
    }else{
      if(  totMCMC!=sum(sapply(MCMCiterPostPredictiveY, function(x) x$dims[2]))){
        stop("totMCMC iterations differ for beta and posterior predictive: may indicate 2 different models")
        
      }
      
    }
    
    #combine the chains: PostPredictiveY
    Output[[treatmentname]][["PostPredictiveY"]]=matrix(NA, nrow = MCMCiterPostPredictiveY[[1]]$dims[1],
                                              ncol = totMCMC, 
                                              dimnames =MCMCiterPostPredictiveY[[1]]$dimname[1] )
    istart=0
    for(c in nameschains){
      iend=istart+MCMCiterPostPredictiveY[[c]]$dims[2]
      Output[[treatmentname]][["PostPredictiveY"]][, (istart+1):iend]=PostMultipleChains[[c]][[treatmentname]]$PostPredictiveY
      istart=iend
    }
    rm(MCMCiterPostPredictiveY, istart, iend,c)
    
  }
  
  
  #
  Output[["OriginalVariablesCoding"]]=PostMultipleChains[[1]]$OriginalVariablesCoding
  
  Output[["Info"]]=lapply(PostMultipleChains,function(x){x$Info})
  Output[["Info"]][["MCMCsamplesize"]]=totMCMC
  
  #names betas, ieprognostic variables effects
   
  return(Output)
}



#Function that select the variable whose MPP is in belonging view v is > tha cut off
SelectVariablesByMPP=function(PostDraws, maxCat=NULL,Cutoff=0.9, View=1){
  if(is.null(maxCat)){maxCat=max(PostDraws)}
  
  if(maxCat<0 || maxCat< max(PostDraws)){stop("maxCat must be greather than 0 and >= max(PostDraws)")}
  
  if(View==0 || View>= maxCat){stop("View must be greather than 0 and <= maxCat")}
  
  
  Result=(rowMeans(PostDraws==View)>= Cutoff)
  #print(table(Result))
  

  return(Result)
  
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

  if(TypeofPlot==2){
    A=ggplot(Result, aes(x = factor(Var1), y = value, fill = Var2)) +
      geom_bar(stat = "identity", width = 0.5) +  # narrower bars
      xlab("gamma_d") +
      ylab("Probability") +  scale_y_continuous(expand = c(0, 0)) +  theme_minimal()
  }


  if(TypeofPlot==1){

    A=ggplot(Result, aes(x = factor(Var1), y = value)) +
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

  plot(A)
  
}

#6. plot betas and prior
PlotBeta = function(BetasArr, Listhyperparam = list(mu = 0, sigma = 2.5, nu = 7),
                    cred_interval = 0.95) {
  require(tidyverse)
  require(ggplot2)
  require(reshape2)
  require(LaplacesDemon)
  
  # Set dimnames if not present
  if(is.null(dimnames(BetasArr))){
    dimnames(BetasArr) = list(
      paste0("p", 1:dim(BetasArr)[1]),
      paste0("r", 1:dim(BetasArr)[2]),
      paste0("mcmc", 1:dim(BetasArr)[3])
    )
  }
  
  # Melt array
  df = reshape2::melt(BetasArr) %>% rename(Progn = Var1, Category = Var2, mcmcindex = Var3)
  
  # Precompute probs for quantile
  probs = c((1 - cred_interval)/2, 1 - (1 - cred_interval)/2)
  
  # Compute credible intervals
  ci_df = df %>%
    group_by(Progn, Category) %>%
    summarize(
      lower = quantile(value, probs = probs[1]),
      upper = quantile(value, probs = probs[2]),
      .groups = "drop"
    )
  
  # Plot
  df %>%
    ggplot(aes(x = value)) +
    geom_density() +
    facet_grid(Progn ~ Category, scales = "free") +
    stat_function(fun = LaplacesDemon::dst, args = Listhyperparam,
                  color = "black", linetype = "dotted") +
    geom_segment(data = ci_df, 
                 aes(x = lower, xend = upper, y = 0, yend = 0), 
                 color = "black", size = 1) 
    
    }


#7. Estimated utility for treatment t
PosteriorUtility=function(PostPredictivetreatment_t, weights=c(0,40,100)){
  if(max(PostPredictivetreatment_t)+1!=length(weights)){stop("n category respone != than weights length")}
  Utility_hat=0
  for( k in 1:length(weights)){
    Utility_hat=Utility_hat+rowSums(PostPredictivetreatment_t==(k-1))*weights[k]
  }
return(Utility_hat/ncol(PostPredictivetreatment_t))
}







#Plot clustering and view with selected variables
PlotViewsandClusterAnalysis=function(Post =PostDraws,
                          view=1,
                          cutoff=0.9,
                          Outcomename="Outcome",
                          plotOutcome = TRUE, 
                          Dataset = dataset,
                   LevelsResp=c("Complete Response"=2,  
                     "Stable Disease"=1,
                     "Progressive Disease"=0),
                          UsePostPred=FALSE# True only if the Post predictive units are that of the sample
                   ){
  
  require(ggtern)
  listresults=list()
  

  

     for( t in names(Post)[grep("Treatment",names(Post))]){
      trtnumb=as.numeric(sub(".*?(\\d+)$", "\\1", t))
      
      listresultsTRT=list()
      
      
      #selected gammas 
      Indselvariables=SelectVariablesByMPP(Post[[paste0("Treatment",1)]]$Gammas,
                                           Cutoff = cutoff,
                                           View =  view) 
      
      NamesselVar=sapply(names(Indselvariables)[Indselvariables], 
                         function(x)strsplit(x,"gamma")[[1]][2] )
      
      
      DP=data.matrix(dataset[dataset$Treatment==trtnumb , NamesselVar])
      rownames(DP)=1:nrow(DP)
      
      #clustering structure
      Z=OptPostPart(Post[[t]]$ZZallviews)[,view]
      
      dfannotation=data.frame(Outcome=dataset[dataset$Treatment==trtnumb, Outcomename],
                              Cluster=as.factor(Z) )
      
      #compute proportion 
      tab=table(dfannotation$Cluster, dfannotation$Outcome)
      proptable=prop.table(tab, margin = 1)
      n_cluster = rowSums(tab)
      se_tab = sqrt(proptable * (1 - proptable) / n_cluster)

      
      fishertest=fisher.test(tab,simulate.p.value = TRUE)
 
      
      listresultsTRT[["Plotsview1"]]=
      
      
      
      listresultsTRT[["proptable"]]=cbind(proptable, n_cluster)
      listresultsTRT[["se_tab"]]=se_tab
      listresultsTRT[["fishertest"]]=fishertest
      
      
      ##Posterior predictive checks 
      if(UsePostPred){
        if(nrow(Dataset)!=nrow(Post[[t]]$PostPredictiveY)){stop("nrow(Dataset)!=nrow(Post[[t]]$PostPredictiveY")}
      }
      
      #mean predictive response 
      probitilde=apply(Post[[t]]$PostPredictiveY[dataset$Treatment==trtnumb,], 1, function(x){
        table(factor(x))/length(x)
      })
      
      
      #Predictive conditional on clusters
      probitilde=apply(Post[[t]]$Thetaus, c(2,3), softmax_ref3)
      dimnames(probitilde)[[1]]=c("PD","SD","CR")
      
      
      probitildetibble=reshape2::melt(probitilde) %>% as_tibble()  %>% 
          pivot_wider(names_from = Var1,  values_from = value  ) %>% 
          mutate(unit=as.integer(sub(".*i", "", Var2)))


      dfternplot=probitildetibble %>% group_by(Var3) %>% mutate(Cluster=Z) %>% ungroup() %>% 
                group_by(Cluster,Var3) %>% 
                summarise(mcmcind=mean(Var3), `PD`=mean(`PD`),
                          `SD`=mean(`SD`), `CR`=mean(`CR`))
                
      #ternary plots
      PostThetatrt=ggtern(data=dfternplot,
             aes(x=`PD`, y=`SD`, z=`CR`, 
                                       colour  = factor(Cluster),fill= (Cluster))) +
        geom_density_tern(bdl=0.001) +  
        theme_showarrows() +    #draw labeled arrows beside axes
        facet_wrap(Cluster~., ncol=5)+
        labs(title = NULL) +
        theme(
          tern.axis.text  = element_blank(),
          tern.axis.ticks = element_blank(),
          tern.axis.title = element_text(),
          strip.text = element_blank(),
          strip.background = element_blank()
        )+
        labs(color = "Cluster", fill = "Cluster")+  theme(
          tern.axis.text = element_blank(),
          tern.axis.title = element_blank()
        )+ggtitle(paste0("Treatment ", trtnumb))
      
      
      
        
      listresultsTRT[[paste0("ternplot")]]=PostThetatrt
        
        listresults[[t]]=listresultsTRT
        
        }
  
return(listresults)
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
ComputeESM=function(PredOptTrt, ActualAssTrt, TrueResp, 
                   responders=c("Complete Response")){
  
  yresp=(TrueResp%in%responders)*1
  

  #prop respondents under rnd all
  P_resp_rnd=mean(yresp) 

  #weights: prop recommended trt 
  passtrt1=mean(PredOptTrt==1)
  passtrt2=mean(PredOptTrt==2)
  
  #prop response in assigned and Recommended trt 
  p_rec_ass1=mean(yresp[(ActualAssTrt==1 & PredOptTrt==1) ]*1)

  
  p_rec_ass2=mean(yresp[(ActualAssTrt==2 & PredOptTrt==2) ])
  

  
  
  cat("increase/decrease in the response rate wrt rnd assignment: ", 
      (p_rec_ass1*passtrt1+p_rec_ass2*passtrt2)/P_resp_rnd, "\n" )
  
  ESM=(p_rec_ass1*passtrt1+p_rec_ass2*passtrt2)-P_resp_rnd
  
  return(ESM)
}

 



ComputeNPC=function(TrueResponse, PredResponse, LevelsResp ){
  if(!is.factor(TrueResponse)){stop("TrueResponse must be factor")}
  
  PredResponseFactor = factor(
    PredResponse,
    levels = LevelsResp,
    labels = names(LevelsResp)
  )
  
  return(mean(PredResponseFactor == TrueResponse))
}



 
#10. Plot estimated utility  (ONLY IF TRUE UTILITY AVAILABLE)
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



#10. Plot estimated utility (ONLY IF TRUE UTILITY AVAILABLE)
PlotTrueveEstdiffUtilityMethods=function(Dataset, namesutility=c("utility1","utility2"), PredUtilitylist){
  require(tidyverse)
  diffU1U2true=Dataset[,namesutility[1]]-Dataset[,namesutility[2]]

  #
  order=order(diffU1U2true)

  #RMSE on trt utility
  lapply(PredUtilitylist, function(PredutMethod)
    {print=apply(Dataset[,c(namesutility[1],namesutility[2])]-PredutMethod, 2,
                     function(x){sqrt(mean(x^2))})
    
    cat("RMSE ", print)}
    
    )
  
  

  diffU1U2pred=lapply(PredUtilitylist, function(l){l[,1]-l[,2]})

  #reordering and melt
  diffU1U2pred= reshape2::melt(lapply(diffU1U2pred, function(l){l[order]}))

  
  #RMSE on trt utility
  cbind(diffU1U2pred, valuetrue=c(diffU1U2true[order],diffU1U2true[order])) %>% 
    mutate(x2=(value-valuetrue)^2) %>% group_by(L1) %>% 
    summarise(RMSE=sqrt(mean(x2))) %>% print()
  

  tibble(DiffUtility=c(diffU1U2true[order],diffU1U2pred[,1]), reorderedindex=c(1:length(diffU1U2true),
                                                                               rep(1:length(diffU1U2true),length(unique(diffU1U2pred[,2])))),
         Legend=c(rep("True", length(diffU1U2true)), diffU1U2pred[,2])) %>%
    ggplot(., aes(x=reorderedindex, y=DiffUtility,color=Legend))+
    geom_line()
}






##### Diagnostic only for betas and alphas
DiagnosticBetaAlpha=function(Post, iter=NULL, 
                             pathtosaveplots="~/"){
  require(coda)
  require(ggmcmc)
  
  #Multiple chains? 
  nch=sum(grepl("Chain",names(Post)))
  
  
  if(nch>1){
    #check is chains have same length    
    Chlen=sapply(Post, function(CH){dim(CH$Beta)[3]})
    if (length(unique(Chlen)) != 1) {
      stop( "Chains of different length")
    }else{
      if(is.null(iter)){
        iter=seq_len(unique(Chlen))
      }}
  }else{
    if(is.null(iter)){
      iter=seq_len(dim(Post[[1]]$Beta)[3])
    }
  }
  
  
  #betas   
  Betacodalist=coda::mcmc.list( 
    lapply(Post, function(CH){
      CH=CH$Beta[,,iter]
      #store names 
      Variablesnames = dimnames(CH)[[1]]
      refindex = dimnames(CH)[[2]]
      
      CH=matrix(
        aperm(CH, c(3, 1, 2)),              # s × (beta × r)
        nrow = dim(CH)[3],
        ncol = dim(CH)[1] * dim(CH)[2]  )
      
      colnames(CH) = as.vector(outer(Variablesnames, refindex, paste, sep = "_"))
      
      return(coda::mcmc(CH))
    }))
  
  Betaggs=ggmcmc::ggs(Betacodalist)
  rm(Betacodalist)
  
  
  #alphas
  AlphasDPs=coda::mcmc.list( 
    lapply(Post, function(CH){
      Trtnames=names(CH)[grepl("Treatment",names(CH))]
      CH=lapply(CH[Trtnames], function(CHtr) t(CHtr$alpha[,iter]))
      
      CH = do.call(
        cbind,
        lapply(names(CH), function(trt) {
          M = CH[[trt]]
          colnames(M) = paste(trt, colnames(M), sep = " ")
          M
        })
      )
      return(coda::mcmc(CH))
    }))
  
  Alphasggs=ggmcmc::ggs(AlphasDPs)
  rm(AlphasDPs)
  
  
  #compute geweke beta alphas
  FinalDiagn=ggmcmc::ggs_geweke(bind_rows(Betaggs, Alphasggs),plot = FALSE) |> 
    dplyr::select(Parameter, Chain, z) |>
    tidyr::pivot_wider(
      names_from  = Chain,
      values_from = z,
      names_prefix = "z_chain"
    )
  
  #compute Rhat and Effective if nchain>1
  if(nch>1){
    RhatEss=ggmcmc::ggs_diagnostics(bind_rows(Betaggs, Alphasggs),
                                    version_rhat = "BG98") |>
      dplyr::filter(Diagnostic %in% c("Rhat", "Effective")) |>
      dplyr::select(Parameter, Diagnostic, value) |>
      tidyr::pivot_wider(
        names_from  = Diagnostic,
        values_from = value )
    
    FinalDiagn=dplyr::left_join(FinalDiagn,RhatEss)
  }
  
  
  
  
  #save plots pdf
  if(!is.null(pathtosaveplots)){
    #plots betas
    ggmcmc::ggmcmc(Betaggs, plot=c("density", "ggs_traceplot","ggs_compare_partial", "ggs_autocorrelation",
                           "ggs_geweke","ggs_Rhat","ggs_effective","ggs_caterpillar"),
           file = paste0(pathtosaveplots,"/Betas.pdf"))
    
    #plots betas
    ggmcmc::ggmcmc(Alphasggs, plot=c("density","ggs_traceplot", "ggs_compare_partial", "ggs_autocorrelation",
                             "ggs_geweke","ggs_Rhat","ggs_effective"),
           file = paste0(pathtosaveplots,"/Alphas.pdf"))
    
  }
  
  
  return(FinalDiagn)
  
  
}




#Function to plot the heatmap of clustering for each view and treatment 
PlotPSM=function(Post, Dataset=NULL, plotOutcome=TRUE, DrawOnePlot=FALSE,main_title=NULL,
                 nameTreatment="Treatment", iters=NULL,
                 nameOutcome="Outcome"){
  require(ComplexHeatmap)
  listplot=list()
  
  if(!is.null(iters)){Post=KeepIters(Post,iters)}
  
  
  for (trt in names(Post)[grep( "Treatment", names(Post))]) {
    #Number of non-null views is the same in all treatments
    nNNV=dim(Post$Treatment1$ZZallviews)[2]
    
    #If plotOutcome=TRUE both the predicted and oberved outcome will be plotted as annotation
    #Note to plot the predicted outcome the predicted dataser in the tmvbpr but be NULL, so
    #that the  Post$Treatment[X]$PostPredictiveY refers to the sample (in sample prediction)
    if(plotOutcome){
      nubtrt=as.numeric(strsplit(trt,"Treatment")[[1]][2])
      indextrt=(Dataset[[nameTreatment]]==nubtrt)
      ObsOutcome=Dataset[[nameOutcome]][indextrt]
      
      #MAP Outcome    
      PredOutcome=OptPostCategorical(Post[[trt]]$PostPredictiveY) #note this numerical coded
      
      PredOutcome=PredOutcome[indextrt]
      
      PredOutcome=names(Post$OriginalVariablesCoding$Outcome)[match(PredOutcome, Post$OriginalVariablesCoding$Outcome)]
      
      
      PredOutcome=factor(PredOutcome, levels = names(Post$OriginalVariablesCoding$Outcome))
    }
    
    listtrt=list()
    for (v in 1:1) {
      PSM=mcclust::comp.psm(t(Post[[trt]]$ZZallviews[,v,]) )
      dimnames(PSM)=list(1:dim(PSM)[1], 1:dim(PSM)[2])
      
      #compute optimal partition 
      Cluster=mcclust::maxpear(PSM)$cl
      Clusterfct=factor(Cluster)
      
      
      # Create histogram annotation
      if(plotOutcome){
        
        outcome_colors = setNames(
          scales::hue_pal()(length(levels(PredOutcome))),
          levels(PredOutcome)
        )
      }
      
      #colors clusters
      if(nlevels(Clusterfct)<8 && nlevels(Clusterfct)>2 ){
        print(nlevels(Clusterfct))
        temp=RColorBrewer::brewer.pal(n = nlevels(Clusterfct), 
                                      name = "Set1")
      }else{
        temp=scales::hue_pal()(nlevels(Clusterfct))
      }
      names(temp)=levels(Clusterfct)
      
      
      
      ha_top = HeatmapAnnotation(
        Cluster = Clusterfct,
        which = "column",
        col=list(Cluster=temp)
      )
      
      if(plotOutcome){
        
        ha_left = HeatmapAnnotation(df = data.frame(
          `Observedoutcome`=as.character(ObsOutcome),
          `Predictedoutcome`=PredOutcome), 
          which = "row",
          col =list(
            `Observedoutcome`  = outcome_colors,
            `Predictedoutcome` = outcome_colors
          )
        )
      }else{
        ha_left=NULL
        
      }
      #Heatmap
      HM=Heatmap(
        PSM,
        name = "PSM", 
        show_row_dend = FALSE,
        show_column_dend = FALSE,
        # histogram = anno_histogram(ObsOutcome_ordered)
        top_annotation  = ha_top,
        left_annotation = ha_left,
        column_title = paste0(trt, "  View ", v),
        column_title_gp = grid::gpar(
          fontsize = 18,
          fontface = "bold"
        ),
        col =viridis::viridis(100,begin = 0,end = 1, direction = 1)   #circlize::colorRamp2(c(0, 1), c("white", "#008080"))
      ) 
      
      listtrt[[paste0("view",v)]]=HM
    }
    
    listplot[[trt]]=listtrt
    
  }
  
  
  if(DrawOnePlot){
    library(ComplexHeatmap)
    library(grid)
    library(gridExtra)
    ht_list <- unlist(listplot, recursive = FALSE)
    
    ncol=length(listplot[[1]])
    
    
    
    grobs = lapply(ht_list, function(ht) {
      grid::grid.grabExpr(draw(ht))
    })
    
    gridExtra::grid.arrange(
      grobs = grobs,
      ncol = ncol,
      top = if (!is.null(main_title)) {
        grid::textGrob(
          main_title,
          gp = grid::gpar(fontsize = 26, fontface = "bold")
        )
      } else NULL
    )
    
    
  }else{
    listplot}
}
 




###Diagnostic variables in relevant view
DiagnosticViewAllocations=function(PostMultipleChains, iter=30000:60000){
  if(!is.null(iter)){
    PostMultipleChains=lapply(PostMultipleChains, function(x)KeepIters(x,iter))
  }
  
  require(mltools)
  

  MPList=lapply(PostMultipleChains, function(CH){
    Trtnames=names(CH)[grepl("Treatment",names(CH))]
    
    Margprobinclusion=lapply(CH[Trtnames], function(CHtr) {
      CHtr=CHtr$Gammas
      
      nobs=dim(CHtr)[1]
      nmcmc=dim(CHtr)[2]

      
      Result=matrix(NA,nrow = nobs,ncol =1 )
      #view to inspect (not due to label switching across irrelevant views only for k=1 
      #the following diagnostic measure make sense )
      k=1
      colnames(Result)=paste0("View",1)
      rownames(Result)=rownames(CHtr)
      #for (k in 0:(Nview-1)) {
        
        Result[,k]=rowSums(CHtr==k)
      #}
      Result=Result/nmcmc
      
      return(Result)              
    })
    
    return(Margprobinclusion)
    
  })
  
  chains = names(MPList)
  treatments = names(MPList[[1]])

 
  
  ##if nchains==2
  if(length(chains)==2){
    
  for(trt in treatments){
    ch1 = (PostMultipleChains$Chain1[[trt]]$Gammas == 1) * 1
    ch2 = (PostMultipleChains$Chain2[[trt]]$Gammas == 1) * 1
    
    gv1 = array(c(ch1, ch2), 
                 dim = c(nrow(ch1), ncol(ch1), 2))
    
    AgreementyVariables=apply(gv1, 1, function(x){
      mean(x[1,] == x[2,])
      
    })

    
    # 
    cat("Treatment ", trt, " agreement between chains \n")
    print(summary(AgreementyVariables))
    
    
    
  }
  }
  
  
  
  
  #to plot
  RESMP = lapply(treatments, function(tr) {
    arr=simplify2array(lapply(chains, function(ch) MPList[[ch]][[tr]]))
    arr
  })
  names(RESMP) = treatments

  
  max_diff = lapply(RESMP, function(arr) {
                    apply(arr, c(1, 2), function(x) max(x) - min(x))})

  
  cat("\n summary statistics of max difference Marginal probabilities: \n ")
  print(lapply(max_diff,summary))
  
  
  df_all <- do.call(rbind, lapply(names(max_diff), function(tr) {
    tmp <- reshape2::melt(max_diff[[tr]])
    colnames(tmp) <- c("Variable", "State", "MaxDiff")
    tmp$Treatment <- tr
    tmp
  }))
  ggplot(df_all, aes(x = State, y = Variable, fill = MaxDiff)) +
    geom_tile() +
    scale_fill_gradientn(
      colors = c("white", "yellow", "orange", "red"),
      values = scales::rescale(c(0, 0.02, 0.05, 0.10, max(df_all$MaxDiff))),
      limits = c(0, max(df_all$MaxDiff)),
      name = "MaxDiff"
    ) +
    facet_wrap(~Treatment, scales = "free_y") +
    theme_minimal() +
    theme(axis.text.y = element_text(size = 6)) +  # adjust for many variables
    labs(title = "Max cross-chain differences by Treatment")
  
  }

###Diagnostic Cluster allocations
DiagnosticClusterallocations=function(PostMultipleChains, iter=NULL){
  require(ggplot2)
  require(tidyverse)
  
  if(!is.null(iter)){
    PostMultipleChains=lapply(PostMultipleChains, KeepIters)
  }
  
  
  NumberofClust=
    lapply(PostMultipleChains, function(CH){
      Trtnames=names(CH)[grepl("Treatment",names(CH))]
      
      nclustrbytrts=lapply(CH[Trtnames], function(CHtr) {
        CHtr=CHtr$ZZallviews
        dimnames(CHtr)=NULL
        ncbytrt=t(apply(CHtr,c(2,3), function(x) length(unique(x)) ) )
        # rownames(ncbytrt)=NULL
        return(ncbytrt)              
      })
      
      return(nclustrbytrts)
      
    })
  
  
  Nclustermelted=reshape2::melt(NumberofClust)
  colnames(Nclustermelted)=c("iter","View","Nclusters","Treatment","Chain")
  Nclustermelted[,"View"]=paste0("View", Nclustermelted[,"View"])
  
  Nclustermelted=Nclustermelted %>% filter(View=="View1")
  
  #trace plot n clusters each treatment view
  NcluPlot=ggplot2::ggplot(Nclustermelted, ggplot2::aes(x=iter, y=Nclusters, color=Chain))+
    ggplot2::geom_line(linewidth=1,alpha=0.7)+
    ggplot2::facet_grid(Treatment~View)
  
  print(NcluPlot)
  
  PostmeanNclust=Nclustermelted %>% dplyr::group_by(Treatment, View, Chain) %>% 
    dplyr::summarise(Mean=mean(Nclusters) ) %>% 
    tidyr::pivot_wider(
      names_from  = Chain,
      values_from = Mean )
  print(PostmeanNclust)
  ##PSM multiple chains
  PSMmultChains=lapply(PostMultipleChains, function(x) PlotPSM(Post =x ,
                                                               plotOutcome =FALSE ,
                                                               DrawOnePlot = FALSE) )
  
  chainsnames=names(PSMmultChains)
  names(chainsnames)=names(PSMmultChains)
  
  PSMmultChains=lapply(chainsnames, function(chname){
    
    PSMslistChain = unlist(PSMmultChains[[chname]], recursive = FALSE)
    PSMslistChain=lapply(PSMslistChain, FUN = function(x, name=chname){
      x@column_title=name
      x
    })
  })
  
  
  #concatenation across chains
  
  chains = names(PSMmultChains)
  keys   = names(PSMmultChains[[1]])   # Treatment1.view1, ...
  
  PSMmultChainsRES = setNames(vector("list", length(keys)), keys)
  
  for (k in keys) {
    PSMmultChainsRES[[k]] <- Reduce(
      `+`,
      lapply(chains, function(ch) PSMmultChains[[ch]][[k]])
    )
  }
  
  library(grid)
  library(gridExtra)
  library(ComplexHeatmap)
  #PSM diagnostic plot
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(
    nrow = length(PSMmultChainsRES), ncol = 1
  )))
  
  for (i in seq_along(PSMmultChainsRES)) {
    pushViewport(viewport(layout.pos.row = i))
    draw(
      PSMmultChainsRES[[i]],
      column_title = names(PSMmultChainsRES)[i],
      newpage = FALSE
    )
    popViewport()
  }
  
}




















































##Compute posterior Y for each treatment STAN (cambiare script non salvare qui)
postTrtinfStan=function(PostY=A$ypred, weights=c(0,40,100)){
  n=nrow(PostY)
  #map yt1 and yt2
  Yt1=OptPostCategorical(PostY[1:(n/2),])
  Yt2=OptPostCategorical(PostY[(n/2+1):n,])

  #estimated utility
  u1=PosteriorUtility(PostY[1:(n/2),],weights )
  u2=PosteriorUtility(PostY[(n/2+1):n,], weights)

  #opt trearment
  Otrt=apply(cbind(u1,u2), 1, which.max)


  #Outcome: y1 if utility 1> unility 2 ,ow y2
  predoutc=Yt1*(Otrt==1)+Yt2*(Otrt==2)

  return(list(y1pred=Yt1, y2pred=Yt2, PredOptTrt=Otrt, PredOutcome=predoutc, Utility=cbind(u1,u2)))
}














###Function only for dm
CreateDatalistDM=function(dataset, predictiveDataset=NULL, 
                Outcomemapping=NULL, 
                NumberofViews=NULL,
                ListPrior=NULL,
                ListInitialValues=NULL,
                nameTreatment="Treatment", 
                nameOutcome="Outcome",
                namesPredictive=DATA$namesCNV,
                namesPrognostic=c("Age", "Morphology","Sex", "Laterality"),
                MCMCfinalSampleSize=NULL,
                MCMCburnin=NULL,
                MCMCThinning=NULL,
                ListMCMCparam=NULL,
                HowgammasInit=NULL
                #ListInit Aggiungere in futuro
){
  #####
  
  listRelabeling=list()
  if(class(dataset)!="data.frame"){stop("dataset must be a data.frame")}
  
  if(!all(c(nameOutcome,namesPrognostic, namesPredictive) %in% colnames(dataset))){
    stop("dataset lacks one or more required columns")}
  
  
  
  #relabel outcome (convert fct in numeric from 0:nlev(outcome)-1)
  Relabtemp=FcttoNum(dataset[[nameOutcome]],Mapping = Outcomemapping)
  dataset[[nameOutcome]]=Relabtemp$Vect
  listRelabeling[[nameOutcome]]=Relabtemp$Mappingused
  
  
  #set n of categories outcome
  ncategoriesOutcome=length(Relabtemp$Mappingused)
  
  #relabel predictive variables
  ncategoriesPred=setNames(rep(NA,length(namesPredictive)), namesPredictive)
  for (predvarnam in namesPredictive) {
    
    Relabtemp=FcttoNum(dataset[[predvarnam]])
    dataset[[predvarnam]]=Relabtemp$Vect
    listRelabeling[[predvarnam]]=Relabtemp$Mappingused
    
    ncategoriesPred[predvarnam]= length(Relabtemp$Mappingused)
    
    
  }
  
  
  #Build dummies for factor prognostivariables
  oldnamesPrognostic=namesPrognostic
  #prognostic variables
  contrastsarg = lapply(
    dataset[oldnamesPrognostic],
    function(x) if (is.factor(x)) contrasts(x)
  )
  contrastsarg=contrastsarg[ !unlist(lapply(contrastsarg, is.null)) ]
  
  
  ModmatProgn=model.matrix(~.,dataset[oldnamesPrognostic], 
                           contrasts.arg =contrastsarg ) 
  Prognvarrecoded = as.data.frame( ModmatProgn[, setdiff(colnames(ModmatProgn),
                                                         "(Intercept)"), 
                                               drop = FALSE]  )
  #reference levels
  listRelabeling[["RefPrognosticfct"]]=lapply(dataset[oldnamesPrognostic], function(x) if (is.factor(x)) levels(x)[1] else NULL)
  
  dataset[namesPrognostic]=NULL
  dataset[colnames(Prognvarrecoded)]=Prognvarrecoded
  
  #change the names of prognostic variables with those recoded
  namesPrognostic=colnames(Prognvarrecoded)
  
  
  
  #Predictive dataset
  if (!is.null(predictiveDataset)) {
    if(!all(c(oldnamesPrognostic, namesPredictive) %in% colnames(predictiveDataset))) {
      stop("predictiveDataset lacks required columns: namesPrognostic, namesPredictive must be the same in predictiveDataset and dataset")
    }
    
    #relabel predictive variables 
    for (predvarnam in namesPredictive) {
      mapping = listRelabeling[[predvarnam]]
      
      predictiveDataset[[predvarnam]] =
        as.integer(mapping[as.character(predictiveDataset[[predvarnam]]) ])
      
      #check if there are more levels in predictive dataset 
      if(any(is.na(predictiveDataset[[predvarnam]]))){
        stop(paste0("In predictiveDataset the variable ", predvarnam, "has more levels than in dataset" ))
        
      }
      
    }
    
    #prognostic variables
    ModmatPredProgn=model.matrix(~., predictiveDataset[oldnamesPrognostic],
                                 contrasts.arg = contrastsarg)
    
    PrognvarrecodedPred = as.data.frame( ModmatPredProgn[, setdiff(colnames(ModmatPredProgn),
                                                                   "(Intercept)"), 
                                                         drop = FALSE]  )
    predictiveDataset[oldnamesPrognostic]=NULL
    predictiveDataset[colnames(PrognvarrecodedPred)]=PrognvarrecodedPred
    
  }
  
  
  
  
  ##
  rm(ModmatProgn,Prognvarrecoded, Relabtemp, oldnamesPrognostic, ModmatPredProgn,PrognvarrecodedPred )
  
  #create suitable datastructure for the cpp function
  Datalist=CreateDataStructureforCpp(Dataset = dataset, DatasetPostPred=predictiveDataset,
                                     nameTreat=nameTreatment,
                                     nameY=nameOutcome,
                                     namesPred =namesPredictive,
                                     namesProg =namesPrognostic )
 
  return(Datalist) 
}

































































