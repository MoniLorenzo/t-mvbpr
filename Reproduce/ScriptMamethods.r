###This code is an adaptation of the code provided by Predone et al., 2024  
#to preform the methods described in Ma et al. 2019 in case of two treatment 

#@author: Lorenzo Moni, Silvia Liverani, Alberto Cassese, Francesco Claudio Stingo  

################################ Functions ########################################
library(ConsensusClusterPlus);
library("mvtnorm");
library(parallel)
library(doParallel)
library(doRNG)



##W
mymultt <- function(Xtrain=trtc1[,-1], X.pred=mycovX){
  #prognostic variables
  myln <- length(Xtrain[,1,drop=FALSE])
  myls <- Xtrain
  mylmu <- apply(myls, 2, mean)
  mymun <- myln/(myln+kappa0)*mylmu
  myS <- cov(myls)*(myln-1)
  
  ## for the covariates
  myd <- length(Xtrain[1,,drop=FALSE])       ## numer of covariates
  nu0 <- length(Xtrain[1,,drop=FALSE])+1       ## numer of covariates +1;
  lambda0 <- diag(myd)                ## identity matrix
  
  kappan <- kappa0 + myln
  nun <- nu0 + myln
  lambdn <- lambda0 + myS + kappa0*myln/(kappa0 + myln)*(mylmu - mu0)%*%t(mylmu - mu0)
  return2 <- mvtnorm::dmvt(x = X.pred, sigma=(kappan + 1)*lambdn/(kappan*(nun - myd + 1)), df = nun - myd + 1, log = FALSE)
}



#second step
con.cluster <- function(cons=con_clu[[2]][["consensusMatrix"]],#similarity matrix cons clustering
                        ymat = s_train_ymat,
                        yvec = s_train_yord,
                        prog = s_train_prog, trt = s_train_trt){
  ni <- dim(cons)[1]
  utpred1 <- matrix(1, nrow = ni, ncol=12)  ## ut1,ut2,trt
  
  for (i in 1:ni){
    totut <- ymat*cons[i,]  ### utility
    totpre <- totut[-i,]
    trtpre <- trt[-i]
    
    l_train_prog <- prog[-i,, drop=FALSE]
    l_yvec <- yvec[-i];
    
    # lik y reweighted using similaity 3cols; trt ; y ordinal; prognostic
    mytempt <- cbind(totpre, trtpre, l_yvec, l_train_prog)
    
    ## for genes
    #lik y reweighted using similaity 3cols subset by trt
    trt1g <- subset(mytempt[,1:3], mytempt[,4] == 0)
    trt2g <- subset(mytempt[,1:3], mytempt[,4] == 1)
    
    #sum col lik+prior
    alphatrt1 <- apply(trt1g, 2, sum) + prior1
    probpre1_temp <- alphatrt1/sum(alphatrt1)
    alphatrt2 <- apply(trt2g, 2, sum) + prior2
    probpre2_temp <- alphatrt2/sum(alphatrt2)
    #These probs are due to predictive
    
    #prognostic unit i
    mycovX <- prog[i,,drop=FALSE];                  ## observed covariates
    
    trtc <- subset(mytempt[,-(1:4)])
    
    trtc0 <- subset(trtc,trtc[,1] == 0)
    trtc0dst <- mymultt(Xtrain = trtc0[,-1, drop=FALSE],X.pred =  mycovX)
    
    trtc1 <- subset(trtc,trtc[,1] == 1)
    trtc1dst <- mymultt(trtc1[,-1, drop=FALSE], mycovX)
    trtc2 <- subset(trtc, trtc[,1] == 2)
    
    trtc2dst <- mymultt(trtc2[,-1, drop=FALSE], mycovX)
    trtpfy <- c(trtc0dst, trtc1dst, trtc2dst)   ##prob with cov
    
    #combine prognostic with predictive
    probpre1 <- probpre1_temp*trtpfy/sum(probpre1_temp*trtpfy)
    probpre2 <- probpre2_temp*trtpfy/sum(probpre2_temp*trtpfy)
    
    ## calculate the utility
    ut1pre <- probpre1%*%wk
    ut2pre <- probpre2%*%wk
    utpred1[i, 1] <- ut1pre
    utpred1[i, 2] <- ut2pre
    if(ut2pre > ut1pre){
      utpred1[i, 3] = 2
    }
    utpred1[i, 4:6] <- probpre1
    utpred1[i, 7:9] <- probpre2
    utpred1[i, 10:12] <- trtpfy/(sum(trtpfy))
    ut1pre<-ut2pre<-totut<-totpre<-mytempt<-trt1g<-trt2g<-alphatrt1<-alphatrt2<-NULL;
  }
  return <- utpred1;
}

PreUt<-function(mth, trt, out.response, SUB.ID){
  #mth<-mth;
  myresults <- cbind(mth, trt+1, out.response, SUB.ID)
  pred1 <- subset(myresults, myresults[,3] == 1)
  table1 <- table(pred1[,14], pred1[,13])
  pred2 <- subset(myresults, myresults[,3] == 2)
  table2 <- table(pred2[,14], pred2[,13])
  p1 <- sum(table1)/(sum(table1) + sum(table2))
  p2 <- sum(table2)/(sum(table1) + sum(table2))
  ## set prob as 0 for cases that none selected patients had response (all 0s);
  ## set prob as 1 for cases that none selected patients had non-response (all 1s);
  if(length(table1) == 4){
    crt1<-table1[2,1]/sum(table1[,1])
  }
  if(length(table1) < 4){
    crt1 <- 0
  }
  #### if table 1 is empty then set crt1=0;
  if(length(table2) == 4){
    crt2 <- table2[2,2]/sum(table2[,2])
  }
  if(length(table2) < 4){
    crt2 <- 0
  }
  
  #### summary meaures
  return <- crt1*p1 + crt2*p2 - sum(out.response)/length(out.response)
}


###################################################


##





#this function need to be seq since inside alredy have a parallelization
Mafunction=function(Datalistrep, Dpred=DP, method="hc",CORES=10){
  source("~/Sim/RFunctionsMultiTreatMVBPR.R")
  
  namerep=Datalistrep$repname
  
  DATA=Datalistrep$Data[,]
  
  ld=MultiTreatmentMVBPR(DATA, DatasetPostPred=DP,
                         nameTreat="asstrt",
                         nameY="outcome")
  
  #n prgognostic vars
  Qprogn=ncol(ld$fitdata[[1]]$WProgn)
  
  #setup Parameters:
  wk <- c(0,40,100)
  prior1 <- c(1/3,1/3,1/3)
  prior2 <- c(1/3,1/3,1/3)
  kappa0 <- 1
  mu0 <- rep(0, Qprogn)#c(0,0) #length = # prognostic var
  n <- length(ld$fitdata[[1]]$Y)+length(ld$fitdata[[2]]$Y)#124 #n train
  #  K <- 50 #replications
  
  
  #ADDItIONAL STEP tO DETERMINE THE NUMBER OF CLUSTERS FOR EACH UNIT USING LOOCV
  #  HC.sum.all <- foreach(k = 1:K) %dorng%
  
  DATA=ld
  
  train_pred <- rbind(ld$fitdata[[1]]$XDiscPred, ld$fitdata[[2]]$XDiscPred)
  #train_pred <- simdata$pred[[k]][1:124,]
  train_prog <- rbind(ld$fitdata[[1]]$WProgn, ld$fitdata[[2]]$WProgn)
  #train_prog <- simdata$prog[[k]][1:124,]
  
  train_yord <- c(ld$fitdata[[1]]$Y,ld$fitdata[[2]]$Y)
  
  #train_ymat <- simdata$ymat[[k]][1:124,]
  train_ymat=matrix(0,ncol = max(ld$fitdata[[1]]$Y)+1,
                    nrow =length(ld$fitdata[[1]]$Y)+length(ld$fitdata[[2]]$Y) )
  
  train_ymat[cbind(1:length(train_yord), train_yord+1)]=1
  
  train_trt <- c(rep(0, length(ld$fitdata[[1]]$Y)), rep(1, length(ld$fitdata[[2]]$Y)))
  
  
  
  #  HC.sum<-matrix(0,nrow=n,ncol=14)
  
  
  n=nrow(train_pred)
  
  ttt=Sys.time()
  # for (mysub in 1:n){
  set.seed(126)
  
  cl <- makeCluster(CORES)  # leave 1 core free
  registerDoParallel(cl)
  #   clusterExport(cl, c("mymultt", "con.cluster", "PreUt"))
  
  #"mymultt","con.cluster","PreUt"
  
  HC.sum <- foreach(mysub = 1:n, .combine=rbind,
                    .packages=c("cluster","ConsensusClusterPlus"),
                    .export=c("mymultt", "con.cluster", "PreUt",
                              "wk", "prior1", "prior2", "kappa0", "mu0")
  ) %dorng% {
    
    #FOR CATEGORICAL DATA
    s_train_pred= as.data.frame(train_pred[-mysub,])
    
    s_train_pred[] = lapply(s_train_pred, factor)
    
    # Compute Gower distance for categorical data
    Dist <- cluster::daisy(s_train_pred, metric = "gower")
    
    
    s_train_prog <- train_prog[-mysub,, drop=FALSE]
    s_train_yord <- train_yord[-mysub]
    s_train_ymat <- train_ymat[-mysub,]
    s_train_trt <- train_trt[-mysub]
    
    ### clustering using CONSENSUS MATRIX method ###################################
    require(ConsensusClusterPlus)
    con_clu <- ConsensusClusterPlus(Dist,maxK=15,pFeature=1,
                                    clusterAlg=method,#distance=Dist[,],
                                    #clusterAlg="km",distance="euclidean",
                                    #clusterAlg="pam",distance="manhattan",
                                    seed=126,)
    
    
    hc2 <- con.cluster(cons = con_clu[[2]][["consensusMatrix"]],
                       ymat = s_train_ymat,  yvec = s_train_yord,
                       prog = s_train_prog, trt = s_train_trt)
    hc3 <- con.cluster(con_clu[[3]][["consensusMatrix"]], ymat = s_train_ymat,  yvec = s_train_yord,
                       prog = s_train_prog, trt = s_train_trt)
    hc4 <- con.cluster(con_clu[[4]][["consensusMatrix"]], ymat = s_train_ymat,  yvec = s_train_yord,
                       prog = s_train_prog, trt = s_train_trt)
    hc5 <- con.cluster(con_clu[[5]][["consensusMatrix"]], ymat = s_train_ymat,  yvec = s_train_yord,
                       prog = s_train_prog, trt = s_train_trt)
    hc6 <- con.cluster(con_clu[[6]][["consensusMatrix"]], ymat = s_train_ymat,  yvec = s_train_yord,
                       prog = s_train_prog, trt = s_train_trt)
    hc7 <- con.cluster(con_clu[[7]][["consensusMatrix"]], ymat = s_train_ymat,  yvec = s_train_yord,
                       prog = s_train_prog, trt = s_train_trt)
    hc8 <- con.cluster(con_clu[[8]][["consensusMatrix"]], ymat = s_train_ymat,  yvec = s_train_yord,
                       prog = s_train_prog, trt = s_train_trt)
    hc9 <- con.cluster(con_clu[[9]][["consensusMatrix"]], ymat = s_train_ymat,  yvec = s_train_yord,
                       prog = s_train_prog, trt = s_train_trt)
    hc10 <- con.cluster(con_clu[[10]][["consensusMatrix"]], ymat = s_train_ymat,  yvec = s_train_yord,
                        prog = s_train_prog, trt = s_train_trt)
    hc11 <- con.cluster(con_clu[[11]][["consensusMatrix"]], ymat = s_train_ymat,  yvec = s_train_yord,
                        prog = s_train_prog, trt = s_train_trt)
    hc12 <- con.cluster(con_clu[[12]][["consensusMatrix"]], ymat = s_train_ymat,  yvec = s_train_yord,
                        prog = s_train_prog, trt = s_train_trt)
    hc13 <- con.cluster(con_clu[[13]][["consensusMatrix"]], ymat = s_train_ymat,  yvec = s_train_yord,
                        prog = s_train_prog, trt = s_train_trt)
    hc14 <- con.cluster(con_clu[[14]][["consensusMatrix"]], ymat = s_train_ymat,  yvec = s_train_yord,
                        prog = s_train_prog, trt = s_train_trt)
    hc15 <- con.cluster(con_clu[[15]][["consensusMatrix"]], ymat = s_train_ymat,  yvec = s_train_yord,
                        prog = s_train_prog, trt = s_train_trt)
    
    out.response <- as.numeric(s_train_yord > 1)
    SUB.ID <- c(1:(n-1))
    
    #  HC.sum[mysub,] =
    res= 
      c(PreUt(hc2, s_train_trt, out.response, SUB.ID),
        PreUt(hc3, s_train_trt, out.response, SUB.ID),
        PreUt(hc4, s_train_trt, out.response, SUB.ID),
        PreUt(hc5, s_train_trt, out.response, SUB.ID),
        PreUt(hc6, s_train_trt, out.response, SUB.ID),
        PreUt(hc7, s_train_trt, out.response, SUB.ID),
        PreUt(hc8, s_train_trt, out.response, SUB.ID),
        PreUt(hc9, s_train_trt, out.response, SUB.ID),
        PreUt(hc10, s_train_trt, out.response, SUB.ID),
        PreUt(hc11, s_train_trt, out.response, SUB.ID),
        PreUt(hc12, s_train_trt, out.response, SUB.ID),
        PreUt(hc13, s_train_trt, out.response, SUB.ID),
        PreUt(hc14, s_train_trt, out.response, SUB.ID),
        PreUt(hc15, s_train_trt, out.response, SUB.ID))
    return(res)
    
  }
  
  stopCluster(cl)
  ttt-Sys.time()
  K=1
  HC.sum.all <- array(unlist(HC.sum), dim = c(n, 14, K))
  
  
  objtokeep=c("HC.sum.all","ld","wk","prior1","prior2","kappa0","mu0","method","n")
  
  rm(list = setdiff(ls(),objtokeep))
  #end parallel foe each replication
  
  ###
  ################################ Functions ########################################
  mymultt <- function(Xtrain, X.pred){
    myln <- tryCatch(expr = length(Xtrain[,1]), error = function(e){return(1)})
    myls <- Xtrain
    mylmu <- tryCatch(expr = apply(myls, 2, mean), error = function(e){return(myls)})
    mymun <- myln/(myln+kappa0)*mylmu
    myS <- tryCatch(expr = cov(myls), error = function(e){return(0)})*(myln-1)
    
    ## for the covariates
    myd <- tryCatch(expr = length(Xtrain[1,]), error = function(e){return(length(Xtrain))})       ## numer of covariates
    nu0 <- myd+1       ## numer of covariates +1;
    lambda0 <- diag(myd)                ## identity matrix
    
    kappan <- kappa0 + myln
    nun <- nu0 + myln
    lambdn <- tryCatch(expr = lambda0 + myS + kappa0*myln/(kappa0 + myln)*(mylmu - mu0)%*%t(mylmu - mu0), error = function(e){return(diag(1, 2, 2))})
    return2 <- dmvt(x = X.pred, sigma=(kappan + 1)*lambdn/(kappan*(nun - myd + 1)), df = nun - myd + 1, log = FALSE)
  }
  
  # calculate the NPC
  countUT <- function(resultsum, myoutot){
    myctut <- array(0, dim = c(3, 3, 100))
    myctutSum <- NULL
    for(i in 1:length(my.pick)){
      mycurdata <- resultsum[,,i]
      mypre <- NULL
      pretrt1 <- apply(mycurdata[,4:6], 1, which.max)
      pretrt2 <- apply(mycurdata[,7:9], 1, which.max)
      mypreTall <- cbind(pretrt1, pretrt2)
      for(j in 1:length(trtsgn)){
        mypre[j] <- mypreTall[j, trtsgn[j]]
      }
      sts <- table(factor(mypre, levels = 1:3), factor(myoutot+1, levels = 1:3))
      mysdls <- as.numeric(rownames(sts))
      str1 <- matrix(0, nrow = 3, ncol = 3)
      str1[mysdls,] <- sts
      
      myctut[,,i] <- str1*diag(3)
      myctutSum[i] <- sum(str1*diag(3))
    }
    return <- cbind(myctutSum)
  }
  
  
  
  
  
  #INFERENCE
  #with test_ we mean out of sample data 
  ################################ setup Parameters ########################################
  #wk <- c(0,40,100)
  #prior1 <- prior2 <- c(1/3,1/3,1/3)
  #kappa0 <- 1
  #mu0 <- c(0, 0)
  
  ##predictive step 
  
  DATA=ld
  #work for two treatment
  Xd <- rbind(ld$fitdata[[1]]$XDiscPred, ld$fitdata[[2]]$XDiscPred)
  ntrain=nrow(Xd)
  
  
  #train_pred <- simdata$pred[[k]][1:124,]
  W <- rbind(ld$fitdata[[1]]$WProgn, ld$fitdata[[2]]$WProgn)
  #train_prog <- simdata$prog[[k]][1:124,]
  
  Y <- c(ld$fitdata[[1]]$Y,ld$fitdata[[2]]$Y)
  
  #train_ymat <- simdata$ymat[[k]][1:124,]
  Ymat=matrix(0,ncol = max(ld$fitdata[[1]]$Y)+1,
              nrow =length(ld$fitdata[[1]]$Y)+length(ld$fitdata[[2]]$Y) )
  
  Ymat[cbind(1:length(Y), Y+1)]=1
  
  trt = c(rep(0, length(ld$fitdata[[1]]$Y)), rep(1, length(ld$fitdata[[1]]$Y)))
  
  ##predictive dataset 
  units=1:nrow(ld$postpredictivedata$XDiscPred) # for debug
  Xdpred =(ld$postpredictivedata$XDiscPred)[units,]
  npred=nrow(Xdpred)[]
  
  Wpred=ld$postpredictivedata$WProgn[units,,drop=FALSE]
  
  
  nrep=1;myrep=1
  
  
  utpred1APT.all <- array(0, dim = c(npred, 19, nrep))
  
  ### clustering using CONSENSUS MATRIX method ###################################
  #FOR CATEGORICAL DATA
  
  #dataset containing all the Xd and the Xdpredictive
  Dftocluster= as.data.frame(rbind(Xd,Xdpred))
  
  Dftocluster[] = lapply(Dftocluster, factor)
  
  # Compute Gower distance for categorical data
  DistMat <- cluster::daisy(Dftocluster, metric = "gower")
  
  
  
  
  
  
  rst.hc<-ConsensusClusterPlus(DistMat,maxK=15,pFeature=1,
                               clusterAlg=method,#distance="pearson",
                               #clusterAlg="km",distance="euclidean",
                               #clusterAlg="pam",distance="manhattan",
                               seed=126);
  
  
  
  utpred1APT<-matrix(1,nrow= npred,ncol=19)  ### ut1,ut2,trt,cluster
  
  ### pick the median rank with the largest summary measure
  max.clus<-apply(HC.sum.all[,,myrep],1,which.max)+1
  max.clus <- median(max.clus)
  
  #for loop  predictive subj
  for (mysub in 1:npred){
    trtAPT <- trt
    outcom <- Y #??
    select.sub.n <- length(outcom)
    myRapp <- W #prognostic variables sample
    trtcSub <- cbind(outcom, myRapp)
    
    
    mycovXSub <- Wpred[mysub,]                 ## observed covariates predictive unit
    
    trtc0Sub <- subset(trtcSub, trtcSub[,1] == 0)
    trtc0dstSub <- mymultt(trtc0Sub[,-1, drop=FALSE], mycovXSub)
    trtc1Sub <- subset(trtcSub, trtcSub[,1] == 1)
    trtc1dstSub <- mymultt(trtc1Sub[,-1, drop=FALSE], mycovXSub)
    trtc2Sub <- subset(trtcSub, trtcSub[,1] == 2)
    trtc2dstSub <- mymultt(trtc2Sub[,-1, drop=FALSE], mycovXSub)
    trtpfySub <- c(trtc0dstSub, trtc1dstSub, trtc2dstSub)   ##prob with cov
    
    
    # npred=length(trt)#togliere
    myyAPT <- matrix(0, nrow = n, ncol = 3) #matrix outcome
    for(m in 1:n){
      myyAPT[m, outcom[m]+1] = 1
    }
    
    
    mycons1APT <- rst.hc[[max.clus]][["consensusMatrix"]]
    
    totutAPT <- myyAPT*mycons1APT[n + mysub, 1:n]               ### utility
    
    totpreAPT <- totutAPT[]
    trtpreAPT <- trtAPT[]
    mytemptAPT <- cbind(totpreAPT, trtpreAPT)
    
    ## calculate alpha hat
    trt1gAPT <- subset(mytemptAPT[,1:3], mytemptAPT[,4] == 0)
    trt2gAPT <- subset(mytemptAPT[,1:3], mytemptAPT[,4] == 1)
    alphatrt1APT <- apply(trt1gAPT, 2, sum) + prior1
    probpre1APT <- alphatrt1APT/sum(alphatrt1APT)
    
    alphatrt2APT <- apply(trt2gAPT, 2, sum) + prior2
    probpre2APT <- alphatrt2APT/sum(alphatrt2APT)
    
    probpre1 <- probpre1APT*trtpfySub/sum(probpre1APT*trtpfySub)
    probpre2 <- probpre2APT*trtpfySub/sum(probpre2APT*trtpfySub)
    
    ## calculate the utility
    ut1preAPT <- probpre1%*%wk
    ut2preAPT <- probpre2%*%wk
    utpred1APT[mysub, 1] <- ut1preAPT
    utpred1APT[mysub, 2] <- ut2preAPT
    if(ut2preAPT > ut1preAPT){
      utpred1APT[mysub, 3] = 2
    }
    utpred1APT[mysub, 4:6] <- probpre1
    utpred1APT[mysub, 7:9] <- probpre2
    utpred1APT[mysub, 10:12] <- probpre1APT
    utpred1APT[mysub, 13:15] <- probpre2APT
    utpred1APT[mysub, 16:18] <- trtpfySub/sum(trtpfySub)
    utpred1APT[mysub, 19] <- max.clus#[mysub]
  }
  
  utpred1APT.all[,,myrep] <- utpred1APT
  
  
  
  ##end 
  return(utpred1APT.all[,,myrep])
  
  
}#end 

######################################################################################################
######################################################################################################
######################################################################################################
########1##############################################################################################
#new function parallel

###########
Scen=c("1a","1b","2a","2b","2c","3") 
for (S in Scen) {
  
  
  LISTDATA=readRDS(file = paste0(   "~/Sim/Scen/S",S,"Data")) # <- CHANGE PATH HERE 
  for (r in names(LISTDATA)) {
    LISTDATA[[r]][["repname"]]=r
  }
  
  LISTDATA$info
  #predictive dataset (out of sample)
  DP=LISTDATA$predData$Data
  

  
  for (rep in 1:50) {
    
    for(MEH in c("hc","pam")){#hc pam
      Estim=Mafunction(Datalistrep =LISTDATA[[paste0("rep",rep)]] ,Dpred =DP ,
                       method =MEH,CORES = 100)
      
      
      mptosaveread=mptosaveread="~/Sim/Scen" #"~/Documents/UNIFI/PhD/Peoject/Results/Scen/"
      scen=S
      namerep=LISTDATA[[paste0("rep",rep)]]$repname
      method=MEH
      ##save results
      cat(paste0(mptosaveread,"/results",scen,"/",namerep,"_ma-",method ))
      saveRDS(Estim,paste0(mptosaveread,"/results", scen,"/",namerep,"_ma-",method ))
    }
    
  }
  
}



































