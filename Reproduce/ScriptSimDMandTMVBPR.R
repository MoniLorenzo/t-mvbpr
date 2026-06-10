#SCRIPT to run sim tmvbpr and dm-reg in parallel 

##@author: Lorenzo Moni, Silvia Liverani, Alberto Cassese, Francesco Claudio Stingo  



source("~/CodeStanDM-reg.r") #load stan script fo dm-reg

source("~/Sim/RFunctionsMultiTreatMVBPR.R") #script helper functions 
library(Rcpp)
sourceCpp("~/cppCodetMVBPR")


#require(tidyverse)

set.seed(10)





###
Scen=c("1a","1b","2a","2b","2c","3") 

for (S in Scen) {
  LISTDATA=readRDS(file = paste0("~/Sim/Scen/S",S,"Data")) #<- change here path


    for (r in names(LISTDATA)) {
    LISTDATA[[r]][["repname"]]=r
  }

  DP=LISTDATA$predData$Data
  

  ##
   

  toparall=function(Datalistrep, Dpred=DP, maxd=LISTDATA$trueparam$MaxcatPredvariables, 
                    method=c(2),
                    scen=S){#1 tmvbpr #2 DM  #3 lr

    namerep=Datalistrep$repname

    DATA=Datalistrep$Data

    ld=MultiTreatmentMVBPR(DATA, DatasetPostPred=Dpred,
                           nameTreat="asstrt",
                           nameY="outcome")

    #MAxxdp=apply(ld$postpredictivedata$XDiscPred,2,max)+1
    #MAxxd1=apply(ld$fitdata[[1]]$XDiscPred,2,max)+1
    #MAxxd2=apply(ld$fitdata[[2]]$XDiscPred,2,max)+1

    MAxxd=maxd#apply(cbind(MAxxdp,MAxxd1,MAxxd2),1,max)
    #MT MVBPR
    Q=ncol(ld$fit[[1]]$WProgn)

    nv=3
    
    ginit=InitializerGammas(ld, p=0.1)*1#matrix(c(rep(0,15),rep(1,5)), ncol=2,nrow = 20)
    if(1 %in% method ){
      ifss=5*10^4 #5*10^4
      bi=2*10^4
      tin=15
      Estim=MVBPRmultitreat(DataList =ld ,
                            ListPrior = genDefaultPrior(Nvie =nv ,Maxcategories =MAxxd ,Ry =3 ,
                                                        nprogn = Q),
                            numbtreatments =2 ,
                            numbprognosticvar =Q,
                            maximumcatDiscrete = MAxxd ,
                            NumberofView = nv,
                            TypeofXmodel =1,
                            Maxcatresponse =3 ,
                            M = 10,
                            updatinggammas =1 ,
                            lambdtheta = 0.1,
                            lambdbeta = 5.0,
                            Ninitclusters = 10,
                            InitGammas =-2,
                            InitGammasmat = ginit,
                            AlphaDPsinit =rep(4,nv-1) ,
                            BetaMatinit = matrix((0), nrow =Q, ncol = 2),
                            MCMCFinalSamplesize =ifss,
                            Burnin = bi,
                            Thinning = tin )
      cat("--Methods--", method, "\n" ,"--Scen--",Scen)
      
      
cat("~/Sim/Scen/results",scen,"/",namerep,"_tmvbpr" )
saveRDS(Estim,paste0("~/Sim/Scen/results",scen,"/",namerep,"_tmvbpr" ))
      rm(Estim);gc()
    }
    if(2 %in% method ){
      stan=StanDR(ld)

      saveRDS(stan,paste0("~/Sim/Scen/results",scen,"/",namerep,"_lr" ))
      rm(stan);gc()}


  
  
   
  
}
  
  



  parallel::mclapply(LISTDATA[1:50],
                     toparall,mc.cores = 50)


}
 












