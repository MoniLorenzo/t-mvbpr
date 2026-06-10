#SCRIPT simulation multi treatment: Sensitivity 
#
#Sys.setenv(OPENBLAS_NUM_THREADS = 1)
#Sys.setenv(OMP_NUM_THREADS = 1)


source("~/Sim/RFunctionsMultiTreatMVBPR.R")
library(Rcpp)
#install.packages("RcppDist")
sourceCpp("~/cppCodetMVBPR")
set.seed(10)


#Sensitivity performed on scenario 2.b
S=15

LISTDATA=readRDS(file = paste0("~/Sim/Scen/S",S,"Data"))


#LISTDATA=readRDS(file = paste0("~/Documents/UNIFI/PhD/Peoject/Results/Scen/S",S,"Data"))


#for (r in names(LISTDATA)) {
#  LISTDATA[[r]][["repname"]]=r
#}
#predictive dataset
DP=LISTDATA$predData$Data
MacCategoriesXd=LISTDATA$trueparam$MaxcatPredvariables


#Grid of hyperparameters and rep name
#Grid=expand.grid(dfbt=c(1,7,100), scalebt=c(1,2.5,10), ad=c("1Rd","1"), nu=c("nu1","nu2","nu3","nu4","nu5"),
#                 Repn=paste0("rep",1:50))


#GrilList=apply(Grid,MARGIN = 1, function(x){
#  list(dfbt=x["dfbt"], scalebt=x["scalebt"], ad=x["ad"], nu=x["nu"], Repn=x["Repn"])} ,simplify = FALSE)

#
#rlab3 1:25
#rlab2 rep31:40  and rep 41:50

GridList=vector("list", 4500);i=0
for (dfbt in c(1,7,100)) { #values degree of fredom beta/theta
    for (scalebt in c(1,2.5,10)) { #values scale beta/theta   
        for (ad in c("1Rd","1")) { # hyperparameter Phi_v,k,d:  a_v,k,d,r= 1/R_d (default) or a_v,k,d,r= 1 
          for (nu in c("nu1","nu2","nu3","nu4","nu5")) {
            for (Repn in paste0("rep",31:45)) {
            i=i+1
            GridList[[i]]=list(dfbt=dfbt, scalebt=scalebt, ad=ad, nu=nu, Repn=Repn)
            
  

            
            }
            
          }
        }
    }
}

###
lf=list.files("~/Sim/Scen/resultsScenSens15/")


Li2 = parallel::mclapply(GridList, function(Parinfo){
  
  namerep  = Parinfo$Repn
  dfbt     = Parinfo$dfbt
  scalebt  = Parinfo$scalebt
  ad       = Parinfo$ad
  nu       = Parinfo$nu
  
  iname = paste0("dfF",dfbt,"-scaleF",scalebt,
                 "-adF",ad,"-nuF",nu )
  
  Name = paste0("~/Sim/Scen/resultsScenSens15","/",namerep,"_",iname)
  
  if (file.exists(Name)) {
    return(NULL)        # <<< better than NA
  } else {
    return(Parinfo)     # keep this job
  }
  })

MissingList <- Filter(Negate(is.null), Li2)


saveRDS(MissingList, "~/Sim/Missinglist")

length(MissingList)



rm(list = c("LISTDATA","dfbt", "scalebt", "ad", "nu","Repn","i")); gc()








toparallSENSITIVITY=function(Parinfo,  Dpred=DP, #ListStructurePrior=ListPrior,
                             maxd=MacCategoriesXd, Sc=S,
                             Infoname=NULL){
          
            #parse the parrallel info
            namerep=Parinfo$Repn
            dfbt=Parinfo$dfbt
            scalebt=Parinfo$scalebt
            ad=Parinfo$ad
            nu=Parinfo$nu
            

            #namerep=Datalistrep$repname
            DATA=readRDS(file = paste0("~/Sim/Scen/S",S,"Data"))[[namerep]]$Data

            #DATA=readRDS(file = paste0("~/Documents/UNIFI/PhD/Peoject/Results/Scen/S",Sc,"Data"))[[namerep]]$Data#Datalistrep$Data
            
            ld=MultiTreatmentMVBPR(DATA, 
                                   DatasetPostPred=Dpred,
                                   nameTreat="asstrt",
                                   nameY="outcome")
            
            
            
            Q=ncol(ld$fit[[1]]$WProgn)
            nv=3
            
  
  
  
            #create the listprior scructure. Note here it contains (mostly) the default hyperparameters 
            if (ad == "1Rd") {
              AD = NULL
            } else {
              AD = 1
            } 
  
  
  
            Listprior=genDefaultPrior(Nvie =nv,singlevalue = AD, #AD change the hyperparameter of the X_D model 
                                      Maxcategories =maxd,
                                      Ry =3 ,
                                      nprogn = 2)
            
            
            #check 
            if(any(is.na(Listprior$a))){stop()}
            
            #change the hyperparameter according to the grid; beta and theta df and scale
            Listprior$HParbetas$df=(Listprior$HParbetas$df*0)+dfbt
            Listprior$HParbetas$scale=(Listprior$HParbetas$scale*0)+scalebt
            
            
            Listprior$HParthetak$df=(Listprior$HParthetak$df*0)+dfbt
            Listprior$HParthetak$scale=(Listprior$HParthetak$scale*0)+scalebt
            
            #change the view allocation (log) prior prob 
            Nuvectors=cbind(nu1=c(0.9,0.05,0.05),
                            nu2=c(0.05,0.95,0.05),
                            nu3=c(0.05,0.05,0.95),
                            nu4=c(1/3,1/3,1/3),
                            nu5=c(0.2,0.4,0.4))
            Nu=Nuvectors[,nu]
            
            Listprior$HParViewallocationlog[] = log(Nu)
            
            
            #sensitivity combination
            iname=paste0("dfF",dfbt,"-scaleF",scalebt,
                          "-adF",ad,"-nuF",nu )
            
            
            ##function to parallel
            



              
              
              ginit=InitializerGammas(ld)*1#matrix(c(rep(0,15),rep(1,5)), ncol=2,nrow = 20)
              ifss=5*10^4  
              bi=2*10^4
              tin=15
              Estim=MVBPRmultitreat(DataList =ld ,
                                    ListPrior = Listprior,
                                    numbtreatments =2 ,
                                    numbprognosticvar =Q,
                                    maximumcatDiscrete = maxd ,
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

              
              cat("~/Sim/Scen/resultsScenSens15","/",namerep,"_",iname )
              saveRDS(Estim,paste0("~/Sim/Scen/resultsScenSens15","/",namerep,"_",iname ))
              rm(Estim);gc()
NULL
            }
            
            
            
            
            



parallel::mclapply(GridList,toparallSENSITIVITY,mc.cores = 100)
            
            
            
            
  
  
  
  
  

  








