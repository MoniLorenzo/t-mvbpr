## Script to reproduce lgg cluster analysis and trt selection in the LGG application  

#t-mvbrp cluster analisys 

source("~/RMultiTreatMVBPR.R")
Seeds=c(14722,24236, 61791, 49038, 33607, 66407, 60562, 80469, 7523, 65871) 

DATA=readRDS("~/Data/DATAlgg1") 

dataset=DATA$data
 
#dataset$C1orf162=relevel(droplevels(factor(dataset$C1orf162)), ref = "0")
dataset$MKLN1=relevel(droplevels(factor(dataset$MKLN1)), ref = "0")
dataset$ADARB2=relevel(droplevels(factor(dataset$ADARB2)), ref = "0")


#dataset[["int"]]=1



progn=c("MKLN1",   "ADARB2" )
a=function(R){
 

set.seed(R)
res=R

paste0("SAVE ~/Data/lgg1_",res)


Estim=tMVBPR(dataset, predictiveDataset=NULL, 
                Outcomemapping=c("Complete Response"=2,"Progressive Disease"=0,"Stable Disease"=1), 
                NumberofViews=4,
                ListPrior=NULL,
                ListInitialValues=list(Numbinitialclusters=10),
                #ListInitialValues=NULL,
                nameTreatment="Treatment", 
                nameOutcome="Outcome",
                namesPredictive=setdiff(DATA$namesCNV,progn),
                namesPrognostic=progn,
                MCMCfinalSampleSize=5*10^4,
                MCMCburnin=10*10^4,
                MCMCThinning=25,
                ListMCMCparam=list(PropVarTheta=0.05, InitVarBeta=2, RelViewNeal8Mparam=15),
                HowgammasInit=hg
                #ListInit Aggiungere in futuro
                )


print("here")
saveRDS(object = Estim,paste0("~/Data/results/tmlgg2g1_",res))
}

parallel::mclapply(c(1:10), a, mc.cores=10)






## code used to perform the CV t-mvbpr and dm-reg for the application
rm(list = ls())
source("~/Data/Code/RFunctionsMultiTreatMVBPR.R")

DATA=readRDS("~/Data/DATAlgg1") 
dataset=DATA$data


dataset$MKLN1=relevel(droplevels(factor(dataset$MKLN1)), ref = "0")
dataset$ADARB2=relevel(droplevels(factor(dataset$ADARB2)), ref = "0")
set.seed(1)
dataset=dataset[sample(1:nrow(dataset)),]
progn=c("MKLN1",   "ADARB2" )

 
#
folds = sort(rep(1:10, length.out = dim(DATA$data)[1]))
ListCV=list()



for(f in (1:10)){
  ListCV[[paste0(f)]]=list(train=dataset[f!=folds,], test=dataset[f==folds,],fold=f)
}




#Estimation tmvbpr
CVparalleltmvbpr=function(Lst, R){
  hg=-2
  
  set.seed(R)
  res=R
  
  Estim=tMVBPR(dataset = Lst$train, predictiveDataset=Lst$test, 
               Outcomemapping=c("Complete Response"=2,"Progressive Disease"=0,"Stable Disease"=1), 
               NumberofViews=4,
               ListPrior=NULL,
               ListInitialValues=list(Numbinitialclusters=20),
               #ListInitialValues=NULL,
               nameTreatment="Treatment", 
               nameOutcome="Outcome",
               namesPredictive=setdiff(DATA$namesCNV,progn),
               namesPrognostic=progn,
               MCMCfinalSampleSize=4*10^4,#,
               MCMCburnin=2*10*10^4,#5*10^4,
               MCMCThinning=25,#25,
               ListMCMCparam=list(PropVarTheta=0.05, InitVarBeta=1, RelViewNeal8Mparam=20),
               HowgammasInit=hg
               #ListInit Aggiungere in futuro
               
               
  )

  print("here")
  saveRDS(object = Estim,file = paste0("~/Data/folds1/Tm4_",Lst$fold,"_aseed",res))
  return(paste0("~/Data/folds/fold_",Lst$fold,"_seed",res))
}

#1:4MCMCfinalSampleSize=4*10^4,#, MCMCburnin=10*10^4,#5*10^4,MCMCThinning=25,#25,

#dm
CVparalleldm=function(Lst){

ld=CreateDatalistDM(dataset = Lst$train, predictiveDataset=Lst$test, 
                    Outcomemapping=c("Complete Response"=2,"Progressive Disease"=0,"Stable Disease"=1), 
                    nameTreatment="Treatment",
                    nameOutcome="Outcome",
                    namesPredictive=setdiff(DATA$namesCNV,progn),
                    namesPrognostic=progn)
stan=StanDR(ld)

saveRDS(object = stan,file = paste0("~/Data/folds/Dm4_",Lst$fold,"_0"))
return(paste0("~/Data/folds/dm_",Lst$fold,"_0"))
}




#
library(parallel)
jobs = list()
for (f in c(1:10)) {
  for (seed in c(14722,61791)){
    jobs[[paste0("tm_",f,"_",seed)]] = mcparallel(CVparalleltmvbpr(ListCV[[f]], R=seed))
  }
 jobs[[paste0("dm_",f,"_0")]] = mcparallel(CVparalleldm(ListCV[[f]]))
  
  
}



res = mccollect(jobs)
#
##stop()
#################################################################################
ResultsESMNPC=parallel::mclapply(ListCV,mc.cores = 10, function(fo){
  library(tidyverse)
   f=fo$fold
   datatest=fo$test
   
   
    ChainsTmvbpr=list(Chain1=readRDS(paste0("~/Data/folds1/Tm4_",f,"chain1")),
                     Chain2=readRDS(paste0("~/Data/folds1/Tm4_",f,"chain2"))
                    )
   
   alpha=DiagnosticBetaAlpha(Post = ChainsTmvbpr, iter=  (200000.00*3):640000.00,
                            pathtosaveplots = NULL) %>%   
     filter(str_detect(Parameter, "^beta_") | str_detect(Parameter, "Treatment.*alpha1"))
   
   Posttmvbpr=combineChains(ChainsTmvbpr,  iter= (200000.00*3):640000.00)
    
   
   OptTrt=EstimOptTreatment(Posttmvbpr,weights = c(0,10,100))
   
    
   #metric dm
   Postdm=readRDS(paste0("~/Data/folds1/Dm4_",f))
   Postdm=postTrtinfStan(Postdm$ypred,weights = c(0,10,100))
   
 
   
   df=data.frame(ActTrt=datatest$Treatment, TrueResp=datatest$Outcome,
        tmvbprPredTrt=OptTrt$EstimetedAssTrt,
        tmvbprPredOutcome=OptTrt$PredOutcome,
        dmPredTrt=Postdm$PredOptTrt,
        dmPredOutcome=Postdm$PredOutcome)
   
   
   return(list(DF=df,conv=alpha))
})


#Diagnostic
DiagnosticCV=reshape2::melt(lapply(ResultsESMNPC, function(x){x$conv})) %>% 
  rename(fold=L1)

DiagnosticCV %>%
  filter(variable == "Effective") %>%
  mutate(value=value) %>% 
  ggplot(.,aes(x = value, y = Parameter, colour = fold)) +
  geom_point() +
  labs(x = "Proportion Effective Sample Size")+
  xlim(0,NA)

DiagnosticCV%>%
  filter(variable == "Rhat") %>%
  ggplot(.,aes(x = value, y = Parameter, colour = fold)) +
  geom_point() +
  labs(x = "Rhat")+
 # xlim(NA,NA)+
  geom_vline(xintercept = 1.15)



##RESULTS 10fold CV ESM and NPC    
Res=bind_rows(lapply(ResultsESMNPC, function(x){x$DF}))
ComputeESM(PredOptTrt=Res$tmvbprPredTrt, 
           #  PredResponse=Postdm$PredOutcome, 
           ActualAssTrt=Res$ActTrt, 
           TrueResp=Res$TrueResp, 
           responders=c("Complete Response"))
ComputeESM(PredOptTrt=Res$dmPredTrt, 
           #  PredResponse=Postdm$PredOutcome, 
           ActualAssTrt=Res$ActTrt, 
           TrueResp=Res$TrueResp, 
           responders=c("Complete Response"))


ComputeNPC(TrueResponse =Res$TrueResp ,
           PredResponse =Res$tmvbprPredOutcome, 
           LevelsResp = c("Complete Response"=2,"Progressive Disease"=0,"Stable Disease"=1) )

ComputeNPC(TrueResponse =Res$TrueResp ,
           PredResponse =Res$dmPredOutcome, 
           LevelsResp = c("Complete Response"=2,"Progressive Disease"=0,"Stable Disease"=1) )





###################################
#sensitivity 
ResultsESMNPCsens=parallel::mclapply(seq(10,90,by=1),mc.cores = 80, function(omega){
  library(tidyverse)
  
  DFs=lapply(ListCV,function(fo){
  f=fo$fold
  datatest=fo$test
  
  
  ChainsTmvbpr=list(Chain1=readRDS(paste0("~/Data/folds1/Tm4_",f,"chain1")),
                    Chain2=readRDS(paste0("~/Data/folds1/Tm4_",f,"chain2")))
  
   Posttmvbpr=combineChains(ChainsTmvbpr,  iter= (200000.00*3):640000.00)
  
  
  OptTrt=EstimOptTreatment(Posttmvbpr,weights = c(0,omega,100))
  
  
  Postdm=readRDS(paste0("~/Data/folds1/Dm4_",f))
  Postdm=postTrtinfStan(Postdm$ypred,weights = c(0,omega,100))
  
  
  
  df=data.frame(ActTrt=datatest$Treatment, TrueResp=datatest$Outcome,
                tmvbprPredTrt=OptTrt$EstimetedAssTrt,
                tmvbprPredOutcome=OptTrt$PredOutcome,
                dmPredTrt=Postdm$PredOptTrt,
                dmPredOutcome=Postdm$PredOutcome)
  
  
  return(list(DF=df,conv=alpha))
})
  
  Res=bind_rows(lapply(DFs, function(x){x$DF}))
  

  esmt=ComputeESM(PredOptTrt=Res$tmvbprPredTrt, 
             #  PredResponse=Postdm$PredOutcome, 
             ActualAssTrt=Res$ActTrt, 
             TrueResp=Res$TrueResp, 
             responders=c("Complete Response"))
  esmd=ComputeESM(PredOptTrt=Res$dmPredTrt, 
             #  PredResponse=Postdm$PredOutcome, 
             ActualAssTrt=Res$ActTrt, 
             TrueResp=Res$TrueResp, 
             responders=c("Complete Response"))
  
  
  npct=ComputeNPC(TrueResponse =Res$TrueResp ,
             PredResponse =Res$tmvbprPredOutcome, 
             LevelsResp = c("Complete Response"=2,"Progressive Disease"=0,"Stable Disease"=1) )
  
  npcd=ComputeNPC(TrueResponse =Res$TrueResp ,
             PredResponse =Res$dmPredOutcome, 
             LevelsResp = c("Complete Response"=2,"Progressive Disease"=0,"Stable Disease"=1) )
  
  return(c(ESMt=esmt,ESMD=esmd, NPCt=npct,NPCt=npcd))
}

)


DFsens=bind_rows(ResultsESMNPCsens) %>% mutate(omega2=seq(10,90,by=1))  %>% reshape2::melt()#mutate(omega2=seq(10,90,by=1)) %>% mel
  ggplot(.,aes(y=ESMt, x=omega2))+geom_line()


  DFsens <- bind_rows(ResultsESMNPCsens) %>% mutate(omega2=seq(10,90,by=1))  %>% 
    mutate(id = row_number()+9) %>%   # x-axis
    pivot_longer(
      cols = c(ESMt, ESMD),
      names_to = "model",
      values_to = "value"
    )
  
  ggplot(DFsens, aes(x = id, y = value, color = model)) +
    geom_line() +
    geom_point() +
    labs(x = "Row", y = "Value", color = "Series") 


  DFsens %>% filter(omega2 %in% c(20,40,60,80))
  
  








