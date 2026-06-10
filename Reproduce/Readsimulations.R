##R Script to read and compute the resulst of the simulations 
#author lorenzo moni



rm()
library("MLmetrics")
postTrtinfMamethods=function(PostInf=utpred1APT.all[,,1],   weights=c(0,40,100),Predtrttouse=NULL){
  
  #MAP outcome under the two treatments
  preouttrt1 <- apply(PostInf[,4:6], 1, which.max)
  preouttrt2 <- apply(PostInf[,7:9], 1, which.max)
  
  #predicted utilities for a generic set of weights
  PU1=PostInf[,4 :6]%*%weights
  PU2=PostInf[,7:9]%*%weights
  
  PU=cbind(PU1,PU2)
  #for check
  
  print(all(PU==cbind(PostInf[,1],PostInf[,2])))
  
  #predicted optimal treatment
  Predopttrt=as.vector((PU2>PU1)*1+1)
  
  if(is.null(Predtrttouse)){Predtrttouse=Predopttrt}
  
  print(Predtrttouse)
  #pred outcome
  Outpred=preouttrt1*(Predtrttouse==1)+preouttrt2*(Predtrttouse==2)
  Outpred=Outpred-1
  
  return(list(PredOptTrt=Predopttrt, PredOutcome=Outpred, Utility=PU, yt1=preouttrt1,yt2=preouttrt2))
  
  
}






# Metrics for treatment inference (ONLY IN SIMULATION SETTINGS) only for 2 treatment as def in pedone 2024
MetricsTrtInf=function(Dataset, namesutility=c("utility1","utility2"),outcomename="outcome", 
                       PredOptTrt, PredResponse, PredUtility ){
  n=nrow(Dataset)
  
  cat("Dim pred Data", n, "-Dims estim: t re ut ", length(PredOptTrt), length(PredResponse), nrow(PredUtility))
  
  #relative gain %MTU
  trueu1=Dataset[,namesutility[1]]
  trueu2=Dataset[,namesutility[2]]
  
  diffutility=trueu1-trueu2
  
  trueopttrt=(diffutility<0)+1
  #delta(trueopttrt, predictedtrt) =-1 if trueopttrt != predictedtrt, =1 o/w
  delta=(trueopttrt==PredOptTrt)-(trueopttrt!=PredOptTrt)
  
  
  DELTAMTU=sum(delta*abs(diffutility))
  
  DELTAMTUoptim=sum(abs(diffutility))
  
  #MOT: #  patients misassigned to their optimal treatment
  mot=sum(delta==-1)
  
  #COT correctly assigned treatment
  cot=sum(delta==1)
  
  #F1cOT
  F1cot=MLmetrics::F1_Score(trueopttrt, PredOptTrt)
  
  
  
  #NPC: # patients whose outcome has been correctly predicted
  npc=sum(Dataset[,"outcome"]==PredResponse)
  
  
  #RMSE utilities
  rmse=apply(cbind(trueu1,trueu2)-PredUtility, 2, function(x){sqrt(mean(x^2))})
  rmse1=rmse[1]
  rmse2=rmse[2]
  #RMSE diff utilities
  # drmse=apply(trueu1-trueu2+PredUtility[,2]-PredUtility[,1], 2, function(x){sqrt(mean(x^2))})
  
  
  
  esm=ComputeESM(PredOptTrt, ActualAssTrt=Dataset[,"asstrt"], 
                TrueResp=Dataset[,"outcome"], responders=c("2"))
    
  
  
  return(list(MOT=mot, COT=cot,  MTU=DELTAMTU, NPC=npc, 
              F1COT=F1cot, percMOT=mot/n, percCOT=cot/n, percMTU=DELTAMTU/DELTAMTUoptim ,
              percNPC=npc/n, RMSE1=rmse1,RMSE2=rmse2,ESM=esm))
}





















##############################################################################################################
mainpath="~/Documents/UNIFI/PhD/Peoject/Results/Scen/"
#scennumb=9
#Truedatapred=readRDS(paste0(mainpath,"/Results/Scen1initial/Scen",scennumb,"Data"))#[[repnamepred]]$Data



#tNOTE to use this to compute the metrics specific for tmvbpr one need to set manually the relevant variables for eac treatment (1 if relevan 0 ow)
ReadDatafunct=function(f, 
                       mainpath="~/Documents/UNIFI/PhD/Peoject/Results/Scen/", scennumb=18,
                       tmvbprINFO=FALSE, namePreddata="predData"){
  source(paste0("~/Documents/UNIFI/PhD/Peoject/Code/RFunctionsMultiTreatMVBPR.R"))
  library(tidyverse)
  
  Metricstmvbpr=NULL
  
  repname=strsplit(f,"_")[[1]][1]
  Method=strsplit(f,"_")[[1]][2]
  addinfo=strsplit(f,"_")[[1]][3]
  
  #remove if not debug or convergence check
  addinforep=addinfo; addinfo=NA;
  if(is.na(addinfo)){addinfo=""}else{addinfo=paste0("_",addinfo)}
  
  
  Datalist=readRDS(paste0(mainpath,"/S",scennumb,"Data"))
  
  TrueParameters=Datalist$trueparam
  
  
  repnamepred=namePreddata
  cat("Predictive dataset names:", repnamepred)

  Truedatapred=Datalist[[repnamepred]]$Data
  
  
  
  Truedatatrainlist=Datalist[[repname]]
  
  
   
  if(tmvbprINFO){
  #true clustering structure view 1
  trtindex=Truedatatrainlist$Data[,"asstrt"]
  truecluster_t1=Truedatatrainlist$clusterallbygroup[trtindex==1, TrueParameters$ClusterStrInfo$RelGroupVars ]
  truecluster_t2=Truedatatrainlist$clusterallbygroup[trtindex==2, TrueParameters$ClusterStrInfo$RelGroupVars] 
  
  ##true gamma: same for both treatment ow need to be modified manually
  truegammat1_zeroone=1*(TrueParameters$ClusterStrInfo$Groupvar %in% TrueParameters$ClusterStrInfo$RelGroupVars[[1]])
  truegammat2_zeroone=1*(TrueParameters$ClusterStrInfo$Groupvar %in% TrueParameters$ClusterStrInfo$RelGroupVars[[1]]) #if different need to be modified manually ie [2]
  }
  
  
  cat(paste0(mainpath,"/results",scennumb,"/",f), "\n")
  Estim=readRDS(paste0(mainpath,"/results",scennumb,"/",f))
  
  
  
  
  #MVBRP
  if(Method=="tmvbpr"){
    OPTtrt=EstimOptTreatment(Estim, weights=c(0,40,100))
    Metric=MetricsTrtInf(Truedatapred,
                         PredOptTrt =OPTtrt$EstimetedAssTrt,
                         PredResponse = OPTtrt$PredOutcome,
                         PredUtility=OPTtrt$Utility )
    
    
    if(tmvbprINFO){
      
    #metrics specific for mvbpr
    Zoptim1=as.vector(OptPostPart(Estim$Treatment1$ZZallviews[,1,,drop=FALSE]))
    arit1=mcclust::arandi(truecluster_t1, Zoptim1)
    Zoptim2=as.vector(OptPostPart(Estim$Treatment2$ZZallviews[,1,,drop=FALSE]))
    arit2=mcclust::arandi(truecluster_t2, Zoptim2)
    #
    AccGamma1=MLmetrics::Accuracy(truegammat1_zeroone, 1*(OptPostCategorical(Estim$Treatment1$Gammas)==1))
    F1SGamma1=MLmetrics::F1_Score(truegammat1_zeroone, 1*(OptPostCategorical(Estim$Treatment1$Gammas)==1),positive = "1")
    fdr1=1-MLmetrics::Precision(truegammat1_zeroone, 1*(OptPostCategorical(Estim$Treatment1$Gammas)==1),positive = "1")

    AccGamma2=MLmetrics::Accuracy(truegammat2_zeroone, OptPostCategorical(Estim$Treatment2$Gammas)==1)
    F1SGamma2=MLmetrics::F1_Score(truegammat2_zeroone, 1*(OptPostCategorical(Estim$Treatment2$Gammas)==1),positive = "1")
    fdr2=1-MLmetrics::Precision(truegammat2_zeroone, 1*(OptPostCategorical(Estim$Treatment2$Gammas)==1),positive = "1")
    
    
    if(is.nan(F1SGamma1)){F1SGamma1=0}
    if(is.nan(F1SGamma2)){F1SGamma2=0}
    # print(cbind(relevantVariables[[1]] ,OptPostCategorical(Estim$Treatment1$Gammas)))
    
    #performances of initial gamma vs true
    if(!is.null(Estim$initvalues$trtt1_Gamma)|!is.null(Estim$initvalues$trtt2_Gamma)){
    initaccg1=MLmetrics::Accuracy(truegammat1_zeroone, Estim$initvalues$trtt1_Gamma[,1])
    initaccg2=MLmetrics::Accuracy(truegammat2_zeroone, Estim$initvalues$trtt2_Gamma[,1])
    initfdr1=1-MLmetrics::Precision(truegammat1_zeroone, Estim$initvalues$trtt1_Gamma[,1],positive = "1")
    initfdr2=1-MLmetrics::Precision(truegammat2_zeroone, Estim$initvalues$trtt2_Gamma[,1],positive = "1")
    }else{
      initaccg1=initaccg2=initfdr1=initfdr2=-100
      
    }
     
    Metricstmvbpr=
    tibble(reshape2::melt(tibble(arit1,arit2, AccGamma1, F1SGamma1, fdr1,AccGamma2, F1SGamma2 ,fdr2, 
                                 initaccg1,initaccg2,initfdr1,initfdr2))) %>% rename(Metric=variable) %>%
      mutate(Model=paste0(Method,addinfo),rep=paste0(repname,"_",addinforep))
    }
  }
  #STAN
  if(Method=="lr"||Method=="lrhs"){
    postSTAN=postTrtinfStan(PostY = Estim$ypred)
    Metric=MetricsTrtInf(Truedatapred,
                         PredOptTrt =postSTAN$PredOptTrt,
                         PredResponse = postSTAN$PredOutcome,
                         PredUtility=postSTAN$Utility )
    
  }
  
  
  #ma
  if(Method=="ma-hc"){
    postMA=postTrtinfMamethods(PostInf = Estim)
    Metric=MetricsTrtInf(Truedatapred,
                         PredOptTrt =postMA$PredOptTrt,
                         PredResponse = postMA$PredOutcome,
                         PredUtility=postMA$Utility )
    
  }
  
  if(Method=="ma-pam"){
    postMA=postTrtinfMamethods(PostInf = Estim)
    Metric=MetricsTrtInf(Truedatapred,
                         PredOptTrt =postMA$PredOptTrt,
                         PredResponse = postMA$PredOutcome,
                         PredUtility=postMA$Utility )
    
  }
  Metrics=tibble(reshape2::melt(Metric)) %>% rename(Metric=L1) %>% mutate(Model=paste0(Method,addinfo),rep=paste0(repname,"_",addinforep))
  
  
  #Metrics=Metricstmvbpr
  
  return(list(Metrics,Metricstmvbpr))
}



#####################################################################
#####################################################################
#check all files and method are present 
library(tidyverse)

files=list.files("~/Documents/UNIFI/PhD/Peoject/Results/Scen/results18/")

files=str_split(files,"_",simplify = TRUE)

#all(sort(as.numeric(str_split(files[,1],"rep",simplify = TRUE)[,2])[files[,2]=="lr"])
#==1:50)
setdiff(1:50,as.numeric(str_split(files[,1],"rep",simplify = TRUE)[,2])[files[,2]=="tmvbpr"] )
setdiff(1:50,as.numeric(str_split(files[,1],"rep",simplify = TRUE)[,2])[files[,2]=="lr"] )
setdiff(1:50,as.numeric(str_split(files[,1],"rep",simplify = TRUE)[,2])[files[,2]=="ma-hc"] )
setdiff(1:50,as.numeric(str_split(files[,1],"rep",simplify = TRUE)[,2])[files[,2]=="ma-pam"] )


aa=readRDS("~/Documents/UNIFI/PhD/Peoject/Results/Scen/St")
aa$rep1$Data==aa$rep2$Data
### 
scennumb=SCENnumb
Td=readRDS(paste0(mainpath,"/Results/Scen1initial/Scen",scennumb,"Data"))[["predData"]]$Data
cat("pred opt trt")
table(1*(Td[,"utility1"]>Td[,"utility2"])+2*(Td[,"utility1"]<Td[,"utility2"]))

Drep=readRDS(paste0(mainpath,"/Results/Scen1initial/Scen",scennumb,"Data"))[["rep20"]]$Data
cat("rep opt trt")
table(1*(Drep[,"utility1"]>Drep[,"utility2"])+2*(Drep[,"utility1"]<Drep[,"utility2"]))


ReadDatafunct("rep1_lr",scennumb = 9)
ReadDatafunct("rep2_tmvbpr",scennumb = 9)
table(Truedatapred=readRDS(paste0(mainpath,"/Results/Scen1initial/Scen",scennumb,"Data"))[["rep1"]]$Data[,"asstrt"])
##



# 
mp="~/Documents/UNIFI/PhD/Peoject"
Estim=readRDS(paste0(mp,"/Results/Scen1initial/results",SCENnumb,"/","rep1_ma-hc_onlyrel"))

Td=readRDS(paste0(mp,"/Results/Scen1initial/Scen",SCENnumb,"Data"))[["rep1"]]$Data


postMA=postTrtinfMamethods(PostInf = Estim)
Metric=MetricsTrtInf(Td,
                     PredOptTrt =postMA$PredOptTrt,
                     PredResponse = postMA$PredOutcome,
                     PredUtility=postMA$Utility )


x11()
PlotMPPViewAlloc(Estim$Treatment1$Gammas[,25000:30000])
PlotMPPViewAlloc(Estim$Treatment2$Gammas)
#####################################################################

##
library(tidyverse)
SCENnumb="1a"

lf=list.files(paste0("~/Documents/UNIFI/PhD/Peoject/Results/Scen/results",SCENnumb,"/")) #change here path

length(grep("38",lf) ); length(grep("lr",lf) )

lf=lf[grep("lr",lf) ][20:38]




reslist=parallel::mclapply(lf,ReadDatafunct,mc.cores =15,scennumb=SCENnumb,tmvbprINFO=FALSE,
                           namePreddata="predData")
RES=bind_rows(reslist)

 


#per gamma trt
xtable(apply(RES,2,summary),digits = ,
       caption = paste0("GAMMA scen", SCENnumb),
       include.rownames = FALSE
) %>% print(.,       type = "latex",
            include.rownames = TRUE,   
            booktabs = TRUE)



RES %>% filter(Model=="lr",Metric=="percMTU") %>% 
  filter(value>0.8)
RES %>% filter(Model=="tmvbpr",Metric=="percMTU") %>% 
  filter(value>0.8)

RES %>% filter(rep=="rep45")
##plot and results
library(ggplot2)
RES %>%filter(Metric %in% c("percCOT", "percMTU", "percNPC", "percMOT","RMSE","ESM"))  %>%
  ggplot(.,aes(x=value, color=Model))+
  facet_wrap(Metric~.,scales = "free")+
  geom_boxplot()



#table results
RESlogformat=RES %>%#  filter(rep %in%goodrep)
 group_by(Metric,Model) %>% summarize(mean=mean(value),median=median(value),sd=sd(value))
RESlogformat%>% print(n=1000)

#metrics to include
RESwidformat=RESlogformat  %>%filter(Metric %in% c("percCOT", "percMTU", "percNPC", "percMOT","RMSE")) %>% 
  pivot_longer(cols = c(mean, median, sd), names_to = "Statistic", values_to = "Value") %>%
  pivot_wider(names_from = Metric, values_from = Value) %>%
  arrange(Model, factor(Statistic, levels = c("mean", "median", "sd")))  

#check 
RESwidformat%>% mutate(c=percCOT+percMOT)



###T-MVBPR metrics convergence diagnostics
RES %>%filter(Model=="tmvbpr") %>%  mutate(r=str_split(rep,"_",simplify = TRUE)) %>% 
 # filter(r[,2]!="NA") %>% filter(value>0.75) %>% 
  filter(Metric%in%c("arit1","AccGamma1","AccGamma2","arit2","percMTU"))%>% 
  ggplot(., aes(y=value,x=rep, color=Metric))+ 
  geom_point(size=2.5)



###T-MVBPR metrics convergence diagnostics

RES %>%filter(Model=="tmvbpr") RES %>% filter(Metric%in%c("AccGamma1","AccGamma2",
                                                      "F1SGamma1","F1SGamma2",
                                                        "fdr1","fdr2")) %>% 
  group_by(Metric) %>% summarise(mean=mean(value))


 

#-------
###plot utility
#true data
repnamepred="rep1"
SCENnumb=9
DP=  readRDS(paste0("~/Documents/UNIFI/PhD/Peoject/Results/Scen/S",SCENnumb,"Data"))[[repnamepred]]$Data
dim(DP)

DATAparam=readRDS(paste0("~/Documents/UNIFI/PhD/Peoject/Results/Scen/S",SCENnumb,"Data"))$trueparam
DATAparam$bt1
DATAparam$bt2
DATAparam$bw
DATAparam$ntreat
table(DATAparam$ClusterStrInfo$Groupvar)

DATAparam$MaxcatPredvariables %>% length()

DATAparam$ClusterStrInfo


readRDS(paste0("~/Documents/UNIFI/PhD/Peoject/Results/Scen/S",SCENnumb,"Data"))$trueparam$bw

readRDS(paste0("~/Documents/UNIFI/PhD/Peoject/Results/Scen1initial/Scen",SCENnumb,"Data"))$trueparam


reptoplot= "rep1"#strsplit(lf[4],"_")[[1]][1]
reptoplot
EstimMV=readRDS(paste0("~/Documents/UNIFI/PhD/Peoject/Results/Scen/results",SCENnumb,"/",reptoplot,"_tmvbpr"))
EstimST=readRDS(paste0("~/Documents/UNIFI/PhD/Peoject/Results/Scen/results",SCENnumb,"/",reptoplot,"_lr"))
EstimMAhc=readRDS(paste0("~/Documents/UNIFI/PhD/Peoject/Results/Scen1initial/results",SCENnumb,"/",reptoplot,"_ma-hc"))
EstimMApam=readRDS(paste0("~/Documents/UNIFI/PhD/Peoject/Results/Scen1initial/results",SCENnumb,"/",reptoplot,"_ma-pam"))




posterior <- rstan::extract(EstimST$Fit)

Posttau=posterior$tau # MCMCiter x R_y
dim(posterior$lambda_tilde) #row MCMCiteration, col variables, Arrdim R_y 







list1=list(tMVBPR=EstimOptTreatment(EstimMV)$Utility, 
           lr_stan=postTrtinfStan(PostY = EstimST$ypred)$Utility)#,
           ma_hc=postTrtinfMamethods(PostInf = EstimMAhc)$Utility)

PlotTrueveEstdiffUtilityMethods(DP,PredUtilitylist = list1)




PlotTrueveEstdiffUtility(DP,PredUtility =EstimOptTreatment(EstimMV)$Utility )

PlotTrueveEstdiffUtility(DP,PredUtility =postTrtinfStan(PostY = EstimST$ypred)$Utility )


plot((DP[,"utility1"]-DP[,"utility2"]))
sum(DP[,"utility1"]-DP[,"utility2"]>0)






reptoplot="rep31"

goodrep=c()
for(reptoplot in paste0("rep" ,)){
  
EstimMV=readRDS(paste0("~/Documents/UNIFI/PhD/Peoject/Results/Scen/results17/",reptoplot,"_tmvbpr"))
EstimMV$initvalues$trtt1_Gamma
PlotMPPViewAlloc(EstimMV$Treatment1$Gammas[,1:2000])
PlotMPPViewAlloc(EstimMV$Treatment2$Gammas[,])


a=(((OptPostCategorical(EstimMV$Treatment1$Gammas)==1)-(groupvar==2) )%>% abs() %>% sum())
b=(((OptPostCategorical(EstimMV$Treatment2$Gammas)==1)-(groupvar==2) )%>% abs() %>% sum())

if(a<2 & b<2){cat("---", reptoplot)
print(a)
print(b)

  goodrep=c(goodrep,reptoplot)
  
}
}
#-------
rowMeans(EstimMV$Beta[1,,])
rowMeans(EstimMV$Beta[2,,])
PlotBeta(EstimMV$Beta)

coda::effectiveSize(EstimMV$Beta[1,1,])
coda::effectiveSize(EstimMV$Beta[1,2,])
coda::effectiveSize(EstimMV$Beta[2,1,])
coda::effectiveSize(EstimMV$Beta[2,1,])


PlotBeta(EstimMV$Beta)

EstimMV$Treatment1$Acc_Theta
EstimMV$Treatment2$Acc_Theta


coda::mcmc(EstimMV$Treatment1$Thetaus[1,1,]) %>% plot()









