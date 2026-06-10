#code 1/2 to build dataset from LGG data
#author lorenzo moni

rm(list = ls())
library(tidyverse)

#Read clinical and follow up data 
myclin=read_tsv("/home/lorenzo/Documents/UNIFI/PhD/Peoject/tcga/TCGALGG/clinical.tsv")
myfoup=read_tsv("/home/lorenzo/Documents/UNIFI/PhD/Peoject/tcga/TCGALGG/follow_up.tsv")


#remove the variables (cols) that bring no info ( all rows are NA) 
myclinNONull=myclin %>% select(where(~ !all(.=="'--")))
myfoupNONull=myfoup %>% select(where(~ !all(.=="'--")))

#NOTE the structure of myclinNONull is 1PATIENT -> n records,
#     each record is corresponds to a single treatment given/performed on patient 
cat("N records ", nrow(myclinNONull), "  N patients", length(unique(myclinNONull$cases.submitter_id)))


#To define TWO competing treatment we consider the 1st treatment.
#We retain only the treatment's records that have:
#TIME Start and TIME End variables are NOT null
#Non-negative TIME Start/ENd
#the treatment target the primary tumor, ie
##if the trt is a radio we consider to be part of the 1st composite trt
##only if the treatments.treatment_anatomic_sites is not:  Locoregional/Regional Site 
myclinNONull1=myclinNONull %>% filter(treatments.days_to_treatment_start!="'--",
                                      treatments.days_to_treatment_end!="'--") %>% 
  mutate(TIMEtStart=as.numeric(treatments.days_to_treatment_start),
         TIMEtEnd=as.numeric(treatments.days_to_treatment_end)) %>% 
  filter(TIMEtStart>=0, TIMEtEnd>=0) %>% 
  filter(!(treatments.treatment_anatomic_sites%in% 
             c("Locoregional Site","Regional Site")))


cat("N records ", nrow(myclinNONull1), "  N patients", length(unique(myclinNONull1$cases.submitter_id)))


myclinNONull1 %>% select(cases.submitter_id, 
                         treatments.treatment_type, 
                         treatments.treatment_intent_type,
                         TIMEtStart,TIMEtEnd, treatments.course_number,
                         treatments.treatment_outcome)

#Note that a treatment may be composite 
myclinNONull2=myclinNONull1 %>%group_by(cases.submitter_id) %>% 
  #create variable minimum start trt time  for each patient
  mutate(MinStart= min(TIMEtStart)  ) %>%
  #  create variable maximum end trt time, considering the treatments started at the beginning (ie TIMEtStart==MinStart)
  mutate(MaxEndStartMin=max(TIMEtEnd[TIMEtStart==MinStart])) %>% 
  # We def. an indicator variables denoting if the record is relative to 
  #the 1st treatment (since a trt can be a composit).   If
  #the start time <= MaxEndStartMin AND End time<= MaxEndStartMin the record is relative to the 
  #1st treatment  
  mutate(FirstTrt= (TIMEtStart <= MaxEndStartMin)*(TIMEtEnd<= MaxEndStartMin)  ) %>% ungroup()

# should be == to that of myclinNONull1
cat("N records ", nrow(myclinNONull2), "  N patients", length(unique(myclinNONull2$cases.submitter_id)))



myclinNONull2%>%    select(cases.submitter_id, 
                           treatments.treatment_type,MinStart,MaxEndStartMin, FirstTrt,
                           treatments.treatment_intent_type,
                           TIMEtStart,TIMEtEnd, treatments.course_number,
                           treatments.treatment_outcome)


#Check1: if every patient has at least one row (treatments) flagged as 1st trt
Check1=myclinNONull2 %>% group_by(cases.submitter_id) %>% summarise(l1=any(FirstTrt==1)) %>% 
  filter(l1==FALSE)
print(Check1)# 0x2  matrix-> OK



##Retain only records relative to 1st trt, ie FirstTrt==1, and re
myclinNONull3=myclinNONull2 %>% filter(FirstTrt==1)

cat("N records ", nrow(myclinNONull3), "  N patients", length(unique(myclinNONull3$cases.submitter_id)))
myclinNONull3 %>% count(treatments.treatment_outcome)

#relable  treatments.treatment_outcome==Unknown as missing "'--"
myclinNONull4=myclinNONull3 %>%
  mutate(treatments.treatment_outcome=ifelse(treatments.treatment_outcome=="Unknown",yes="'--",no=treatments.treatment_outcome))
myclinNONull4 %>% count(treatments.treatment_outcome)



#Def. my temp trt outcome variable
#
myclinNONull5=myclinNONull4 %>% 
  # we define the trt response as the latest  treatments.treatment_outcome, for each parients
  group_by(cases.submitter_id) %>% 
  mutate(TempTrtOUTCOME = if (n_distinct(treatments.treatment_outcome) == 1) {
    # if all records have the same outcome  take the 1st one 
    first(treatments.treatment_outcome)
  } else {
    
    #latest trt 
    treatments.treatment_outcome[which.max(TIMEtEnd)]
  }) %>% ungroup()

cat("N records ", nrow(myclinNONull5), "  N patients", length(unique(myclinNONull5$cases.submitter_id)))


#count TempTrtOUTCOME 
myclinNONull5 %>%  group_by(cases.submitter_id) %>% 
  summarise(patientoutcome=first(TempTrtOUTCOME)) %>% count(patientoutcome)


#check
#myclinNONull5 %>%  group_by(cases.submitter_id) %>% summarise(max(TIMEtEnd)==
#                                                             mean(MaxEndStartMin))


myclinNONull5%>% select(cases.submitter_id, TIMEtEnd, treatments.treatment_outcome) %>% print(n=1000)


### since we have a lot of missing (139) we try to recover the outcome information  
## using the followUP data 

#save: patient id and the when the 1st trt if ended (ie MaxEndStartMin ) if
#treatment outcome is missing  

#trt info 
InfoPatientsmissingoutcome=myclinNONull5 %>%  #%>%filter(TempTrtOUTCOME=="'--") %>% 
  select(cases.submitter_id, MaxEndStartMin) %>% 
  group_by(cases.submitter_id) %>% 
  summarise(cases.submitter_id=first(cases.submitter_id), TimeTrtEnded=first(MaxEndStartMin))

dim(InfoPatientsmissingoutcome) #check OK 139



##Now we extract information using FUdata
dim(myfoupNONull)
myfoupNONull1=myfoupNONull %>% filter(cases.submitter_id %in% InfoPatientsmissingoutcome$cases.submitter_id)

dim(myfoupNONull1)


#cpunt  for how many patients (with missing outcome),   we have at least one followup record
length(InfoPatientsmissingoutcome$cases.submitter_id)-length(setdiff(InfoPatientsmissingoutcome$cases.submitter_id,myfoupNONull1$cases.submitter_id ))
#ok



#oJint to the followup information 
#the TIME trt end variable (ie, TimeTrtEnded) [parsed by InfoPatientsmissingoutcome] 
###select only relevant variables (cols), and
myfoupNONull2=myfoupNONull1%>% left_join(InfoPatientsmissingoutcome) 


myfoupNONull2 %>%   count(other_clinical_attributes.timepoint_category)
#check if the fu records with other_clinical_attributes!= NA are prior to the trt
myfoupNONull2 %>%   filter(other_clinical_attributes.timepoint_category!="'--") %>% 
  count(follow_ups.days_to_follow_up) #-> yes: so these info are not relevant to the trt outcome 

#so we can safely remove the entries with other_clinical_attributes.timepoint_category NOT null, and
myfoupNONull3=myfoupNONull2 %>% filter(other_clinical_attributes.timepoint_category=="'--")  

cat("N records ", nrow(myfoupNONull3), "  N patients", length(unique(myfoupNONull3$cases.submitter_id)))



#keep only entries records with follow_ups.days_to_follow_up>=TimeTrtEnded or NA
myfoupNONull4=myfoupNONull3 %>% 
  mutate(TimeFU=parse_double(myfoupNONull3$follow_ups.days_to_follow_up,na = "'--")) %>% 
  filter(TimeFU>=TimeTrtEnded| is.na(TimeFU) )


#check:  what info we have   
myfoupNONull4%>% count(follow_ups.days_to_follow_up,follow_ups.days_to_recurrence,follow_ups.days_to_progression,
                       follow_ups.progression_or_recurrence,follow_ups.progression_or_recurrence_type,follow_ups.timepoint_category,
                       follow_ups.disease_response) %>% print(n=1000)


myfoupNONull4 %>% filter(is.na(TimeFU)) %>% 
  count(follow_ups.days_to_follow_up,follow_ups.days_to_recurrence,follow_ups.days_to_progression,
        follow_ups.progression_or_recurrence,follow_ups.progression_or_recurrence_type,follow_ups.timepoint_category,
        follow_ups.disease_response) %>% count(follow_ups.days_to_recurrence,follow_ups.days_to_progression,
                                               follow_ups.progression_or_recurrence,follow_ups.progression_or_recurrence_type)

cat("N records ", nrow(myfoupNONull4), "  N patients", length(unique(myfoupNONull4$cases.submitter_id)))


##remove the records with follow_ups.days_to_follow_up==`'--`, since  the other Time variables
#have no info to  allocate the record in the time continuum  
myfoupNONull5=myfoupNONull4 %>% filter(!is.na(TimeFU))
cat("N records ", nrow(myfoupNONull5), "  N patients", length(unique(myfoupNONull5$cases.submitter_id)))


#check:
myfoupNONull5 %>%  select( cases.submitter_id,
                           TimeFU,TimeTrtEnded, 
                           follow_ups.timepoint_category, 
                           follow_ups.days_to_progression, 
                           follow_ups.days_to_recurrence,
                           follow_ups.progression_or_recurrence,
                           follow_ups.progression_or_recurrence_type,
                           follow_ups.disease_response,
                           project.project_id,
                           other_clinical_attributes.timepoint_category) %>% 
  filter(follow_ups.progression_or_recurrence!="'--") %>% print(n=1000)

myfoupNONull5 %>%  filter(follow_ups.days_to_recurrence!="'--") %>% 
  mutate(l1=follow_ups.days_to_recurrence==follow_ups.days_to_follow_up) %>% 
  count(l1)

myfoupNONull5 %>%  filter(follow_ups.days_to_progression!="'--") %>% 
  mutate(l1=follow_ups.days_to_progression==follow_ups.days_to_follow_up) %>% 
  count(l1)


#the time of follow_ups.days_to_progression OR follow_ups.days_to_recurrence, IF not NA,
#if equal to follow_ups.days_to_follow_up... so we consider the  indicator follow_ups.progression_or_recurrence
##also we add a variable indicate sequential followup #   
myfoupNONull6= myfoupNONull5 %>%  select( cases.submitter_id,
                                          TimeFU,TimeTrtEnded, 
                                          follow_ups.timepoint_category, 
                                          follow_ups.progression_or_recurrence,
                                          #    follow_ups.progression_or_recurrence_type,
                                          follow_ups.disease_response) %>% 
  group_by(cases.submitter_id) %>% 
  mutate(FUorder=rank(TimeFU,ties.method = "min")) #rank the treatments regimes by time of subministration
cat("N records ", nrow(myfoupNONull6), "  N patients", length(unique(myfoupNONull6$cases.submitter_id)))



###
#Defining our trt response as follows :


#SIDE NOTE: 
myfoupNONull6 %>% filter(FUorder==1) %>% group_by(cases.submitter_id) %>% 
  summarise(nd=n_distinct(follow_ups.timepoint_category)) %>% filter(nd>1)
#for some patients we may have multiple regord for the same, 1st, followup

## PD
myfoupNONull6 %>% filter(FUorder==1) %>% filter(follow_ups.progression_or_recurrence=="Yes",
                                                follow_ups.disease_response!="'--")
# 0x9 matrix-> there are no recorde where the 
#follow_ups.progression_or_recurrence==YES and the follow_ups.disease_response in not NA
#

myfoupNONull6 %>% filter(FUorder==1) %>% group_by(cases.submitter_id) %>% 
  filter(any(follow_ups.progression_or_recurrence=="Yes")) %>% 
  summarise(nd1=n_distinct(follow_ups.disease_response)) %>% filter(nd1>1) 
# patients with at least one  follow_ups.progression_or_recurrence=="Yes"   records 
#should not have other records (at same time) where follow_ups.disease_response is Tumor Free:
myfoupNONull6 %>% filter(FUorder==1) %>% filter(cases.submitter_id=="TCGA-S9-A6TX")  
myfoupNONull6 %>% filter(FUorder==1) %>% filter(cases.submitter_id=="TCGA-VM-A8CE") 
#checked: these repeated record are coherent  wrt progression_or_recurrence and disease_response


myfoupNONull6 %>% filter(FUorder==1) %>% group_by(cases.submitter_id) %>% 
  filter(any(follow_ups.progression_or_recurrence=="Yes")) %>% 
  summarise(nd1=n_distinct(follow_ups.progression_or_recurrence)) %>% filter(nd1>1) 
# patients with at least one  follow_ups.progression_or_recurrence=="Yes"   records 
#typically have only one record, if multiple records are present the other records 
#should  have missing follow_ups.progression_or_recurrence , for the coherence  
#should NOT  have other records where follow_ups.disease_response is Tumor Free:
myfoupNONull6 %>% filter(FUorder==1) %>% filter(cases.submitter_id=="TCGA-DU-5854") 
myfoupNONull6 %>% filter(FUorder==1) %>% filter(cases.submitter_id=="TCGA-DU-7300") 
myfoupNONull6 %>% filter(FUorder==1) %>% filter(cases.submitter_id=="TCGA-DU-7301") 
myfoupNONull6 %>% filter(FUorder==1) %>% filter(cases.submitter_id=="TCGA-S9-A6TX") #same as before expected
myfoupNONull6 %>% filter(FUorder==1) %>% filter(cases.submitter_id=="TCGA-VM-A8CE")  #same as before expected
#checked: these repeated record are coherent  wrt progression_or_recurrence 
#interpretation: IF ANY[records] progression_or_recurrence (YES)  --> MYtreatmentresp= Progressive 



#CR 
myfoupNONull6 %>% filter(FUorder==1) %>% group_by(cases.submitter_id)  %>% 
  filter(any(follow_ups.disease_response=="TF-Tumor Free" )) %>% 
  summarise(nd1=n_distinct(follow_ups.disease_response)) %>% filter(nd1>1)   
# patients with at least one follow_ups.disease_response=="TF-Tumor Free"  records,
#typically have only one record, if multiple records are present the other records
#should have missing follow_ups.disease_response, for the coherence  
#should NOT  have other records where follow_ups.disease_response is With Tumor:
myfoupNONull6 %>%  filter(FUorder==1) %>% filter(cases.submitter_id=="TCGA-DB-5275") #ok
myfoupNONull6 %>%  filter(FUorder==1) %>% filter(cases.submitter_id=="TCGA-FG-A87Q") #ok
myfoupNONull6 %>%  filter(FUorder==1) %>% filter(cases.submitter_id=="TCGA-HT-7473") #ok
myfoupNONull6 %>%  filter(FUorder==1) %>% filter(cases.submitter_id=="TCGA-HT-7481") # non-coherent
myfoupNONull6 %>%  filter(FUorder==1) %>% filter(cases.submitter_id=="TCGA-HT-7603") #ok
myfoupNONull6 %>%  filter(FUorder==1) %>% filter(cases.submitter_id=="TCGA-HT-7606") #ok
myfoupNONull6 %>%  filter(FUorder==1) %>% filter(cases.submitter_id=="TCGA-HT-7692") #ok
myfoupNONull6 %>%  filter(FUorder==1) %>% filter(cases.submitter_id=="TCGA-HT-7695") #ok
myfoupNONull6 %>%  filter(FUorder==1) %>% filter(cases.submitter_id=="TCGA-QH-A65X") #ok
myfoupNONull6 %>%  filter(FUorder==1) %>% filter(cases.submitter_id=="TCGA-QH-A6X4") #ok
myfoupNONull6 %>%  filter(FUorder==1) %>% filter(cases.submitter_id=="TCGA-QH-A6XC") #ok
# the TCGA-HT-7481 is not coherent have 2 records where the response is simultaneously 
# TF and WT, this patient will be removed 
#interpretation: IF ANY  disease_response (TF-Tumor Free) AND  ALL disease_response (NOT WT -with tumor)
# -> MYtreatmentresp= Complete [if disease_response==WT then MYtreatmentresp=MISSING ]  


#S or PD: 
myfoupNONull6 %>% filter(FUorder==1) %>% group_by(cases.submitter_id)  %>% 
  filter(  any(follow_ups.disease_response=="WT-With Tumor" )) %>% 
  summarise(nd1=n_distinct(follow_ups.disease_response)) %>% filter(nd1>1)  
# patients with at least one follow_ups.disease_response=="WT-With Tumor"  records,
#typically have only one record, if multiple records are present the other records
#should have missing follow_ups.disease_response, for the coherence  
#should NOT  have other records where follow_ups.disease_response is Tumor free:
myfoupNONull6 %>%  filter(FUorder==1) %>% filter(cases.submitter_id=="TCGA-DH-A66B")  #ok
myfoupNONull6 %>%  filter(FUorder==1) %>% filter(cases.submitter_id=="TCGA-FG-6689")  #ok
myfoupNONull6 %>%  filter(FUorder==1) %>% filter(cases.submitter_id=="TCGA-HT-7481")  # non-coherent! [expected, as before]
myfoupNONull6 %>%  filter(FUorder==1) %>% filter(cases.submitter_id=="TCGA-S9-A6TX")  #ok 
myfoupNONull6 %>%  filter(FUorder==1) %>% filter(cases.submitter_id=="TCGA-VM-A8CE")  #ok
#interpretation: IF ANY  disease_response (WT-With Tumor) AND  ALL disease_response (NOT TF-Tumor Free),
#THEN we consider 2 cases:
#  (a) if   ANY  follow_ups.timepoint_category==Last Contact --> MYtreatmentresp= Progressive
#[if a patient still have a tumor after treatment, and the parient Lost to follow-up just  after the treatment
# we assume that the disease has progressed ]
#  (b)   ALL  follow_ups.timepoint_category!=Last Contact --> MYtreatmentresp= STABLE
#[if a patient still have a tumor after treatment, but the patient has subsequent  followups  
# we assume that the disease is stable]


#We extrapolate the treatment outcome 
myfoupNONull7=myfoupNONull6 %>% filter(FUorder==1) %>% group_by(cases.submitter_id) %>% 
  mutate(PD=any(follow_ups.progression_or_recurrence=="Yes")*1,
         CR=1*(any(follow_ups.disease_response=="TF-Tumor Free") && all(follow_ups.disease_response!="WT-With Tumor")),
         SPDtemp=any(follow_ups.disease_response=="WT-With Tumor") && all(follow_ups.disease_response!="TF-Tumor Free"),
         SPDlastconct=any(follow_ups.timepoint_category=="Last Contact"),
         SPD_PD= (SPDtemp==1)*(SPDlastconct==1),
         SPD_S= (SPDtemp==1)*(SPDlastconct==0)) %>% 
  summarise(PD=first(PD),    CR=first(CR), first(SPDtemp), first(SPDlastconct) ,
            SPD_PD=first(SPD_PD) ,SPD_S=first(SPD_S))

#check
myfoupNONull7 %>% count(PD,CR, SPD_PD,SPD_S)



#dataset containing the recovered trt from fudata
ResponseFromFU=myfoupNONull7 %>% mutate(PD1=PD+SPD_PD) %>% select(cases.submitter_id,PD1,CR, SPD_S, PD, SPD_PD) %>% 
  mutate(FUTrtOutcome=if_else(PD1==1,"Progressive Disease", 
                              if_else(SPD_S==1, "Stable Disease",
                                      if_else(CR==1, "Complete Response", "'--" ))))


##how may response we have ``recovered'' from followup
ResponseFromFU %>% count(FUTrtOutcome)
#104 recovered  







##COMBINE clinical and FollowUP data 

myclinNONull6=myclinNONull5 %>%#  select(cases.submitter_id, treatments.treatment_outcome) %>% 
  left_join(ResponseFromFU %>% select(cases.submitter_id, FUTrtOutcome))

#Potential unit using the outcome in clinical data
myclinNONull6 %>% group_by(cases.submitter_id) %>% summarise(resp=first(TempTrtOUTCOME)) %>% 
  ungroup() %>% count(resp) #35+11+25+110 =181
#Potential unit using the outcome derived from followup  data
myclinNONull6 %>% group_by(cases.submitter_id) %>% summarise(resp=first(FUTrtOutcome)) %>% 
  ungroup() %>% count(resp) #35+11+25+110=263


#check the difference between my trtoutcome is somehow coherent with that recorded in the clinical data for 
myclinNONull6 %>% group_by(cases.submitter_id) %>% 
  summarise(resporiginal=first(TempTrtOUTCOME), respFU=first(FUTrtOutcome)) %>% 
  filter(resporiginal!="'--", respFU!="'--") %>% select(!cases.submitter_id) %>% 
  count(resporiginal,respFU)  

myclinNONull6 %>%
  group_by(cases.submitter_id) %>%
  summarise(resporiginal = first(TempTrtOUTCOME),
            respFU = first(FUTrtOutcome)) %>%
  filter(resporiginal != "'--", respFU != "'--") %>%
  select(!cases.submitter_id) %>%
  count(resporiginal, respFU) %>%
  pivot_wider(names_from = respFU, values_from = n, values_fill = 0)



#check
myclinNONull6 %>% count(TrtOutcome1)
myclinNONull6 %>% select(cases.submitter_id, treatments.treatment_outcome,FUTrtOutcome) %>%  filter(cases.submitter_id=="TCGA-HT-8114") 
myclinNONull6 %>%  select(cases.submitter_id, treatments.treatment_outcome) %>% filter(cases.submitter_id=="TCGA-DU-8158") 
myclinNONull6 %>% select(cases.submitter_id, treatments.treatment_outcome) %>%  filter(cases.submitter_id=="TCGA-DU-6542") 
myclinNONull6 %>%  select(cases.submitter_id, treatments.treatment_outcome) %>%  filter(cases.submitter_id=="TCGA-DU-6542") 



###remove the patients with no response
#check 
myclinNONull7= myclinNONull6 %>% filter(FUTrtOutcome!="'--")



##Count 
myclinNONull7%>% group_by(cases.submitter_id) %>%
  summarise(uniqtrt=first(FUTrtOutcome)) %>% count(uniqtrt)




##Def competing trt
myclinNONull8 %>% count(treatments.treatment_type, treatments.treatment_intent_type)
#263 patients



###Define treatment profiles 
myclinNONull7 %>% count(treatments.treatment_type, treatments.treatment_intent_type)


myclinNONull8=myclinNONull7 %>%# select(cases.submitter_id,TIMEtStart,TIMEtEnd,treatments.treatment_type,treatments.treatment_intent_type) %>% 
  group_by(cases.submitter_id) %>% #mutate(any(grepl("Adjuvant", treatments.treatment_intent_type)))
  mutate(OnlyChemo=all(grepl("Chemotherapy", treatments.treatment_type)),
         OnlyRadio=all(grepl("Radiation", treatments.treatment_type)),
         HasRadioAdj=any(grepl("Radiation", treatments.treatment_type)*grepl("Adjuvant", treatments.treatment_intent_type)),
         HasChemo=any(grepl("Chemotherapy", treatments.treatment_type)),
         HasRadio=any(grepl("Radiation", treatments.treatment_type))
  )


cat("N records ", nrow(myclinNONull8), "  N patients", length(unique(myclinNONull8$cases.submitter_id)))
#check possible trt
myclinNONull8 %>% ungroup() %>% count(treatments.treatment_type)

##remove patients with Immunotherapy (Including Vaccines)  , Pharmaceutical Therapy, NOS   Brachytherapy, NOS
casetorm1=myclinNONull8 %>% filter(treatments.treatment_type %in% c(
  "Immunotherapy (Including Vaccines)", "Pharmaceutical Therapy, NOS","Brachytherapy, NOS"
)) %>% pull(cases.submitter_id)


myclinNONull9=myclinNONull8 %>%filter(!cases.submitter_id %in% casetorm1) 
cat("N records ", nrow(myclinNONull9), "  N patients", length(unique(myclinNONull9$cases.submitter_id)))


myclinNONull9%>%  summarise(OnlyChemo=dplyr::first(OnlyChemo),
                            OnlyRadio=dplyr::first(OnlyRadio),
                            HasRadioAdj=dplyr::first(HasRadioAdj),
                            HasChemo=dplyr::first(HasChemo)) %>% 
  dplyr::count( HasRadioAdj)
#->Possible treatments (A):  HasRadioAdj vs NOT HasRadioAdj # 263 patients (74 VS 184)
#THIS   trt 1 std treatment vs 

myclinNONull9%>%  summarise(TRToutcome=first(FUTrtOutcome),
                            HasRadioAdj=first(HasRadioAdj)) %>% 
  count(HasRadioAdj, TRToutcome) %>%
  group_by(HasRadioAdj) %>%
  mutate(prop = n / sum(n))


cat("N records ", nrow(myclinNONull9), "  N patients", length(unique(myclinNONull9$cases.submitter_id)))





##Save patient id to match with omics data
IdPatientstoinclude1=myclinNONull9 %>% pull(cases.submitter_id) %>% unique()
length(  IdPatientstoinclude1)
#here 258 patients other that do not have record in omics data need to be removed




######################################################################
library(TCGAbiolinks)
##
#query <- GDCquery(
#  project = "TCGA-LGG",
#  data.category = "Copy Number Variation",
#  data.type = "Gene Level Copy Number",
#  sample.type = "Primary Tumor",
##  workflow.type = "GISTIC2_CopyNumber_Gistic2_all_thresholded.by_genes"
#)
#GDCdownload(query = query,files.per.chunk = 10)
#
#library(SummarizedExperiment)
#
#data=GDCprepare(query)
#cnv_mat <- assay(data, "copy_number")
#
#cnv_matmin <- assay(data, "min_copy_number")
#cnv_matmax <- assay(data, "max_copy_number")
#
#saveRDS(list(CNV_mat=cnv_mat, Min=cnv_matmin,Max=cnv_matmax), 
#        "~/Documents/UNIFI/PhD/Peoject/tcga/myDatalggCNVCopyGeneLevel")
#
#
#CNV_mat=readRDS(
#  "~/Documents/UNIFI/PhD/Peoject/tcga/myDatalggCNVCopyGeneLevel")
#CNV_mat=CNV_mat$CNV_mat
#
#dim(CNV_mat)
#
#colnamesCNVmat=colnames(CNV_mat)
#length(colnamesCNVmat)==length(unique(colnamesCNVmat))
#

###----
library(RTCGAToolbox)
datasets <- getFirehoseDatasets()
datasets  # check that "LGG" is in the list

run_dates <- getFirehoseRunningDates(last = 3)
analyze_dates <- getFirehoseAnalyzeDates(last = 3)

run_dates
analyze_dates  # these are the GISTIC2 “analyze” run dates

# Download Firehose data for LGG, including GISTIC
firehose_obj <- getFirehoseData(
  dataset = "LGG",
  runDate = run_dates[1],
  gistic2Date = analyze_dates[1],
  GISTIC = TRUE
)

# Extract GISTIC‑by-gene data (thresholded calls)
gistic_thresh <- getData(firehose_obj, type = "GISTIC", platform = "ThresholdedByGene")
sum(gistic_thresh[,-c(1,2,3
                      )] %>% as.numeric(),na.rm = TRUE)
# Or extract all-by-gene GISTIC scores
#gistic_all <- getData(firehose_obj, type = "GISTIC", platform = "AllByGene")


dim(gistic_thresh)
GenesNames=gistic_thresh$Gene.Symbol
PatientIDCNV=colnames(gistic_thresh)

#create patient id to match that of Clinicaldata
NewpatId=c(PatientIDCNV[c(1:3)],
           apply(str_split(PatientIDCNV[-c(1:3)],"\\.",simplify = TRUE),1,
                 function(x) str_c(x[1:3],collapse = "-")) )
#length(NewpatId)==length(colnames(gistic_thresh))
#length(unique(NewpatId))==length(colnames(gistic_thresh)) OK

CNVg=gistic_thresh %>% t() 
colnames(CNVg)=GenesNames
rownames(CNVg)=NewpatId

#NOte the rows 1 2 3 are useless...only used for checks
CNVgtibble=CNVg %>% as_tibble(., rownames = "PatientID")
dim(CNVgtibble)

CNVgtibble1=CNVgtibble %>% filter(! PatientID %in% c("Gene.Symbol","Locus.ID", "Cytoband") ) %>% 
  mutate(across(-PatientID, as.numeric)) 


 

##select only patient with clinical info: length(IdPatientstoinclude1)=258
CNVgtibble2=CNVgtibble1 %>% filter(PatientID %in% IdPatientstoinclude1)

dim(CNVgtibble2) #257 patients, need to remove the patient in   CNV but no in clinical


#clinical data
myclinNONull10 = myclinNONull9 %>%  filter(cases.submitter_id %in%  CNVgtibble2$PatientID)

#check
myclinNONull10 %>% pull(cases.submitter_id) %>% unique() %>% length() #257 OK



###some preliminary check: 24777 genes

#select the variable with more variability (scaled)
#Entropy=apply(CNV_mat3,2, 
 
#myclinNONull11=myclinNONull10 %>% summarise(prognMALE=dplyr::first(ifelse(demographic.gender =="female",0,1)),
#                                            prognTUMORgrade3=dplyr::first(ifelse(diagnoses.tumor_grade=="G2",0,1)),
#                                            prognAGE=dplyr::first(as.numeric(demographic.age_at_index)),
#                                            Treat=1+dplyr::first(HasRadioAdj)*1,
#                                            Outcome=dplyr::first(ifelse(FUTrtOutcome=="Complete Response",2,
#                                                                        ifelse(FUTrtOutcome=="Stable Disease",1, 
#                                                                               ifelse(FUTrtOutcome=="Progressive Disease",0,NA))))
#)


#Create dataset of clinical data with only  one row per patient
invariant_cols= myclinNONull10 %>%
  group_by(cases.submitter_id) %>%
  summarise(across(everything(), ~ n_distinct(.x, na.rm = TRUE))) %>%
  # A column is invariant if for ALL patients it has 1 unique value
  summarise(across(everything(), ~ all(.x <= 1))) %>%
  ungroup() %>%
  # convert to vector of column names
  pivot_longer(everything(), names_to = "col", values_to = "keep") %>%
  filter(keep) %>%
  pull(col)

clean_myclin = myclinNONull10 %>%
  select(all_of(invariant_cols)) %>%
  distinct(cases.submitter_id, .keep_all = TRUE) %>% ungroup()

clean_myclin=clean_myclin %>% select(c("cases.submitter_id",                            
                          "cases.primary_site",                            
                          "demographic.age_at_index",                      
                          "demographic.ethnicity",                         
                          "demographic.gender",
                          "demographic.race",                              
                          "diagnoses.icd_10_code",                         
                          "diagnoses.laterality",                          
                          "diagnoses.morphology",                          
                          "diagnoses.sites_of_involvement",                
                          "diagnoses.supratentorial_localization",         
                          "diagnoses.tissue_or_organ_of_origin",           
                          "diagnoses.tumor_grade",                         
                          "diagnoses.year_of_diagnosis",
                          "FUTrtOutcome",                                  
                          "OnlyChemo",                                
                          "OnlyRadio",                                     
                          "HasRadioAdj",                                  
                          "HasChemo",                                      
                          "HasRadio" ) )





DatasetComplete=clean_myclin %>% rename(PatientID=cases.submitter_id) %>% 
  left_join(CNVgtibble2 %>%  #CNV data
              mutate(across(
                .cols = -PatientID,     # all columns except the ID
                .fns  = ~ .x + 2                 # shift values: -2→0, -1→1, 0→2, 1→3, 2→4
              )) )



#save dataset complete
LIST=list(Dataset=DatasetComplete, namesCNVvar=colnames(CNVgtibble2)[-1] )

base::saveRDS(LIST, file = "/home/lorenzo/Documents/UNIFI/PhD/Peoject/COMPLETELGGdata")










###END OF FILE
ScaledEntropy=function(x){
  p=table(x)
  p=p/sum(p)
  return((-sum(p*log(p,base = 2)))/log(length(p),2))
}

EntropyPerGene=CNVgtibble2 %>%
  summarise(across(-PatientID, ScaledEntropy))



plot(density(EntropyPerGene %>% as.numeric()))

#
Entropythesh=0.70
sum((EntropyPerGene %>% as.numeric())>Entropythesh)


vartoselect=(EntropyPerGene %>% as.numeric())>Entropythesh


dplot=CNVgtibble2[,-1][,vartoselect] %>% as.matrix()
rownames(dplot)=CNVgtibble2$PatientID#=mutation_count2[,1] %>% pull()


temptrt=myclinNONull10 %>% group_by(cases.submitter_id)%>%
  summarize(Trt=1+dplyr::first(HasRadioAdj*1)*1)
Trt=temptrt$Trt
names(Trt)=temptrt$cases.submitter_id

intersect(rownames(dplot), temptrt$cases.submitter_id)


tempresp=myclinNONull10%>% group_by(cases.submitter_id)%>%
  summarize(Resp=dplyr::first(FUTrtOutcome))
Resp=tempresp$Resp
names(Resp)=tempresp$cases.submitter_id


annot = data.frame(trt = Trt[rownames(dplot)],  resp= Resp[rownames(dplot)])
rownames(annot) = rownames(dplot)

pheatmap::pheatmap(dplot[,],annotation_row  = annot[,])

#trt 1
pheatmap::pheatmap(dplot[annot$trt==1,],annotation_row  = annot[annot$trt==1,])
#trt 2
pheatmap::pheatmap(dplot[annot$trt==2,],annotation_row  = annot[annot$trt==2,])

#resp 1
pheatmap::pheatmap(dplot[annot$resp=="Stable Disease",],annotation_row  = annot[annot$resp=="Stable Disease",])
#trt 2
x11()
pheatmap::pheatmap(dplot[annot$resp=="Complete Response",],annotation_row  = annot[annot$resp=="Complete Response",])

#trt 2
x11()
pheatmap::pheatmap(dplot[annot$resp=="Progressive Disease",],annotation_row  = annot[annot$resp=="Progressive Disease",])






#trt 1 reps 1 
pheatmap::pheatmap(dplot[annot$trt==1 & annot$resp=="Stable Disease",],annotation_row  = annot[annot$trt==1 & annot$resp=="Stable Disease",])
#trt 1 reps 2

pheatmap::pheatmap(dplot[annot$trt==2 & annot$resp=="Stable Disease",],annotation_row  = annot[annot$trt==2 & annot$resp=="Stable Disease",])




#trt 1 reps 1 
pheatmap::pheatmap(dplot[annot$trt==1 & annot$resp=="Complete Response",],annotation_row  = annot[annot$trt==1 & annot$resp=="Complete Response",])
#trt 1 reps 2

pheatmap::pheatmap(dplot[annot$trt==2 & annot$resp=="Complete Response",],annotation_row  = annot[annot$trt==2 & annot$resp=="Complete Response",])




Predictivevarnames=colnames(CNVgtibble2)[-1]


##join
myclinNONull11=myclinNONull10 %>% summarise(prognMALE=dplyr::first(ifelse(demographic.gender =="female",0,1)),
                                            prognTUMORgrade3=dplyr::first(ifelse(diagnoses.tumor_grade=="G2",0,1)),
                                            prognAGE=dplyr::first(as.numeric(demographic.age_at_index)),
                                            Treat=1+dplyr::first(HasRadioAdj)*1,
                                            Outcome=dplyr::first(ifelse(FUTrtOutcome=="Complete Response",2,
                                                                        ifelse(FUTrtOutcome=="Stable Disease",1, 
                                                                               ifelse(FUTrtOutcome=="Progressive Disease",0,NA))))
)






DatasetComplete


CLINICALVARtoExclude= colnames(DatasetComplete)  #setdiff( colnames(myclinNONull10),
                              c("PatientId","demographic.age_at_index", 
                                "demographic.gender","diagnoses.sites_of_involvement",
                                "diagnoses.tumor_grade","FUTrtOutcome", "HasRadioAdj",
                                "cases.submitter_id"
                            ))


NOTselected_genes70 =EntropyPerGene %>%
  pivot_longer(everything(), names_to = "Gene", values_to = "Entropy") %>%
  filter(Entropy <= 0.7) %>%
  pull(Gene)

NOTselected_genes65 =EntropyPerGene %>%
  pivot_longer(everything(), names_to = "Gene", values_to = "Entropy") %>%
  filter(Entropy <= 0.65) %>%
  pull(Gene)


Predictive70names=setdiff(Predictivevarnames, NOTselected_genes70)
Predictive65names=setdiff(Predictivevarnames, NOTselected_genes65)

VARtoExclude70=c(CLINICALVARtoExclude, NOTselected_genes70)
VARtoExclude65=c(CLINICALVARtoExclude, NOTselected_genes65)


DATA70=DatasetComplete %>% select(-NOTselected_genes70)

# %>% 
  mutate(prognAGE=as.numeric(demographic.age_at_index),
         prognTUMORgrade3=ifelse(diagnoses.tumor_grade=="G2",0,1),
         prognMALE=ifelse(diagnoses.tumor_grade=="female",0,1),
         outcome=ifelse(FUTrtOutcome=="Complete Response",2,
                        ifelse(FUTrtOutcome=="Stable Disease",1, 
                               ifelse(FUTrtOutcome=="Progressive Disease",0,NA)))
  )

DATA65=DatasetComplete %>% select(-NOTselected_genes65)
  mutate(prognAGE=as.numeric(demographic.age_at_index),
         prognTUMORgrade3=ifelse(diagnoses.tumor_grade=="G2",0,1),
         prognMALE=ifelse(diagnoses.tumor_grade=="female",0,1),
         outcome=ifelse(FUTrtOutcome=="Complete Response",2,
                        ifelse(FUTrtOutcome=="Stable Disease",1, 
                               ifelse(FUTrtOutcome=="Progressive Disease",0,NA)))
  )



#need to convert into matrix to be 

DATAtoSAVE=list(Entropy070=list(Data=DATA70, predictivenames=Predictive70names ),
     Entropy065=list(Data=DATA65, predictivenames=Predictive65names))





saveRDS( DATAtoSAVE,file = "~/Documents/UNIFI/PhD/Peoject/tcga/MYData2")











