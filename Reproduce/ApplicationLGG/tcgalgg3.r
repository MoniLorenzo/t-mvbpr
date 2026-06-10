#code 2/2 to build dataset from LGG data
#author lorenzo moni


rm(list = ls())
library(tidyverse)
#main dataset
Data=readRDS("/home/lorenzo/Documents/UNIFI/PhD/Peoject/COMPLETELGGdata")


ScaledEntropy=function(x){
  p=table(x)
  p=p/sum(p)
  return((-sum(p*log(p,base = 2)))/log(length(p),2))
}

#Datasets:
#1. block of 50 genes (only genes with names and chromosome are considered): in each block the gene with max entropy is selected  
#2.a/b only genes with (a)entropy>0.5 and (b)entropy>0.65 are considered,
# then 10 dataset where we pick at rnd genes for a final ngenes ~500 


##########
#Clean data 
set.seed(10)
Clinical=Data$Dataset[,1:20] %>% dplyr::select(PatientID, 
                                               FUTrtOutcome,
                                               HasRadioAdj,
                                               demographic.age_at_index, 
                                               diagnoses.morphology,
                                               demographic.gender,
                                               #diagnoses.year_of_diagnosis
                                               diagnoses.laterality,diagnoses.tumor_grade
                                               ) %>% 
  rename(Outcome=FUTrtOutcome, Treatment=HasRadioAdj, 
         #YYPD=diagnoses.year_of_diagnosis,
         Tumorgrade=diagnoses.tumor_grade,
         Age=demographic.age_at_index, 
         Sex=demographic.gender,
         Morphology=diagnoses.morphology,
         Laterality=diagnoses.laterality) %>% 
  mutate(Age=as.numeric(Age), Morphology=as.factor(Morphology),
       #  YYPD=as.numeric(YYPD),
         Sex=as.factor(Sex)) %>% 
  mutate(
    Laterality = if_else(
      Laterality %in% c("Midline","'--"),
      sample(c("Left", "Right"), size = n(), replace = TRUE),
      Laterality),
      Tumorgrade = if_else(
        Tumorgrade %in% c("'--"),
        sample(c("G2", "G3"), size = n(), replace = TRUE),
        Tumorgrade)) %>% 
  mutate(Laterality=as.factor(Laterality),Tumorgrade=as.factor(Tumorgrade),
         Outcome=as.factor(Outcome),             
         Treatment=Treatment*1+1)



#---------------------------------------------------------------------
# 1.
library(biomaRt)

genes =Data$namesCNVvar[]
length(genes)


mart = useMart(
  "ENSEMBL_MART_ENSEMBL",
  dataset = "hsapiens_gene_ensembl",
  host = "https://grch37.ensembl.org"   
)

positions = getBM(
  attributes = c("external_gene_name","chromosome_name","start_position","end_position","strand"),
  filters = "external_gene_name",
  values = genes,
  mart = mart
)

#keep only records with   chromosome_name  c(1:22, "X", "Y", "MT")
positions_clean <- positions %>%
  dplyr::filter(chromosome_name %in% c(1:22, "X", "Y", "MT")) %>%   # keep real chromosomes
  dplyr::group_by(external_gene_name) %>%
  dplyr::filter(n() == 1) %>%                                      # keep only if 1 record
  dplyr::ungroup() %>% 
  filter(external_gene_name %in% genes)                           # keep only if genes is actually in original genes vec


length(positions_clean$external_gene_name)

#genes witn no position info
genesNOpos=setdiff(genes,positions_clean$external_gene_name)
length(genesNOpos)

#24776==1422+23354
#setdiff(positions_clean$external_gene_name, genes) #check


#create block of 50 genes for each chromosome 
positions_clean1=positions_clean %>%
  group_by(chromosome_name) %>%
  arrange(start_position) %>%
  mutate(block = ceiling(row_number() / 50)) %>% ungroup() %>% 
  group_by(chromosome_name, block)


GenesMaxentropy=positions_clean1 %>%
  group_by(chromosome_name, block) %>%
  summarise(namesblocklist = list(external_gene_name), .groups = "drop") %>% 
  mutate(
    cnv_matrix = lapply(
      namesblocklist,
      function(g, Xcnv=Data$Dataset) Xcnv[,g,drop=FALSE]
    )
  ) %>%   mutate(
    cnv_matrix = map(namesblocklist, ~ Data$Dataset[, .x, drop = FALSE]),
    
    Maxentropygenes = map(cnv_matrix, ~ {
      .x %>%
        summarise(across(everything(), ScaledEntropy)) %>%
        pivot_longer(everything(), names_to = "gene", values_to = "entropy") %>%
        filter(entropy == max(entropy)) %>%
        pull(gene)
    })
  ) %>% 
  mutate(
    MaxEntropyGenes_identical =
      map2_lgl(cnv_matrix, Maxentropygenes, ~ {
        submat <- .x[, .y, drop = FALSE]
        
        # if single gene → TRUE
        if (ncol(submat) <= 1) return(TRUE)
        
        # your elegant summarise() + across() + pick()
        submat %>%
          summarise(
            all_identical = all(across(everything(), ~ all(. == pick(1)[[1]])))
          ) %>%
          pull(all_identical)
      })
  ) %>% mutate(first_MaxEntropyGene = map_chr(Maxentropygenes, ~ .x[1]))

#extract names of genes with max entropy. Note that if MaxEntropyGenes_identical=TRUE
#all the genes in the block are identical
ExtractedGenesMaxent=GenesMaxentropy %>% pull(first_MaxEntropyGene)

#heatmap
Data$Dataset[,ExtractedGenesMaxent] %>% as.matrix() %>% pheatmap::pheatmap(.)



#create dataset
DF1=list(data=data.frame(Clinical, (Data$Dataset[,ExtractedGenesMaxent]-2) %>%   
            mutate(across(everything(), 
                          ~ factor(as.character(.x), 
                                   levels = c("-2","-1","0","1","2"))))),
         namesCNV=ExtractedGenesMaxent)

str(DF1)
saveRDS(object = DF1,"~/Documents/UNIFI/PhD/Peoject/tcga/DATAlgg1")

###############################################################
#2.   
entropybygenes=Data$Dataset[Data$namesCNVvar] %>%
    summarise(across(everything(), ScaledEntropy)) %>%
    pivot_longer(everything(), names_to = "gene", values_to = "entropy") 


Geneswithentrgre50 =entropybygenes%>% filter(entropy >0.5)  %>% arrange(entropy)
Geneswithentrgre65 =entropybygenes%>% filter(entropy >=0.65)  %>% arrange(entropy)

#50: 6758 genes, 65: 1,735

#create 10 database 
Inda=rep(FALSE,nrow(Geneswithentrgre50))
Indb=rep(FALSE,nrow(Geneswithentrgre65))
Inda[1:500]=Indb[1:500]=TRUE


DataTwo=list(TwoA=list(), TwoB=list())

for (se in 1:10) {
  set.seed(1154+se)
  
  DataTwo$TwoA[[paste0(se)]]=Data$Dataset[,Geneswithentrgre50$gene[sample(Inda)]]
  DataTwo$TwoB[[paste0(se)]]=Data$Dataset[,Geneswithentrgre65$gene[sample(Indb)]]
}


#############









